use std::{io::Read, sync::Arc};
use std::io::{self, BufRead};
use arrow::array::RecordBatch;
use arrow::error::ArrowError;
use xz2::bufread::XzDecoder;
use zstd::stream::Decoder as ZstdDecoder;
use bzip2::bufread::BzDecoder;
use flate2::bufread::GzDecoder;
use tar::Archive;


pub enum ReportType {
    BEDGRAPH,
    COVERAGE,
    CGMAP,
    BISMARK,
}

pub enum TarBasedFmt {
    TAR, 
    BZ2, 
    GZ, 
    XZ, 
    ZSTD, 
    /// If file is not compressed
    NONE
}

/// Create decoder for specified compression format
pub fn create_tar_decoder(
    dat: impl BufRead + 'static,
    fmt: TarBasedFmt,
) -> io::Result<Archive<Box<dyn Read>>> {
    let r: Box<dyn Read> = match fmt {
        TarBasedFmt::TAR => Box::new(dat),
        TarBasedFmt::BZ2 => Box::new(BzDecoder::new(dat)),
        TarBasedFmt::GZ => Box::new(GzDecoder::new(dat)),
        TarBasedFmt::XZ => Box::new(XzDecoder::new(dat)),
        TarBasedFmt::ZSTD => Box::new(ZstdDecoder::with_buffer(dat)?),
        TarBasedFmt::NONE => Box::new(dat)
    };

    Ok(Archive::new(r))
}


impl ReportType {
    pub fn schema(&self) -> arrow::datatypes::Schema {
        use { arrow::datatypes::Field, arrow::datatypes::DataType };
        match self {
            ReportType::BEDGRAPH => {
                arrow::datatypes::Schema::new(vec![
                    Field::new("chr", DataType::Utf8, false),
                    Field::new("start", DataType::UInt64, false),
                    Field::new("end", DataType::UInt64, false),
                    Field::new("density", DataType::Float32, false),
                ])
            },
            ReportType::COVERAGE => {
                arrow::datatypes::Schema::new(vec![
                    Field::new("chr", DataType::Utf8, false),
                    Field::new("start", DataType::UInt64, false),
                    Field::new("end", DataType::UInt64, false),
                    Field::new("density", DataType::Float32, false),
                    Field::new("count_m", DataType::UInt32, false),
                    Field::new("count_um", DataType::UInt32, false),
                ])
            },
            ReportType::CGMAP => {
                arrow::datatypes::Schema::new(vec![
                    Field::new("chr", DataType::Utf8, false),
                    Field::new("nuc", DataType::Utf8, false),
                    Field::new("position", DataType::UInt64, false),
                    Field::new("context", DataType::Utf8, false),
                    Field::new("dinuc", DataType::Utf8, false),
                    Field::new("count_m", DataType::UInt32, false),
                    Field::new("count_total", DataType::UInt32, false),
                ])
            },
            ReportType::BISMARK => arrow::datatypes::Schema::new(vec![
                Field::new("chr", DataType::Utf8, false),
                Field::new("position", DataType::UInt64, false),
                Field::new("strand", DataType::Utf8, false),
                Field::new("count_m", DataType::UInt32, false),
                Field::new("count_um", DataType::UInt32, false),
                Field::new("context", DataType::Utf8, false),
                Field::new("trinuc", DataType::Utf8, false),
            ]),
        }
    }
    pub fn decoder(&self, batch_size: usize) -> arrow::csv::reader::Decoder {
        match self {
            ReportType::BISMARK | ReportType::CGMAP | ReportType::COVERAGE | ReportType::BEDGRAPH => {
                arrow::csv::ReaderBuilder::new(Arc::new(self.schema()))
                    .with_delimiter(b'\t')
                    .with_comment(b'#')
                    .with_header(false)
                    .with_batch_size(batch_size)
                    .build_decoder()
            }
        }
    }
    pub fn get_iterator(
        self, path: &str, batch_size: usize, compression: TarBasedFmt
    ) -> Result<impl Iterator<Item = Result<RecordBatch, ArrowError>>, ArrowError> {
        let mut reader = {
            let file = std::fs::File::open(path).unwrap();
            // let buf_reader = io::BufReader::with_capacity(1_000_000_000, file);
            // let tar_decoder = create_tar_decoder(buf_reader, compression)
            //     .unwrap()
            //     .into_inner();
            // io::BufReader::new(tar_decoder)
            io::BufReader::with_capacity(2 << 25, file)
        };
        let mut decoder = self.decoder(batch_size);
        let mut next = move || {
            loop {
                let buf = reader.fill_buf()?;
                let decoded = decoder.decode(buf)?;
                if decoded == 0 {
                    break;
                }
                reader.consume(decoded);
            }
            decoder.flush()
        };
        Ok(std::iter::from_fn(move || next().transpose()))
    }
}