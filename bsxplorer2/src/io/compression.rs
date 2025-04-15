#[cfg(feature = "compression")]
mod inner {
    use polars::io::mmap::MmapBytesReader;
    use std::io::Write;
    use std::{fs::File, io::Read};

    pub enum Compression {
        None,
        Gz,
        Zstd,
        Lz4,
        Xz2,
        Bzip2,
        Zip
    }

    impl Compression {
        pub fn name(&self) -> &str {
            match self {
                Compression::None => "none",
                Compression::Gz => "gzip",
                Compression::Zstd => "zstd",
                Compression::Lz4 => "lz4",
                Compression::Xz2 => "xz2",
                Compression::Bzip2 => "bzip2",
                Compression::Zip => "zip",
            }
        }

        pub fn get_decoder(&self, handle: File) -> Box<dyn MmapBytesReader> {
            match self {
                Compression::Gz => Box::new(flate2::read::GzDecoder::new(handle)),
                Compression::Zstd => Box::new(zstd::Decoder::new(handle)),
                Compression::Lz4 => Box::new(lz4::Decoder::new(handle)),
                Compression::Xz2 => Box::new(xz2::read::XzDecoder::new(handle)),
                Compression::Bzip2 => Box::new(bzip2::read::BzDecoder::new(handle)),
                Compression::Zip => Box::new(zip::read::ZipArchive::new(handle)),
                Compression::None => Box::new(handle),
            }
        }

        pub fn get_encoder(
            &self,
            handle: File,
            compression_level: u32
        ) -> anyhow::Result<Box<dyn Write>> {
            let encoder: Box<dyn Write> = match self {
                Compression::Gz => Box::new(flate2::write::GzEncoder::new(
                    handle,
                    flate2::Compression(compression_level)
                ))
                ,
                Compression::Zstd => Box::new(zstd::Encoder::new(
                    handle,
                    compression_level as i32
                )?),
                Compression::Lz4 => {
                    let encoder = lz4::EncoderBuilder::new()
                        .level(compression_level)
                        .build(handle)?;
                    Box::new(encoder)
                },
                Compression::Xz2 => Box::new(xz2::write::XzEncoder::new(
                    handle,
                    compression_level
                )),
                Compression::Bzip2 => Box::new(bzip2::write::BzEncoder::new(
                    handle, bzip2::Compression(compression_level)
                )),
                Compression::Zip => Box::new(zip::write::ZipWriter::new(handle)),
                Compression::None => Box::new(handle),
            };
            Ok(encoder)
        }
    }

}

#[cfg(feature = "compression")]
pub use inner::*;
