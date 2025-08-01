#[cfg(feature = "compression")]
mod inner {
    use std::convert::Infallible;
    use std::fmt::Display;
    use std::fs::File;
    use std::io::{
        copy,
        Seek,
        SeekFrom,
        Write,
    };
    use std::str::FromStr;

    // Added copy, Seek,
    // SeekFrom
    use polars::io::mmap::MmapBytesReader;
    use tempfile::tempfile; // Added tempfile

    #[derive(Clone, Debug, Copy, PartialEq, Eq)]
    pub enum Compression {
        None,
        Gz,
        Zstd,
        Lz4,
        Xz2,
        Bzip2,
        Zip,
    }

    impl FromStr for Compression {
        type Err = Infallible;

        fn from_str(s: &str) -> Result<Self, Self::Err> {
            match s {
                "none" => Ok(Compression::None),
                "gz" => Ok(Compression::Gz),
                "zstd" => Ok(Compression::Zstd),
                "lz4" => Ok(Compression::Lz4),
                "xz2" => Ok(Compression::Xz2),
                "bzip2" => Ok(Compression::Bzip2),
                "zip" => Ok(Compression::Zip),
                _ => unimplemented!(),
            }
        }
    }

    impl Display for Compression {
        fn fmt(
            &self,
            f: &mut std::fmt::Formatter<'_>,
        ) -> std::fmt::Result {
            let str = match self {
                Compression::None => String::from("none"),
                Compression::Gz => String::from("gz"),
                Compression::Zstd => String::from("zstd"),
                Compression::Lz4 => String::from("lz4"),
                Compression::Xz2 => String::from("xz2"),
                Compression::Bzip2 => String::from("bzip2"),
                Compression::Zip => String::from("zip"),
            };
            write!(f, "{}", str)
        }
    }

    impl Compression {
        /// Returns the name of the compression algorithm.
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

        /// Returns a `MmapBytesReader` wrapped in a decompressor.
        pub fn get_decoder(
            &self,
            handle: File, // Added mut
        ) -> anyhow::Result<Box<dyn MmapBytesReader>> {
            // Changed return type
            let mut temp_file = tempfile()?; // Create a temp file

            match self {
                Compression::Gz => {
                    let mut decoder = flate2::read::GzDecoder::new(handle);
                    copy(&mut decoder, &mut temp_file)?;
                },
                Compression::Zstd => {
                    // zstd::Decoder::new itself returns a Result
                    let mut decoder = zstd::Decoder::new(handle)?;
                    copy(&mut decoder, &mut temp_file)?;
                },
                Compression::Lz4 => {
                    // lz4::Decoder::new itself returns a Result
                    let mut decoder = lz4::Decoder::new(handle)?;
                    copy(&mut decoder, &mut temp_file)?;
                },
                Compression::Xz2 => {
                    let mut decoder = xz2::read::XzDecoder::new(handle);
                    copy(&mut decoder, &mut temp_file)?;
                },
                Compression::Bzip2 => {
                    let mut decoder = bzip2::read::BzDecoder::new(handle);
                    copy(&mut decoder, &mut temp_file)?;
                },
                Compression::Zip => {
                    // zip::ZipArchive::new itself returns a Result
                    let mut archive = zip::ZipArchive::new(handle)?;
                    if !archive.is_empty() {
                        // Extract the first file
                        let mut file_in_zip = archive.by_index(0)?;
                        copy(&mut file_in_zip, &mut temp_file)?;
                    }
                    else {
                        // Handle empty zip file - return empty temp file
                    }
                },
                Compression::None => {
                    // Directly copy if no compression
                    return Ok(Box::new(handle)); // Return the original handle
                },
            }

            // Rewind the temporary file to the beginning before returning
            temp_file.seek(SeekFrom::Start(0))?;

            Ok(Box::new(temp_file)) // Return the temp file handle
        }

        /// Returns a `Write` wrapped in a compressor.
        pub fn get_encoder<W: Write + Seek + 'static>(
            &self,
            handle: W,
            compression_level: u32,
        ) -> anyhow::Result<Box<dyn Write>> {
            let encoder: Box<dyn Write> = match self {
                Compression::Gz => {
                    Box::new(flate2::write::GzEncoder::new(
                        handle,
                        flate2::Compression::new(compression_level),
                    ))
                },
                Compression::Zstd => {
                    Box::new(zstd::Encoder::new(handle, compression_level as i32)?)
                },
                Compression::Lz4 => {
                    let encoder = lz4::EncoderBuilder::new()
                        .level(compression_level)
                        .build(handle)?;
                    Box::new(encoder)
                },
                Compression::Xz2 => {
                    Box::new(xz2::write::XzEncoder::new(handle, compression_level))
                },
                Compression::Bzip2 => {
                    Box::new(bzip2::write::BzEncoder::new(
                        handle,
                        bzip2::Compression::new(compression_level),
                    ))
                },
                Compression::Zip => Box::new(zip::write::ZipWriter::new(handle)),
                Compression::None => Box::new(handle),
            };
            Ok(encoder)
        }
    }
}

#[cfg(feature = "compression")]
pub use inner::*;
