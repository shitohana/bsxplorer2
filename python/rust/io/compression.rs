use bsxplorer2::io::compression::Compression;
use pyo3::prelude::*; // Assuming bsxplorer2 is the crate name

#[pyclass(name = "Compression", eq, eq_int)]
#[derive(Clone, Debug, Ord, PartialOrd, PartialEq, Eq, Hash)]
pub enum PyCompression {
    // Had to name it no instead of none, because python interprets
    // None incorrectly
    No,
    Gz,
    Zstd,
    Lz4,
    Xz2,
    Bzip2,
    Zip,
}

impl From<PyCompression> for Compression {
    fn from(py: PyCompression) -> Self {
        match py {
            PyCompression::No => Compression::None,
            PyCompression::Gz => Compression::Gz,
            PyCompression::Zstd => Compression::Zstd,
            PyCompression::Lz4 => Compression::Lz4,
            PyCompression::Xz2 => Compression::Xz2,
            PyCompression::Bzip2 => Compression::Bzip2,
            PyCompression::Zip => Compression::Zip,
        }
    }
}

impl From<Compression> for PyCompression {
    fn from(compression: Compression) -> Self {
        match compression {
            Compression::None => PyCompression::No,
            Compression::Gz => PyCompression::Gz,
            Compression::Zstd => PyCompression::Zstd,
            Compression::Lz4 => PyCompression::Lz4,
            Compression::Xz2 => PyCompression::Xz2,
            Compression::Bzip2 => PyCompression::Bzip2,
            Compression::Zip => PyCompression::Zip,
        }
    }
}
