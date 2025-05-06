use bsxplorer2::io::compression::Compression;
use pyo3::prelude::*; // Assuming bsxplorer2 is the crate name

#[pyclass(name = "Compression")]
#[derive(Clone)]
pub struct PyCompression(pub Compression);

impl From<Compression> for PyCompression {
    fn from(compression: Compression) -> Self { PyCompression(compression) }
}
impl From<PyCompression> for Compression {
    fn from(py_compression: PyCompression) -> Self { py_compression.0 }
}
