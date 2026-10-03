use std::fs::File;
use std::io::{
    Read,
    Seek,
    Write,
};
use std::os::fd::{
    AsRawFd,
    RawFd,
};
use std::path::PathBuf;

use polars::export::rayon::prelude::*;
use pyo3::exceptions::{
    PyIOError,
    PyValueError,
};
use pyo3::prelude::*;
use pyo3_file::PyFileLikeObject;

pub trait ReadHandle: Read + Seek + AsRawFd {}
impl<T: Read + Seek + AsRawFd> ReadHandle for T {}

pub trait SinkHandle: Write + Seek + 'static {}
impl<T: Write + Seek + 'static> SinkHandle for T {}

#[derive(Debug)]
pub enum FileOrFileLike {
    File(PathBuf),
    ROnlyFileLike(PyFileLikeObject),
    RWFileLike(PyFileLikeObject),
}

impl FileOrFileLike {
    pub fn get_reader(self) -> PyResult<Box<dyn ReadHandle + 'static>> {
        Ok(match self {
            FileOrFileLike::File(path) => Box::new(File::open(path)?),
            FileOrFileLike::ROnlyFileLike(handle) => Box::new(handle),
            FileOrFileLike::RWFileLike(handle) => Box::new(handle),
        })
    }

    pub fn get_writer(self) -> PyResult<Box<dyn SinkHandle>> {
        Ok(match self {
            FileOrFileLike::File(path) => Box::new(File::create(path)?),
            FileOrFileLike::RWFileLike(handle) => Box::new(handle),
            FileOrFileLike::ROnlyFileLike(handle) => {
                return Err(PyIOError::new_err("File does not support writing"))
            },
        })
    }
}

impl<'py> FromPyObject<'py> for FileOrFileLike {
    fn extract_bound(ob: &Bound<'py, PyAny>) -> PyResult<Self> {
        // is a path
        if let Ok(string) = ob.extract::<PathBuf>() {
            return Ok(FileOrFileLike::File(string));
        }

        // is a file-like
        if let Ok(f) =
            PyFileLikeObject::py_with_requirements(ob.clone(), false, true, true, false)
        {
            return Ok(Self::RWFileLike(f));
        };

        let f = PyFileLikeObject::py_with_requirements(
            ob.clone(),
            true,
            false,
            true,
            false,
        )?;
        Ok(FileOrFileLike::ROnlyFileLike(f))
    }
}

#[pyfunction]
pub fn merge_metagene_values(
    py: Python<'_>,
    positions: Vec<Vec<f64>>,
    densities: Vec<Vec<f64>>,
) -> PyResult<(Vec<f64>, Vec<f64>)> {
    if positions.len() != densities.len() {
        return Err(PyValueError::new_err(
            "Positions and densities arrays lengths differ",
        ));
    }
    if positions
        .iter()
        .zip(&densities)
        .any(|(p, d)| p.len() != d.len())
    {
        return Err(PyValueError::new_err(
            "Each positions/densities pair must have equal lengths",
        ));
    }
    if positions.iter().flatten().any(|p| !p.is_finite()) {
        return Err(PyValueError::new_err("Positions must be finite"));
    }
    let zipped_positions = positions.concat();
    let zipped_densities = densities.concat();
    let mut zipped_points = zipped_positions
        .into_iter()
        .zip(zipped_densities.into_iter())
        .collect::<Vec<_>>();
    Ok(py.allow_threads(move || {
        zipped_points.par_sort_unstable_by(|(p1, _), (p2, _)| p1.total_cmp(p2));
        itertools::multiunzip(zipped_points.into_iter())
    }))
}

// Validate Python descriptors before entering the infallible AsRawFd interface.
// In particular, BytesIO.fileno() raises UnsupportedOperation.
pub struct MmapDescriptor {
    fd:     RawFd,
    _owner: Py<PyAny>,
}
impl AsRawFd for MmapDescriptor {
    fn as_raw_fd(&self) -> RawFd {
        self.fd
    }
}

pub enum MmapSource {
    Path(PathBuf),
    Descriptor(MmapDescriptor),
}
impl<'py> FromPyObject<'py> for MmapSource {
    fn extract_bound(ob: &Bound<'py, PyAny>) -> PyResult<Self> {
        if let Ok(path) = ob.extract::<PathBuf>() {
            return Ok(Self::Path(path));
        }
        let fd = ob.call_method0("fileno")?.extract::<RawFd>()?;
        if fd < 0 {
            return Err(PyValueError::new_err("File descriptor must be nonnegative"));
        }
        Ok(Self::Descriptor(MmapDescriptor {
            fd,
            _owner: ob.clone().unbind(),
        }))
    }
}
impl MmapSource {
    pub fn open(self) -> PyResult<bsxplorer2::io::bsx::BsxFileReader> {
        Ok(match self {
            Self::Path(path) => {
                bsxplorer2::io::bsx::BsxFileReader::try_new(File::open(path)?)?
            },
            Self::Descriptor(handle) => {
                bsxplorer2::io::bsx::BsxFileReader::try_new(handle)?
            },
        })
    }
}
