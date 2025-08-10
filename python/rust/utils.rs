use polars::export::rayon::prelude::*;
use pyo3::exceptions::{PyIOError, PyValueError};
use pyo3::prelude::*;
use pyo3_file::PyFileLikeObject;
use std::fs::File;
use std::io::{Read, Seek, Write};
use std::os::fd::AsRawFd;

pub trait ReadHandle: Read + Seek + AsRawFd {}
impl<T: Read + Seek + AsRawFd> ReadHandle for T {}

pub trait SinkHandle: Write + Seek + 'static {}
impl<T: Write + Seek + 'static> SinkHandle for T {}

#[derive(Debug)]
pub enum FileOrFileLike {
    File(String),
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
        if let Ok(string) = ob.extract::<String>() {
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
    positions: Vec<Vec<f64>>,
    densities: Vec<Vec<f64>>,
) -> PyResult<(Vec<f64>, Vec<f64>)> {
    if positions.len() != densities.len() {
        return Err(PyValueError::new_err(
            "Positions and densities arrays lengths differ",
        ));
    }
    let zipped_positions = positions.concat();
    let zipped_densities = densities.concat();
    let mut zipped_points = zipped_positions
        .into_iter()
        .zip(zipped_densities.into_iter())
        .collect::<Vec<_>>();
    zipped_points.par_sort_unstable_by(|(p1, _d1), (p2, _d2)| {
        p1.partial_cmp(p2).expect("Unexpected NaN")
    });

    Ok(itertools::multiunzip(zipped_points.into_iter()))
}
