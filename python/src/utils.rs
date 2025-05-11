use std::fs::File;
use std::io::{Read, Seek, Write};

use pyo3::exceptions::PyIOError;
use pyo3::prelude::*;
use pyo3_file::PyFileLikeObject;

pub trait ReadHandle: Read + Seek {}
impl<T: Read + Seek> ReadHandle for T {}

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
        if let Ok(f) = PyFileLikeObject::py_with_requirements(
            ob.clone(),
            false,
            true,
            true,
            false,
        ) {
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
