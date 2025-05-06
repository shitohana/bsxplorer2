use std::ops::{Add, Sub};

use bsxplorer2::data_structs::coords::{Contig, GenomicPosition};
use bsxplorer2::data_structs::enums::Strand as RsStrand;
use pyo3::exceptions::{PyNotImplementedError, PyValueError};
use pyo3::prelude::*;

use super::utils::Strand as PyStrand;


#[pyclass(get_all, set_all)]
#[derive(Debug, Clone)]
pub struct PyGenomicPosition {
    seqname:  String,
    position: u32,
}

#[pymethods]
impl PyGenomicPosition {
    #[new]
    fn new(
        seqname: String,
        position: u32,
    ) -> Self {
        Self { seqname, position }
    }

    // Add method for display
    fn __str__(&self) -> String {
        // Leverage the existing Display implementation
        format!(
            "{}",
            GenomicPosition::new(self.seqname.clone(), self.position)
        )
    }

    fn __repr__(&self) -> String {
        format!(
            "PyGenomicPosition(seqname='{}', position={})",
            self.seqname, self.position
        )
    }

    // Add comparison methods. PartialOrd returns Option<Ordering>.
    // If None, it means seqnames are different, which the Rust Ord
    // implementation panics on. In Python, we can raise a ValueError or
    // NotImplemented. Let's raise ValueError for clarity for ordering ops.
    // Equality/Inequality is allowed between different seqnames.
    fn __richcmp__(
        &self,
        other: &PyGenomicPosition,
        op: pyo3::basic::CompareOp,
    ) -> PyResult<bool> {
        // Convert to Rust type for comparison
        let self_rust: GenomicPosition<String, u32> = self.into();
        let other_rust: GenomicPosition<String, u32> = other.into();


        match op {
            pyo3::basic::CompareOp::Eq => Ok(self_rust == other_rust),
            pyo3::basic::CompareOp::Ne => Ok(self_rust != other_rust),
            _ => {
                // Use the Rust PartialOrd implementation for ordering
                match self_rust.partial_cmp(&other_rust) {
                    Some(ordering) => Ok(op.matches(ordering).into()),
                    None => {
                        // For ordering ops (<, <=, >, >=), different seqnames
                        // are not comparable
                        Err(PyValueError::new_err(
                            "Cannot order genomic positions on different \
                             chromosomes/seqnames",
                        ))
                    },
                }
            },
        }
    }

    // Add arithmetic methods (__add__, __sub__)
    fn __add__(
        &self,
        other: &PyGenomicPosition,
    ) -> PyResult<Option<Self>> {
        // Convert to Rust type for operation
        let self_rust: GenomicPosition<String, u32> = self.into();
        let other_rust: GenomicPosition<String, u32> = other.into();

        // Use the Rust Add implementation
        match self_rust.add(other_rust) {
            Some(result_rust) => Ok(Some(result_rust.into())), /* Convert result back to PyO3 type */
            None => Ok(None), // Different seqnames
        }
    }

    fn __sub__(
        &self,
        other: &PyGenomicPosition,
    ) -> PyResult<Option<Self>> {
        // Convert to Rust type for operation
        let self_rust: GenomicPosition<String, u32> = self.into();
        let other_rust: GenomicPosition<String, u32> = other.into();

        // Use the Rust Sub implementation
        match self_rust.sub(other_rust) {
            Some(result_rust) => Ok(Some(result_rust.into())), /* Convert result back to PyO3 type */
            None => Ok(None), // Different seqnames or rhs > lhs
        }
    }
}

// Implement conversions between Rust and PyO3 structs
impl From<GenomicPosition<String, u32>> for PyGenomicPosition {
    fn from(gp: GenomicPosition<String, u32>) -> Self {
        Self {
            seqname:  gp.seqname(),
            position: gp.position(),
        }
    }
}

impl From<&PyGenomicPosition> for GenomicPosition<String, u32> {
    fn from(py_gp: &PyGenomicPosition) -> Self {
        Self::new(py_gp.seqname.clone(), py_gp.position)
    }
}


#[pyclass(get_all, set_all)] // Allows access to seqname, start, end directly
#[derive(Debug, Clone)] // Need these for conversion
pub struct PyContig {
    seqname: String,
    start:   u32,
    end:     u32,
    strand:  PyStrand, // Store the Rust enum internally
}

#[pymethods]
impl PyContig {
    #[new]
    #[pyo3(signature = (
        seqname,
        start,
        end,
        strand = "."
    ))] // Default strand to "None" if not provided
    fn new(
        seqname: String,
        start: u32,
        end: u32,
        strand: &str,
    ) -> PyResult<Self> {
        if start > end {
            return Err(PyValueError::new_err(
                "Start position must be less than or equal to end position",
            ));
        }
        let strand = match strand.to_lowercase().as_str() {
            "+" | "forward" => PyStrand::Forward,
            "-" | "reverse" => PyStrand::Reverse,
            "." | "none" => PyStrand::None,
            _ => {
                return Err(PyValueError::new_err(format!(
                    "Invalid strand value: '{}'. Must be '+', '-', '.', \
                     'Forward', 'Reverse', or 'None' (case-insensitive).",
                    strand
                )))
            },
        };
        Ok(Self {
            seqname,
            start,
            end,
            strand,
        })
    }

    // Expose strand as a string
    #[getter]
    fn strand_str(&self) -> String {
        match self.strand {
            PyStrand::Forward => "+".to_string(),
            PyStrand::Reverse => "-".to_string(),
            PyStrand::None => ".".to_string(),
        }
    }

    // Length method
    fn length(&self) -> u32 { self.end - self.start }

    // GenomicPosition methods
    fn start_gpos(&self) -> PyGenomicPosition {
        GenomicPosition::new(self.seqname.clone(), self.start).into()
    }

    fn end_gpos(&self) -> PyGenomicPosition {
        GenomicPosition::new(self.seqname.clone(), self.end).into()
    }

    // Extend methods
    fn extend_upstream(
        &mut self,
        length: u32,
    ) {
        self.start = self.start.saturating_sub(length);
    }

    fn extend_downstream(
        &mut self,
        length: u32,
    ) {
        self.end = self.end.saturating_add(length);
    }

    // is_in method
    fn is_in(
        &self,
        other: &PyContig,
    ) -> bool {
        // Convert to Rust type for the method call
        let self_rust: Contig<String, u32> = self.into();
        let other_rust: Contig<String, u32> = other.into();
        self_rust.is_in(&other_rust)
    }

    // Add method for display
    fn __str__(&self) -> String {
        // Leverage the existing Display implementation
        format!(
            "{}",
            Contig::new(
                self.seqname.clone(),
                self.start,
                self.end,
                RsStrand::from(self.strand)
            )
        )
    }

    fn __repr__(&self) -> String {
        format!(
            "PyContig(seqname='{}', start={}, end={}, strand='{}')",
            self.seqname,
            self.start,
            self.end,
            self.strand_str()
        )
    }

    // Add comparison methods. PartialOrd returns Option<Ordering> and returns
    // None if seqnames differ or regions intersect. For PartialEq, it
    // checks equality including strand.
    fn __richcmp__(
        &self,
        other: &PyContig,
        op: pyo3::basic::CompareOp,
    ) -> PyResult<bool> {
        // Convert to Rust type for comparison
        let self_rust: Contig<String, u32> = self.into();
        let other_rust: Contig<String, u32> = other.into();

        match op {
            pyo3::basic::CompareOp::Eq => Ok(self_rust == other_rust),
            pyo3::basic::CompareOp::Ne => Ok(self_rust != other_rust),
            _ => {
                // Use the Rust PartialOrd implementation for ordering
                match self_rust.partial_cmp(&other_rust) {
                    Some(ordering) => Ok(op.matches(ordering)),
                    None => {
                        // Rust PartialOrd returns None for different seqnames
                        // OR intersecting regions.
                        // Python comparison should ideally return False or
                        // NotImplemented for incomparable items.
                        // Returning NotImplemented allows chaining comparisons.
                        Err(PyNotImplementedError::new_err(
                            "Contigs on different seqnames or with \
                             intersecting regions are not strictly comparable.",
                        ))
                    },
                }
            },
        }
    }
}

// Implement conversions between Rust and PyO3 structs
impl From<Contig<String, u32>> for PyContig {
    fn from(c: Contig<String, u32>) -> Self {
        Self {
            seqname: c.seqname(),
            start:   c.start(),
            end:     c.end(),
            strand:  PyStrand::from(c.strand()),
        }
    }
}

impl From<&PyContig> for Contig<String, u32> {
    fn from(py_c: &PyContig) -> Self {
        Self::new(
            py_c.seqname.clone(),
            py_c.start,
            py_c.end,
            RsStrand::from(py_c.strand),
        )
    }
}


// Function to add these classes to the module
pub fn add_coords_module(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyGenomicPosition>()?;
    m.add_class::<PyContig>()?;
    Ok(())
}
