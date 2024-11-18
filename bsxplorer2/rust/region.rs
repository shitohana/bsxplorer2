use bsx_rs::region::{
    FillNullStrategy as FillNullStrategyRust, NullHandleStrategy as NullHandleStrategyRust,
    RegionCoordinates as RegionCoordinatesRust, RegionData as RegionDataRust,
};
use polars::prelude::*;
/// Python port of bsx_rs::region. Implemented objects:
/// -   [RegionCoordinates] ([bsx_rs::region::RegionCoordinates])
///     - Initializer
///     - Getters for region coordinates
/// -   [FillNullStrategy] ([bsx_rs::region::FillNullStrategy])
///     - Initializer
/// -   [NullHandleStrategy] ([bsx_rs::region::NullHandleStrategy])
///     - Initializer
/// -   [RegionData] ([bsx_rs::region::RegionData])
///     - Initializer
///     - drop_null
///     - [RegionCoordinates] getter
///     - discretize
use pyo3::prelude::*;
use pyo3_polars::*;

#[pyclass]
pub struct RegionCoordinates {
    inner: RegionCoordinatesRust,
}

#[pymethods]
impl RegionCoordinates {
    #[new]
    fn new(chr: String, start: u32, end: u32) -> Self {
        Self {
            inner: RegionCoordinatesRust::new(chr, start, end),
        }
    }

    #[getter]
    fn chr(&self) -> &str {
        self.inner.chr()
    }

    #[getter]
    fn start(&self) -> u32 {
        self.inner.start()
    }

    #[getter]
    fn end(&self) -> u32 {
        self.inner.end()
    }

    #[getter]
    fn length(&self) -> u32 {
        self.inner.length()
    }
}

// Python wrapper for RegionData
#[pyclass]
pub struct RegionData {
    inner: RegionDataRust,
}

#[pymethods]
impl RegionData {
    #[staticmethod]
    fn new(data: PyDataFrame, coordinates: &RegionCoordinates) -> PyResult<Self> {
        let data: DataFrame = data.into(); // convert from PyDataFrame to DataFrame
        let coordinates = coordinates.inner.clone();

        Ok(Self {
            inner: RegionDataRust::new(data, coordinates),
        })
    }

    fn drop_null(&mut self, null_strategy: NullHandleStrategy) -> Option<()> {
        self.inner.drop_null(null_strategy.inner)
    }

    fn get_coordinates(&self) -> RegionCoordinates {
        RegionCoordinates {
            inner: self.inner.get_coordinates().clone(),
        }
    }

    fn discretize(&self, resolution: usize) -> PyDataFrame {
        let df = self.inner.discretize(resolution);
        PyDataFrame(df)
    }

    fn __str__(&self) -> String {
        self.inner.to_string()
    }
}

#[pyclass]
#[derive(Debug, Clone)]
pub enum FillNullStrategy {
    /// mean value in array
    Mean,
    /// minimal value in array
    Min,
    /// maximum value in array
    Max,
    /// replace with the value zero
    Zero,
    /// replace with the value one
    One,
    /// replace with the maximum value of that data type
    MaxBound,
    /// replace with the minimal value of that data type
    MinBound,
}

impl FillNullStrategy {
    fn to_rust(&self) -> FillNullStrategyRust {
        match self {
            FillNullStrategy::Mean => FillNullStrategyRust::Mean,
            FillNullStrategy::Min => FillNullStrategyRust::Min,
            FillNullStrategy::Max => FillNullStrategyRust::Max,
            FillNullStrategy::Zero => FillNullStrategyRust::Zero,
            FillNullStrategy::One => FillNullStrategyRust::One,
            FillNullStrategy::MaxBound => FillNullStrategyRust::MaxBound,
            FillNullStrategy::MinBound => FillNullStrategyRust::MinBound,
        }
    }
}

#[pyclass]
#[derive(Clone)]
pub struct NullHandleStrategy {
    inner: NullHandleStrategyRust,
}

#[pymethods]
impl NullHandleStrategy {
    #[new]
    pub fn new(skip: bool, fill: FillNullStrategy) -> Self {
        Self {
            inner: NullHandleStrategyRust::new(skip, fill.to_rust()),
        }
    }
    #[staticmethod]
    pub fn default() -> Self {
        Self {
            inner: NullHandleStrategyRust::default(),
        }
    }

    pub fn with_fill(&self, fill: FillNullStrategy) -> Self {
        let new_inner = self.inner.clone().with_fill(fill.to_rust());
        Self { inner: new_inner }
    }

    pub fn with_skip(&self, skip: bool) -> Self {
        let new_inner = self.inner.clone().with_skip(skip);
        Self { inner: new_inner }
    }
}

#[cfg(test)]
mod tests {
    use super::{
        FillNullStrategy as FillNullStrategyPython, FillNullStrategyRust, NullHandleStrategy,
        RegionCoordinates, RegionData,
    };
    use polars::prelude::*;

    #[test]
    fn test_python_to_rust_enum_consistency() {
        // Verify each variant in the Python enum maps to the corresponding variant in the Rust enum.
        assert_eq!(
            FillNullStrategyPython::Mean.to_rust(),
            FillNullStrategyRust::Mean
        );
        assert_eq!(
            FillNullStrategyPython::Min.to_rust(),
            FillNullStrategyRust::Min
        );
        assert_eq!(
            FillNullStrategyPython::Max.to_rust(),
            FillNullStrategyRust::Max
        );
        assert_eq!(
            FillNullStrategyPython::Zero.to_rust(),
            FillNullStrategyRust::Zero
        );
        assert_eq!(
            FillNullStrategyPython::One.to_rust(),
            FillNullStrategyRust::One
        );
        assert_eq!(
            FillNullStrategyPython::MaxBound.to_rust(),
            FillNullStrategyRust::MaxBound
        );
        assert_eq!(
            FillNullStrategyPython::MinBound.to_rust(),
            FillNullStrategyRust::MinBound
        );
    }

    #[test]
    fn test_enum_variant_coverage() {
        // Ensure that every variant in the Rust enum has a corresponding Python mapping
        let rust_variants = vec![
            FillNullStrategyRust::Mean,
            FillNullStrategyRust::Min,
            FillNullStrategyRust::Max,
            FillNullStrategyRust::Zero,
            FillNullStrategyRust::One,
            FillNullStrategyRust::MaxBound,
            FillNullStrategyRust::MinBound,
        ];

        let python_variants = vec![
            FillNullStrategyPython::Mean.to_rust(),
            FillNullStrategyPython::Min.to_rust(),
            FillNullStrategyPython::Max.to_rust(),
            FillNullStrategyPython::Zero.to_rust(),
            FillNullStrategyPython::One.to_rust(),
            FillNullStrategyPython::MaxBound.to_rust(),
            FillNullStrategyPython::MinBound.to_rust(),
        ];

        // Check that all Rust variants are covered by Python enum mappings.
        assert_eq!(rust_variants, python_variants);
    }
    #[test]
    fn test_region_coordinates_properties() {
        let coordinates = RegionCoordinates::new("chr1".to_string(), 100, 200);

        // Test chr property
        assert_eq!(coordinates.chr(), "chr1");

        // Test start property
        assert_eq!(coordinates.start(), 100);

        // Test end property
        assert_eq!(coordinates.end(), 200);

        // Test length property
        assert_eq!(coordinates.length(), 100);
    }

    #[test]
    fn test_region_coordinates_length_calculation() {
        let coordinates = RegionCoordinates::new("chr2".to_string(), 300, 500);
        assert_eq!(coordinates.length(), 200); // 500 - 300
    }
}
