/// ***********************************************************************
/// *****
/// * Copyright (c) 2025
/// The Prosperity Public License 3.0.0
///
/// Contributor: [shitohana](https://github.com/shitohana)
///
/// Source Code: https://github.com/shitohana/BSXplorer
/// ***********************************************************************
/// ****
use num::Float;

enum DistanceMetric {
    Eucledian,
    Cosine,
    Canberra,
    Minkowsky(usize),
}

trait MetricAccumulator<N: Float>: Clone {
    fn finalize(&self) -> N;
    fn add_ppoints(
        &mut self,
        xi: N,
        yi: N,
    );
}

#[derive(Clone)]
pub struct MinkowskiAccum<N: Float> {
    power: usize,
    powsum: N,
}

impl<N: Float> MinkowskiAccum<N> {
    fn new(power: usize) -> Self {
        Self {
            power,
            powsum: N::from(0.0).unwrap(),
        }
    }
}

impl<N: Float> MetricAccumulator<N> for MinkowskiAccum<N> {
    fn finalize(&self) -> N {
        self.powsum
            .powf(N::from(1f64 / self.power as f64).unwrap())
    }

    fn add_ppoints(
        &mut self,
        xi: N,
        yi: N,
    ) {
        self.powsum = self.powsum + (xi - yi).abs().powi(self.power as i32)
    }
}

#[derive(Clone)]
pub struct CanberraAccum<N: Float> {
    sum: N,
}

impl<N: Float> CanberraAccum<N> {
    fn new() -> Self {
        Self {
            sum: N::from(0.0).unwrap(),
        }
    }
}

impl<N: Float> MetricAccumulator<N> for CanberraAccum<N> {
    fn finalize(&self) -> N {
        self.sum
    }

    fn add_ppoints(
        &mut self,
        xi: N,
        yi: N,
    ) {
        self.sum = self.sum + (xi - yi).abs() / (xi.abs() + yi.abs());
    }
}

#[derive(Clone)]
pub struct EucledianAccum<N: Float> {
    sqsum: N,
}

impl<N: Float> EucledianAccum<N> {
    fn new() -> Self {
        Self {
            sqsum: N::from(0.0).unwrap(),
        }
    }
}

impl<N: Float> MetricAccumulator<N> for EucledianAccum<N> {
    fn finalize(&self) -> N {
        self.sqsum.sqrt()
    }

    fn add_ppoints(
        &mut self,
        xi: N,
        yi: N,
    ) {
        self.sqsum = self.sqsum + (xi - yi).powi(2)
    }
}

#[derive(Clone)]
pub struct CosineAccum<N: Float> {
    dot_prod: N,
    right_sqsum: N,
    left_sqsum: N,
}

impl<N: Float> CosineAccum<N> {
    fn new() -> Self {
        Self {
            dot_prod: N::from(0.0).unwrap(),
            right_sqsum: N::from(0.0).unwrap(),
            left_sqsum: N::from(0.0).unwrap(),
        }
    }
}

impl<N: Float> MetricAccumulator<N> for CosineAccum<N> {
    fn finalize(&self) -> N {
        N::from(1.0).unwrap()
            - self.dot_prod / (self.right_sqsum.sqrt() * self.left_sqsum.sqrt())
    }

    fn add_ppoints(
        &mut self,
        xi: N,
        yi: N,
    ) {
        self.dot_prod = self.dot_prod + xi * yi;
        self.right_sqsum = self.right_sqsum + xi.powi(2);
        self.left_sqsum = self.left_sqsum + yi.powi(2);
    }
}

#[cfg(test)]
mod tests {
    use assert_approx_eq::assert_approx_eq;

    use super::*;

    #[test]
    fn test_minkowski_accum() {
        // Test with p=1 (Manhattan distance)
        let mut accum = MinkowskiAccum::<f64>::new(1);
        accum.add_ppoints(1.0, 4.0);
        accum.add_ppoints(2.0, 6.0);
        assert_approx_eq!(accum.finalize(), 7.0); // |1-4| + |2-6| = 3 + 4 = 7

        // Test with p=2 (Euclidean distance)
        let mut accum = MinkowskiAccum::<f64>::new(2);
        accum.add_ppoints(0.0, 3.0);
        accum.add_ppoints(0.0, 4.0);
        assert_approx_eq!(accum.finalize(), 5.0); // sqrt((0-3)^2 + (0-4)^2) =
                                                  // sqrt(9 + 16) = 5
    }

    #[test]
    fn test_canberra_accum() {
        let mut accum = CanberraAccum::<f64>::new();
        accum.add_ppoints(1.0, 3.0);
        accum.add_ppoints(2.0, 5.0);
        // |1-3|/(|1|+|3|) + |2-5|/(|2|+|5|) = 2/4 + 3/7 = 0.5 + 0.428... = 0.928...
        assert_approx_eq!(accum.finalize(), 0.5 + 3.0 / 7.0, 1e-10);
    }

    #[test]
    fn test_eucledian_accum() {
        let mut accum = EucledianAccum::<f64>::new();
        accum.add_ppoints(1.0, 4.0);
        accum.add_ppoints(2.0, 6.0);
        // sqrt((1-4)^2 + (2-6)^2) = sqrt(9 + 16) = sqrt(25) = 5
        assert_approx_eq!(accum.finalize(), 5.0);
    }

    #[test]
    fn test_cosine_accum() {
        let mut accum = CosineAccum::<f64>::new();
        accum.add_ppoints(1.0, 0.0);
        accum.add_ppoints(0.0, 1.0);
        // dot product = 0, |a| = 1, |b| = 1, cosine similarity = 0, distance =
        // 1-0 = 1
        assert_approx_eq!(accum.finalize(), 1.0);

        // Test case for identical vectors (should give distance 0)
        let mut accum = CosineAccum::<f64>::new();
        accum.add_ppoints(2.0, 2.0);
        accum.add_ppoints(3.0, 3.0);
        assert_approx_eq!(accum.finalize(), 0.0, 1e-10);
    }

    #[test]
    fn test_with_different_numeric_types() {
        // Test with f32
        let mut accum = EucledianAccum::<f32>::new();
        accum.add_ppoints(1.0, 4.0);
        assert_approx_eq!(accum.finalize(), 3.0);

        // Test with f64
        let mut accum = EucledianAccum::<f64>::new();
        accum.add_ppoints(1.0, 4.0);
        assert_approx_eq!(accum.finalize(), 3.0);
    }

    #[test]
    fn test_clone_functionality() {
        let mut accum = EucledianAccum::<f64>::new();
        accum.add_ppoints(1.0, 4.0);

        let cloned = accum.clone();
        assert_approx_eq!(cloned.finalize(), 3.0);

        // Ensure original is unchanged
        assert_approx_eq!(accum.finalize(), 3.0);

        // Modify original after cloning
        accum.add_ppoints(2.0, 6.0);
        assert_approx_eq!(accum.finalize(), 5.0);

        // Ensure clone is still the same
        assert_approx_eq!(cloned.finalize(), 3.0);
    }
}
