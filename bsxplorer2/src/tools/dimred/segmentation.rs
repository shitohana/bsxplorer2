use core::f64;

use crate::data_structs::typedef::CountType;

/// Trait defining the interface for data structures used in segmentation
/// algorithms.
pub trait SegmentationData {
    /// Calculates the cost of a segment from index r to t (inclusive).
    fn cost_function(
        &self,
        r: usize,
        t: usize,
    ) -> f64;
    /// Returns the total number of data points.
    fn len(&self) -> usize;

    fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

/// Data structure for binomial data segmentation, storing cumulative sums.
pub struct MethDataBinom {
    /// Cumulative sum of methylated counts.
    count_m_cumsum:     Vec<u32>,
    /// Cumulative sum of total counts.
    count_total_cumsum: Vec<u32>,
}

impl MethDataBinom {
    pub fn new(
        count_m: &[CountType],
        count_total: &[CountType],
    ) -> Self {
        assert_eq!(
            count_m.len(),
            count_total.len(),
            "count_m and count_total must have the same length"
        );
        // Calculate cumulative sum of methylated counts.
        let count_m_cumsum = count_m.iter().fold(vec![0], |mut acc, new| {
            acc.push(acc.last().unwrap() + *new as u32);
            acc
        });

        // Calculate cumulative sum of total counts.
        let count_total_cumsum = count_total.iter().fold(vec![0], |mut acc, new| {
            acc.push(acc.last().unwrap() + *new as u32);
            acc
        });

        Self {
            count_m_cumsum,
            count_total_cumsum,
        }
    }
}

impl SegmentationData for MethDataBinom {
    /// Calculates the binomial negative log-likelihood cost for segment [r, t].
    fn cost_function(
        &self,
        r: usize,
        t: usize,
    ) -> f64 {
        // Calculate sum of methylated counts in segment [r, t].
        let m = self.count_m_cumsum[t + 1] - self.count_m_cumsum[r];
        // Calculate sum of total counts in segment [r, t].
        let n = self.count_total_cumsum[t + 1] - self.count_total_cumsum[r];

        if n == 0 {
            return 0.0; // Avoid division by zero.
        }

        let p = m as f64 / n as f64;
        // Avoid log(0) or log(1).
        if p == 0.0 || p == 1.0 {
            return 0.0;
        }

        // Binomial negative log-likelihood cost (scaled).
        -2.0 * (m as f64 * p.ln() + (n - m) as f64 * (1.0 - p).ln())
    }

    /// Returns the number of data points.
    fn len(&self) -> usize {
        self.count_m_cumsum.len() - 1
    }
}

/// PELT algorithm for segmentation based on doi:10.1080/01621459.2012.737745
///
/// # Arguments
///
/// * `data` - A [SegmentationData] object
/// * `beta` - The penalty parameter
/// * `min_size` - The minimum size of a segment
///
/// # Returns
///
/// A vector of change points
pub fn pelt<T: SegmentationData>(
    data: &T,
    beta: f64,
    min_size: usize,
) -> (Vec<usize>, f64) {
    let n = data.len();

    #[allow(non_snake_case)]
    // F[i] stores the minimum cost to segment data up to index i-1.
    let mut F = vec![f64::INFINITY; n + 1];
    // prev[i] stores the index of the last changepoint before i-1 in the optimal
    // segmentation.
    let mut prev = vec![-1; n + 1];
    // candidate_set stores potential previous changepoints due to pruning.
    let mut candidate_set = vec![0]; // Start with a virtual changepoint at index 0.

    // Base case: cost to segment an empty set before index 0.
    F[0] = -beta; // Subtract beta for consistency with the recursive definition.

    // Iterate through all possible end points 't' (0 to n-1). F[t+1] corresponds to
    // data up to index t.
    for t in 0..n {
        let mut best_cost = f64::INFINITY;
        let mut best_r = -1;

        // Iterate through candidate previous changepoints 'r'.
        for r in candidate_set.iter() {
            // Ensure minimum segment size requirement is met.
            if (t + 1 - *r) >= min_size {
                // Calculate cost of segmenting up to t, with a changepoint at r.
                // Cost = Min cost up to r + Cost of segment [r, t] + Penalty.
                let c = F[*r] + data.cost_function(*r, t) + beta;
                if c < best_cost {
                    best_cost = c;
                    best_r = *r as i32;
                }
            }
        }

        // Store the minimum cost and the corresponding previous changepoint for index
        // t+1.
        F[t + 1] = best_cost;
        prev[t + 1] = best_r;

        // Pruning step: Remove candidates 'r' that cannot be the optimal previous
        // changepoint for any future end point t' >= t, based on the pruning
        // condition.
        candidate_set.retain(|r| {
            // F[r] + cost(r, t) <= F[t + 1] - This is a simplified version of the
            // pruning inequality as implemented here. The full PELT
            // condition involves costs for future points.
            // This implementation uses the efficient form F[r] + cost(r,t) <= F[t+1].
            // More robust pruning requires comparing F[r] + cost(r, t_prime) vs F[s] +
            // cost(s, t_prime). This specific condition prunes candidates
            // that are worse than the current optimum ending at t+1.
            F[*r] + data.cost_function(*r, t) <= F[t + 1]
        });
        // Add the current index 't' as a potential candidate for future iterations.
        candidate_set.push(t + 1); // The candidate represents a possible split
                                   // point *after* index t.
    }

    // Traceback to find the changepoints.
    let mut cps = Vec::new();
    let mut t = n as i32;
    while t > 0 {
        let r = prev[t as usize]; // Get the previous changepoint index.
        if r < 0 {
            break; // Stop when we reach the virtual changepoint at index 0.
        }
        // The stored index 'r' is the start of the last segment ending at t-1.
        // This 'r' is the changepoint *before* the segment starting at 'r'.
        cps.push(r as usize); // Add the changepoint index.
        t = r; // Move to the previous changepoint.
    }

    // Sort changepoints in ascending order. Remove the initial 0 if it's added by
    // the loop logic. The initial candidate_set = vec![0] can lead to 0 being
    // added to cps if prev[n] traces back to 0. Filter out 0 if it is added as
    // a changepoint when it's just the start boundary.
    cps.retain(|&x| x > 0); // Retain only indices > 0 as actual changepoints.

    cps.sort();
    // Return the sorted changepoints and the minimum cost for the full
    // segmentation.
    (cps, F[n])
}

pub enum SegmentAlgorithm {
    /// Pelt algorithm with params beta and minimum segment size.
    /// If beta is None it is set to BIC: ln(len(data))
    Pelt(Option<f64>, usize),
}
