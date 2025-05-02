use core::f64;
use std::cmp::Ordering;


trait SegmentationData {
    fn cost_function(
        &self,
        r: usize,
        t: usize,
    ) -> f64;
    fn len(&self) -> usize;
}

struct MethDataBinom {
    count_m_cumsum: Vec<u32>,
    count_total_cumsum: Vec<u32>,
}

impl MethDataBinom {
    pub fn new(
        count_m: &[u32],
        count_total: &[u32],
    ) -> Self {
        assert_eq!(
            count_m.len(),
            count_total.len(),
            "count_m and count_total must have the same length"
        );
        let count_m_cumsum = count_m
            .iter()
            .fold(vec![0], |mut acc, new| {
                acc.push(acc.last().unwrap() + *new);
                acc
            });

        let count_total_cumsum =
            count_total
                .iter()
                .fold(vec![0], |mut acc, new| {
                    acc.push(acc.last().unwrap() + *new);
                    acc
                });

        Self {
            count_m_cumsum,
            count_total_cumsum,
        }
    }
}

impl SegmentationData for MethDataBinom {
    fn cost_function(
        &self,
        r: usize,
        t: usize,
    ) -> f64 {
        let m = self.count_m_cumsum[t + 1] - self.count_m_cumsum[r];
        let n = self.count_total_cumsum[t + 1] - self.count_total_cumsum[r];

        if n == 0 {
            return 0.0;
        }

        let p = m as f64 / n as f64;
        if p == 0.0 || p == 1.0 {
            return 0.0;
        }

        -2.0 * (m as f64 * p.ln() + (n - m) as f64 * (1.0 - p).ln())
    }

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
fn pelt<T: SegmentationData>(
    data: &T,
    beta: f64,
    min_size: usize,
) -> (Vec<usize>, f64) {
    let n = data.len();

    #[allow(non_snake_case)]
    let mut F = vec![f64::INFINITY; n + 1];
    let mut prev = vec![-1; n + 1];
    let mut candidate_set = vec![0];

    F[0] = -beta;

    for t in 0..n {
        let mut best_cost = f64::INFINITY;
        let mut best_r = -1;

        for r in candidate_set.iter() {
            if (t + 1 - *r) >= min_size {
                let c = F[*r] + data.cost_function(*r, t) + beta;
                if c < best_cost {
                    best_cost = c;
                    best_r = *r as i32;
                }
            }
        }

        F[t + 1] = best_cost;
        prev[t + 1] = best_r;

        candidate_set.retain(|r| F[*r] + data.cost_function(*r, t) <= F[t + 1]);
        candidate_set.push(t);
    }

    let mut cps = Vec::new();
    let mut t = n as i32;
    while t > 0 {
        let r = prev[t as usize];
        if r < 0 {
            break;
        }
        cps.push(r as usize);
        t = r;
    }

    cps.sort();
    (cps, F[n])
}

enum MethSegmentor {
    Pelt(Option<f64>, usize),
}

struct SegmentBoundary {
    chr_num: usize,
    pos: u32,
}

impl PartialOrd for SegmentBoundary {
    fn partial_cmp(
        &self,
        other: &Self,
    ) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for SegmentBoundary {
    fn cmp(
        &self,
        other: &Self,
    ) -> Ordering {
        self.chr_num
            .cmp(&other.chr_num)
            .then(self.pos.cmp(&other.pos))
    }
}

impl PartialEq for SegmentBoundary {
    fn eq(
        &self,
        other: &Self,
    ) -> bool {
        self.chr_num == other.chr_num && self.pos == other.pos
    }
}

impl Eq for SegmentBoundary {}
