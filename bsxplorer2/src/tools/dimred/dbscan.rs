//! This is a slightly modified copy of dbscan crate to avoid allocating
//! vector to just a single point. So this is 1D DBSCAN implementation

use rayon::slice::ParallelSliceMut;
// MIT License
//
// Copyright (c) 2018 Michael Lazear
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
use Classification::{
    Core,
    Edge,
    Noise,
};

use crate::data_structs::typedef::PosType;

#[derive(Debug, Copy, Clone, PartialEq, PartialOrd)]
pub enum Classification {
    /// A point with at least `min_points` neighbors within `eps` diameter
    Core(usize),
    /// A point within `eps` of a core point, but has less than `min_points`
    /// neighbors
    Edge(usize),
    /// A point with no connections
    Noise,
}

/// DBSCAN parameters
pub struct Model {
    /// Epsilon value - maximum distance between points in a cluster
    pub eps: PosType,
    /// Minimum number of points in a cluster
    pub mpt: usize,

    c: Vec<Classification>,
    v: Vec<bool>,
}

impl Clone for Model {
    fn clone(&self) -> Self {
        Self::new(self.eps, self.mpt)
    }
}

impl Model {
    /// Create a new `Model` with a set of parameters
    ///
    /// # Arguments
    /// * `eps` - maximum distance between datapoints within a cluster
    /// * `min_points` - minimum number of datapoints to make a cluster
    pub fn new(
        eps: PosType,
        min_points: usize,
    ) -> Model {
        Model {
            eps,
            mpt: min_points,
            c: Vec::new(),
            v: Vec::new(),
        }
    }

    fn expand_cluster(
        &mut self,
        population: &[PosType],
        queue: &mut Vec<usize>,
        cluster: usize,
    ) -> bool {
        let mut new_cluster = false;
        while let Some(ind) = queue.pop() {
            let neighbors = self.range_query(population[ind], population);
            if neighbors.len() < self.mpt {
                continue;
            }
            new_cluster = true;
            self.c[ind] = Core(cluster);
            for n_idx in neighbors {
                // n_idx is at least an edge point
                if self.c[n_idx] == Noise {
                    self.c[n_idx] = Edge(cluster);
                }

                if self.v[n_idx] {
                    continue;
                }

                self.v[n_idx] = true;
                queue.push(n_idx);
            }
        }
        new_cluster
    }

    #[inline]
    fn range_query(
        &self,
        sample: PosType,
        population: &[PosType],
    ) -> Vec<usize> {
        let value_idx = population.binary_search(&sample).unwrap_or_else(|e| e);

        let len = population.len();
        let (min, max) = unsafe {
            let mut cur_idx = value_idx;
            let min_idx = loop {
                if cur_idx == 0 {
                    break 0;
                }
                let val = population.get_unchecked(cur_idx);
                if val.abs_diff(sample) < self.eps {
                    cur_idx -= 1;
                }
                else {
                    break cur_idx + 1;
                }
            };
            let mut cur_idx = value_idx;
            let max_idx = loop {
                if cur_idx >= len {
                    break len - 1;
                }
                let val = population.get_unchecked(cur_idx);
                if val.abs_diff(sample) < self.eps {
                    cur_idx += 1;
                }
                else {
                    break cur_idx - 1;
                }
            };

            (min_idx, max_idx)
        };

        (min..=max).collect()
    }

    pub fn run(
        mut self,
        population: &[PosType],
    ) -> Vec<Classification> {
        self.c = vec![Noise; population.len()];
        self.v = vec![false; population.len()];

        let population = if !population.is_sorted() {
            let mut pop_clone = population.to_vec();
            pop_clone.par_sort_unstable();
            pop_clone
        }
        else {
            population.to_vec()
        };

        let mut cluster = 0;
        let mut queue = Vec::new();

        for idx in 0..population.len() {
            if self.v[idx] {
                continue;
            }

            self.v[idx] = true;

            queue.push(idx);

            if self.expand_cluster(&population, &mut queue, cluster) {
                cluster += 1;
            }
        }
        self.c
    }
}
