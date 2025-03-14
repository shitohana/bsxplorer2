/// ****************************************************************************
/// * Copyright (c) 2025
/// The Prosperity Public License 3.0.0
///
/// Contributor: [shitohana](https://github.com/shitohana)
///
/// Source Code: https://github.com/shitohana/BSXplorer
/// ***************************************************************************

/// ****************************************************************************
/// * Copyright (c) 2025
/// ***************************************************************************

#[cfg(feature = "plots")]
mod inner {
    use std::collections::BTreeMap;

    use anyhow::anyhow;
    use histogram::{Bucket, Histogram};
    use itertools::{izip, Itertools};
    use num::ToPrimitive;
    use plotly::common::{ErrorData, Mode};
    use plotly::layout::Axis;
    use plotly::{Plot, Scatter};
    use statrs::distribution::{ContinuousCDF, Normal};

    use crate::data_structs::bsx_batch::EncodedBsxBatch;
    use crate::data_structs::region_data::RegionData;
    use crate::plots::BsxPlot;
    use crate::utils::types::{PosNum, RefId, Strand};
    use crate::utils::{f64_to_u64_scaled, u64_to_f64_scaled, GROUPING_POWER};

    fn bucket_middle(bucket: &Bucket) -> u64 {
        bucket.start() / 2 + bucket.end() / 2
    }

    type BoxPlotPoints = (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>);

    pub trait LinePlot: BsxPlot {
        fn hist_iter(&self) -> impl Iterator<Item = &Histogram>;
        fn size(&self) -> usize;

        fn mean(&self) -> Vec<f64> {
            self.hist_iter()
                .map(|hist| {
                    let sum_density: f64 = hist
                        .iter()
                        .map(|bucket| {
                            let u64_val =
                                (bucket.start() + 1) / 2 + bucket.end() / 2;
                            u64_to_f64_scaled(u64_val) * bucket.count() as f64
                        })
                        .sum();
                    let counts: u64 = hist
                        .iter()
                        .map(|bucket| bucket.count())
                        .sum();
                    sum_density / counts as f64
                })
                .collect_vec()
        }
        fn std(
            &self,
            mean: &[f64],
        ) -> Vec<f64> {
            self.hist_iter()
                .zip(mean.iter())
                .map(|(hist, mean)| {
                    let mut count_sum = 0f64;
                    let mut dev_sum = 0f64;
                    hist.iter().for_each(|bin| {
                        if bin.count() > 0 {
                            let count = bin.count() as f64;
                            dev_sum += count
                                * (u64_to_f64_scaled(bucket_middle(&bin))
                                    - mean)
                                    .powi(2);
                            count_sum += count
                        }
                    });
                    dev_sum / count_sum
                })
                .collect()
        }
        fn confidence_interval(
            &self,
            confidence: f64,
        ) -> (Vec<f64>, Vec<f64>) {
            let mean = self.mean();
            let std = self.std(&mean);
            let sample_size: Vec<u64> = self
                .hist_iter()
                .map(|hist| hist.iter().map(|b| b.count()).sum())
                .collect_vec();
            let ci = izip!(mean.iter(), std.iter(), sample_size.iter())
                .map(|(m, s, size)| {
                    let z_score = Normal::new(0.0, 1.0)
                        .unwrap()
                        .inverse_cdf(1.0 - (1.0 - confidence) / 2.0);
                    let margin_of_error = z_score * (s / (*size as f64).sqrt());
                    margin_of_error
                })
                .collect_vec();
            (mean, ci)
        }

        fn percentile(
            &self,
            percentile: f64,
        ) -> Vec<f64> {
            self.hist_iter()
                .map(|h| {
                    h.percentile(percentile)
                        .expect("Failed to get percentile")
                        .map(|val| u64_to_f64_scaled(bucket_middle(&val)))
                        .unwrap_or(0.0)
                })
                .collect_vec()
        }
        fn boxplot_points(&self) -> BoxPlotPoints {
            self.hist_iter()
                .map(|hist| {
                    let median = u64_to_f64_scaled(
                        hist.percentile(50.0)
                            .unwrap()
                            .map(|b| bucket_middle(&b))
                            .unwrap_or(0),
                    );
                    let q1 = u64_to_f64_scaled(
                        hist.percentile(25.0)
                            .unwrap()
                            .map(|b| bucket_middle(&b))
                            .unwrap_or(0),
                    );
                    let q3 = u64_to_f64_scaled(
                        hist.percentile(75.0)
                            .unwrap()
                            .map(|b| bucket_middle(&b))
                            .unwrap_or(0),
                    );
                    let iqr = q3 - q1;
                    let min = q1 - iqr;
                    let max = q3 + iqr;
                    (min, q1, median, q3, max)
                })
                .multiunzip()
        }

        fn x_ticks(&self) -> Vec<f64>;

        fn plot(
            &self,
            plot: &mut Plot,
            label: impl AsRef<str>,
            ci_prob: Option<f64>,
        ) {
            let (mean, ci) = if let Some(confidence) = ci_prob {
                let res = self.confidence_interval(confidence);
                (res.0, Some(res.1))
            }
            else {
                (self.mean(), None)
            };
            let ticks = self.x_ticks();
            let mut main_trace = Scatter::new(ticks.clone(), mean.clone())
                .name(&label)
                .mode(Mode::Lines);

            if let Some(ci) = ci {
                let error_data = ErrorData::default()
                    .array(ci)
                    .symmetric(true)
                    .visible(true)
                    .copy_ystyle(true);

                main_trace = main_trace.error_y(error_data)
            }

            plot.add_trace(main_trace)
        }
    }

    #[derive(Debug)]
    pub struct LinePlotData {
        discrete_hist: Vec<Histogram>,
    }

    impl BsxPlot for LinePlotData {
        fn add_region_data<R: RefId, N: PosNum>(
            &mut self,
            region_data: &RegionData<R, N, EncodedBsxBatch>,
        ) -> anyhow::Result<()> {
            let mut discrete_density =
                region_data.discrete_density(self.discrete_hist.len())?;
            if matches!(region_data.strand(), Strand::Reverse) {
                discrete_density.reverse()
            };
            self.add_density(&discrete_density)
        }
    }

    impl LinePlot for LinePlotData {
        fn hist_iter(&self) -> impl Iterator<Item = &Histogram> {
            self.discrete_hist.iter()
        }

        fn size(&self) -> usize { self.discrete_hist.len() }

        fn x_ticks(&self) -> Vec<f64> {
            (0..self.discrete_hist.len())
                .map(|x| x as f64 + 0.5)
                .collect_vec()
        }
    }

    impl LinePlotData {
        pub fn new(resolution: usize) -> LinePlotData {
            LinePlotData {
                discrete_hist: (0..resolution)
                    .map(|_| Histogram::new(GROUPING_POWER, 64).unwrap())
                    .collect(),
            }
        }

        pub fn add_density<V: Into<f64> + Copy>(
            &mut self,
            density_discrete: &[V],
        ) -> anyhow::Result<()> {
            if density_discrete.len() > self.discrete_hist.len() {
                return Err(anyhow!(histogram::Error::IncompatibleParameters));
            }
            density_discrete
                .iter()
                .enumerate()
                .map(|(idx, val)| {
                    let conv: f64 = <V as Into<f64>>::into(*val);
                    self.discrete_hist[idx].increment(f64_to_u64_scaled(conv))
                })
                .collect::<Result<Vec<_>, histogram::Error>>()?;
            Ok(())
        }
    }

    pub fn add_line(
        plot: &mut Plot,
        x: &[f64],
        y: &[f64],
        label: impl AsRef<str>,
    ) {
        let trace = Scatter::new(x.to_vec(), y.to_vec())
            .name(label)
            .mode(Mode::Lines);
        plot.add_trace(trace);
    }

    pub struct GenomeWideLinePlotData<R: RefId, N: PosNum> {
        bins_hist: BTreeMap<(R, N), Histogram>,
        bin_width: u32,
    }

    impl<R: RefId, N: PosNum> BsxPlot for GenomeWideLinePlotData<R, N> {
        fn add_region_data<A: RefId, B: PosNum>(
            &mut self,
            data: &RegionData<A, B, EncodedBsxBatch>,
        ) -> anyhow::Result<()> {
            let ref_id = data.chr();
            for (pos, density) in data
                .get_positions()?
                .into_iter()
                .zip(data.meth_density()?.into_iter())
            {
                let entry = self
                    .bins_hist
                    .entry((
                        format!("{ref_id}").into(),
                        N::from(
                            (pos.to_f64().unwrap() / self.bin_width as f64)
                                .floor(),
                        )
                        .unwrap(),
                    ))
                    .or_insert(Histogram::new(GROUPING_POWER, 64)?);
                entry.increment(f64_to_u64_scaled(density as f64))?
            }
            Ok(())
        }
    }

    impl<R: RefId, N: PosNum> LinePlot for GenomeWideLinePlotData<R, N> {
        fn hist_iter(&self) -> impl Iterator<Item = &Histogram> {
            self.bins_hist.values()
        }

        fn size(&self) -> usize { self.bins_hist.len() }

        fn x_ticks(&self) -> Vec<f64> {
            (0..self.bins_hist.len())
                .into_iter()
                .map(|x| x as f64)
                .collect_vec()
        }
    }

    impl<R: RefId, N: PosNum> GenomeWideLinePlotData<R, N> {
        pub fn new(bin_width: u32) -> Self {
            Self {
                bin_width,
                bins_hist: BTreeMap::new(),
            }
        }

        pub fn x_labels(&self) -> Vec<(&R, f64)> {
            self.bins_hist
                .keys()
                .map(|(ref_id, v)| (ref_id, v.to_f64().unwrap()))
                .collect_vec()
        }

        pub fn draw_x_labels(
            &self,
            plot: &mut Plot,
        ) {
            let mut x_ticks = Vec::new();
            let mut seen_chr = Vec::new();

            for (chr, tick) in self.bins_hist.keys() {
                if !seen_chr.contains(chr) {
                    x_ticks.push(tick.clone());
                    seen_chr.push(chr.clone());
                }
            }

            x_ticks.push(
                self.bins_hist
                    .last_key_value()
                    .unwrap()
                    .0
                    .clone()
                    .1,
            );
            let x_ticks = x_ticks
                .windows(2)
                .map(|arr| (arr[0] + arr[1]).to_f64().unwrap() / 2.0)
                .collect_vec();
            let mut layout = plot.layout().clone();
            layout = layout.x_axis(
                Axis::default()
                    .tick_values(x_ticks)
                    .tick_text(
                        seen_chr
                            .into_iter()
                            .map(|chr| chr.to_string())
                            .collect(),
                    ),
            );

            plot.set_layout(layout);
        }
    }
}

pub use inner::{add_line, LinePlot, LinePlotData};
