use crate::data_structs::typedef::*;
use crate::{
    with_field_fn,
    Context,
};

#[derive(Debug, Clone)]
pub struct DmrConfig {
    pub context:        Context,
    pub n_missing:      usize,
    pub min_coverage:   CountType,
    pub diff_threshold: DensityType,
    pub min_cpgs:       usize,
    pub max_dist:       PosType,
    pub initial_l:      f64,
    pub l_min:          f64,
    pub l_coef:         f64,
    pub seg_tolerance:  DensityType,
    pub merge_pvalue:   f64,
    pub seg_pvalue:     f64,
}

impl DmrConfig {
    with_field_fn!(context, Context);

    with_field_fn!(n_missing, usize);

    with_field_fn!(min_coverage, CountType);

    with_field_fn!(diff_threshold, DensityType);

    with_field_fn!(min_cpgs, usize);

    with_field_fn!(max_dist, PosType);

    with_field_fn!(initial_l, f64);

    with_field_fn!(l_min, f64);

    with_field_fn!(l_coef, f64);

    with_field_fn!(seg_tolerance, DensityType);

    with_field_fn!(merge_pvalue, f64);

    with_field_fn!(seg_pvalue, f64);
}

impl Default for DmrConfig {
    fn default() -> Self {
        Self {
            context:        Context::CG,
            n_missing:      0,
            min_coverage:   5,
            diff_threshold: 0.1,
            min_cpgs:       10,
            max_dist:       100,
            initial_l:      2.0,
            l_min:          1e-3,
            l_coef:         1.5,
            seg_tolerance:  1e-6,
            merge_pvalue:   1e-3,
            seg_pvalue:     1e-2,
        }
    }
}
