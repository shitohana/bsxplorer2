use std::cmp::Ordering;

use hashbrown::HashMap;
use polars::prelude::*;
use rstest::{
    fixture,
    rstest,
};
use statrs::statistics::{
    OrderStatistics,
    Statistics,
};

use super::BsxColumns as BsxCol;
use crate::data_structs::batch::{
    create_caregorical_dtype,
    *,
};
use crate::prelude::*;
#[cfg(feature = "tools")]
use crate::tools::dimred::*;
use crate::utils::get_categorical_dtype;

mod batch_tests {
    use super::*;
    use crate::data_structs::typedef::{
        DensityType,
        PosType,
    };

    #[fixture]
    fn real_batch() -> BsxBatch {
        use std::fs::File;
        use std::path::PathBuf;

        use crate::io::bsx::BsxFileReader;

        let reader = BsxFileReader::try_new(
            File::open(
                PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data/report.bsx"),
            )
            .unwrap(),
        )
        .unwrap();
        reader.into_iter().next().unwrap().unwrap()
    }

    #[fixture]
    fn real_batch_slice_1k(real_batch: BsxBatch) -> BsxBatch {
        let batch = real_batch.lazy()
            .filter_coverage_gt(0) // Ensure no nulls in density for partition check
            .collect()
            .unwrap();
        let sliced_batch = batch.slice(0, 1000);
        sliced_batch
    }

    #[fixture]
    pub fn test_batch() -> BsxBatch {
        BsxBatch::try_from_columns(
            "chr1",
            Some(create_caregorical_dtype(vec!["chr1".into()])),
            vec![3, 5, 9, 12, 15],
            vec![true, false, true, true, false],
            vec![Some(true), Some(false), Some(true), None, None],
            vec![5, 10, 15, 10, 5],
            vec![10, 30, 20, 10, 10],
        )
        .unwrap()
    }

    #[test]
    fn test_empty_batch() {
        let batch = BsxBatch::empty(None);
        assert!(batch.is_empty());
    }

    #[rstest]
    fn test_partition(real_batch: BsxBatch) -> anyhow::Result<()> {
        let real_batch = real_batch
            .lazy()
            .filter_coverage_gt(0)
            .collect()
            .unwrap()
            .slice(0, 1000);
        assert_eq!(real_batch.len(), 1000);

        let batch_size = real_batch.len();
        let agg_fn = Box::new(|slice: &[f32]| {
            slice.iter().sum::<f32>() as f64 / slice.len() as f64
        });

        // Case 1: Empty breakpoints
        let (rel_pos, densities) = real_batch.partition(vec![], agg_fn.clone())?;
        assert!(rel_pos.is_empty());
        assert!(densities.is_empty());

        // Case 2: Breakpoints covering the whole batch (excluding 0)
        // Expect one segment, relative position 1.0, density is mean of whole batch
        let breakpoints_full = vec![batch_size];
        let (rel_pos_full, densities_full) =
            real_batch.partition(breakpoints_full, agg_fn.clone())?;
        assert_eq!(rel_pos_full, vec![1.0]);
        assert_eq!(densities_full.len(), 1);
        let expected_density_full =
            real_batch.density().drop_nulls().mean().unwrap_or(f64::NAN) as f64;
        // Using assert_approx_eq due to float calculations
        use assert_approx_eq::assert_approx_eq;
        assert_approx_eq!(densities_full[0], expected_density_full, 1e-6);

        // Case 3: Valid breakpoints (e.g., splitting into two halves)
        let mid_pos_index = batch_size / 2;
        let breakpoints_half = vec![mid_pos_index];
        let (rel_pos_half, densities_half) =
            real_batch.partition(breakpoints_half, agg_fn.clone())?;
        assert_eq!(rel_pos_half.len(), 2);
        assert_eq!(densities_half.len(), 2);
        assert_approx_eq!(rel_pos_half.last().copied().unwrap(), 1.0, 1e-6);

        // Check relative position calculation
        let first_pos = real_batch.first_pos().unwrap() as f64;
        let last_pos = real_batch.last_pos().unwrap() as f64;
        let total_length = last_pos - first_pos + 1.0;
        let mid_pos_value = real_batch.position().get(mid_pos_index).unwrap() as f64;
        let expected_rel_mid = (mid_pos_value - first_pos) / total_length;
        assert_approx_eq!(rel_pos_half[0], expected_rel_mid, 1e-6);

        // Check density calculations (mean of each half)
        let first_half_density_mean = real_batch
            .slice(0, mid_pos_index)
            .density()
            .drop_nulls()
            .mean()
            .unwrap_or(f64::NAN) as f64;
        let second_half_density_mean = real_batch
            .slice(mid_pos_index as i64, batch_size - mid_pos_index)
            .density()
            .drop_nulls()
            .mean()
            .unwrap_or(f64::NAN) as f64;
        assert_approx_eq!(densities_half[0], first_half_density_mean, 1e-6);
        assert_approx_eq!(densities_half[1], second_half_density_mean, 1e-6);


        // Case 4: Multiple breakpoints
        let breakpoints_multi =
            vec![batch_size / 4, batch_size / 2, 3 * batch_size / 4];
        let (rel_pos_multi, densities_multi) =
            real_batch.partition(breakpoints_multi, agg_fn.clone())?;
        assert_eq!(rel_pos_multi.len(), 4);
        assert_eq!(densities_multi.len(), 4);
        assert_approx_eq!(rel_pos_multi.last().copied().unwrap(), 1.0, 1e-6);

        // Case 5: Breakpoints that already include batch_size
        let breakpoints_with_size = vec![batch_size / 3, batch_size];
        let (rel_pos_with_size, densities_with_size) =
            real_batch.partition(breakpoints_with_size, agg_fn.clone())?;
        assert_eq!(rel_pos_with_size.len(), 2);
        assert_eq!(densities_with_size.len(), 2);
        assert_approx_eq!(rel_pos_with_size.last().copied().unwrap(), 1.0, 1e-6);


        // Case 6: Invalid breakpoint (index 0) - should bail
        let breakpoints_zero = vec![0, 100];
        let result_zero = real_batch.partition(breakpoints_zero, agg_fn.clone());
        assert!(result_zero.is_err());
        assert_eq!(
            result_zero.unwrap_err().to_string(),
            "Partition index can not be 0"
        );

        // Case 7: Invalid breakpoint (out of bounds) - should bail
        let breakpoints_oob = vec![batch_size / 2, batch_size + 1];
        let result_oob = real_batch.partition(breakpoints_oob, agg_fn.clone());
        assert!(result_oob.is_err());
        assert_eq!(
            result_oob.unwrap_err().to_string(),
            "Partition index out of bounds"
        );
        Ok(())
    }

    #[rstest]
    fn test_partition_otherfn(test_batch: BsxBatch) -> anyhow::Result<()> {
        use assert_approx_eq::assert_approx_eq;

        let agg_sum_fn = Box::new(|slice: &[f32]| slice.iter().sum::<f32>() as f64);
        let batch_size = test_batch.len(); // Use test_batch size
        let mid_pos_index = batch_size / 2; // 5 / 2 = 2
        let breakpoints_sum = vec![mid_pos_index]; // breakpoints = [2] -> segments [0..2), [2..5)
        let (rel_pos_sum, densities_sum) =
            test_batch.partition(breakpoints_sum, agg_sum_fn)?;

        assert_eq!(rel_pos_sum.len(), 2);
        assert_eq!(densities_sum.len(), 2);

        // Calculate expected densities (sums) for test_batch
        // Segment 1: indices 0..2 (exclusive end). Densities: [0.5, 0.33333334]
        let first_half_density_sum: f64 = test_batch.slice(0, mid_pos_index) // Use slice(0, mid_pos_index) for the first segment
            .density()
            .drop_nulls()
            .iter()
            .filter_map(|v| v)
            .map(|v| v as f64)
            .sum();
        // Segment 2: indices 2..5 (exclusive end). Densities: [0.75, 1.0, 0.5]
        let second_half_density_sum: f64 = test_batch.slice(mid_pos_index as i64, batch_size - mid_pos_index) // Use slice from mid_pos_index
            .density()
            .drop_nulls()
            .iter()
            .filter_map(|v| v)
            .map(|v| v as f64)
            .sum();

        assert_approx_eq!(densities_sum[0], first_half_density_sum, 1e-6);
        assert_approx_eq!(densities_sum[1], second_half_density_sum, 1e-6);

        // Check relative positions for test_batch
        let first_pos = test_batch.first_pos().unwrap() as f64;
        let last_pos = test_batch.last_pos().unwrap() as f64;
        let total_length = last_pos - first_pos + 1.0; // 15 - 3 + 1 = 13
        let breakpoint_pos_value =
            test_batch.position().get(mid_pos_index).unwrap() as f64; // pos at index 2 is 9
        let expected_rel_breakpoint = (breakpoint_pos_value - first_pos) / total_length; // (9 - 3) / 13 = 6 / 13
        assert_approx_eq!(rel_pos_sum[0], expected_rel_breakpoint, 1e-6);
        assert_approx_eq!(rel_pos_sum[1], 1.0, 1e-6);
        Ok(())
    }

    #[rstest]
    fn test_partition_empty_batch() -> anyhow::Result<()> {
        let empty_batch = BsxBatch::empty(None);
        let agg_fn = Box::new(|slice: &[f32]| {
            slice.iter().sum::<f32>() as f64 / slice.len() as f64
        });
        let breakpoints = vec![1, 2];

        let (rel_pos, densities) = empty_batch.partition(breakpoints, agg_fn)?;
        assert!(rel_pos.is_empty());
        assert!(densities.is_empty());

        Ok(())
    }

    #[rstest]
    fn test_discretise_basic_case(real_batch_slice_1k: BsxBatch) -> anyhow::Result<()> {
        use assert_approx_eq::assert_approx_eq;

        let batch_size = real_batch_slice_1k.len();
        let agg_mean = AggMethod::Mean;
        let n_fragments_basic = 10;

        if batch_size >= n_fragments_basic {
            let (rel_pos, densities) =
                real_batch_slice_1k.discretise(n_fragments_basic, agg_mean.clone())?;
            assert_eq!(rel_pos.len(), n_fragments_basic);
            assert_eq!(densities.len(), n_fragments_basic);
            assert_approx_eq!(rel_pos.last().copied().unwrap(), 1.0, 1e-6);

            // Check relative positions are increasing
            let mut last_pos = 0.0;
            for pos in &rel_pos {
                assert!(*pos >= last_pos); // Should be strictly increasing unless multiple breakpoints are the
                                           // same index
                last_pos = *pos;
            }
        }
        else {
            eprintln!(
                "Skipping test_discretise_basic_case: Batch size {} is less than {} \
                 fragments",
                batch_size, n_fragments_basic
            );
        }

        Ok(())
    }

    // FIXME
    #[rstest]
    fn test_discretise_n_fragments_equals_len(
        real_batch_slice_1k: BsxBatch
    ) -> anyhow::Result<()> {
        use assert_approx_eq::assert_approx_eq;
        // No need for statrs::statistics::Data here, just using the get_fn
        // use statrs::statistics::Data;

        let batch = &real_batch_slice_1k; // Use a shorter name
        let batch_size = batch.len();
        let agg_method = AggMethod::Mean; // Use Mean for the test as in the original
        let agg_fn = agg_method.get_fn(); // Get the actual aggregation function
        let n_fragments_eq_len = batch_size;

        // The test only makes sense for a non-empty batch with more than one element,
        // as n_fragments=1 case is handled separately and tested elsewhere.
        if batch_size <= 1 {
            eprintln!(
                "Skipping test_discretise_n_fragments_equals_len: Batch size {} <= 1",
                batch_size
            );
            return Ok(()); // Skip if batch size is too small
        }
        // Ensure n_fragments is not 1 (handled in another test)
        if n_fragments_eq_len == 1 {
            eprintln!(
                "Skipping test_discretise_n_fragments_equals_len: n_fragments is 1"
            );
            return Ok(()); // Skip if n_fragments is 1
        }
        // Ensure n_fragments matches batch_size as per the test name's intent
        if n_fragments_eq_len != batch_size {
            eprintln!(
                "Skipping test_discretise_n_fragments_equals_len: n_fragments {} != \
                 batch_size {}",
                n_fragments_eq_len, batch_size
            );
            return Ok(()); // Skip if n_fragments doesn't match batch_size
        }

        // Get actual data for calculation
        let positions_series = batch.position();
        let positions: Vec<PosType> = positions_series.into_no_null_iter().collect();
        let densities_f32: Vec<f32> =
            batch.density().drop_nulls().into_no_null_iter().collect();
        // Check if densities_f32 has the same size as positions - it should because
        // filter_coverage_gt(0) was used
        assert_eq!(positions.len(), densities_f32.len());
        let size = positions.len();

        let start = batch.first_pos().unwrap();
        let end = batch.last_pos().unwrap();
        let genomic_length = (end + 1 - start) as f64;

        let fragment_genomic_length = genomic_length / n_fragments_eq_len as f64;
        let mut target_positions: Vec<PosType> =
            Vec::with_capacity(n_fragments_eq_len - 1);
        for i in 1..n_fragments_eq_len {
            let target_pos_f64 = start as f64 + i as f64 * fragment_genomic_length;
            let target_pos = target_pos_f64.round() as PosType;
            target_positions.push(target_pos);
        }

        // Simulate breakpoint finding from discretise
        let mut expected_partition_breakpoints: Vec<usize> =
            Vec::with_capacity(n_fragments_eq_len - 1);
        let mut current_search_start_index = 0;
        for target_pos in target_positions {
            let search_slice = &positions[current_search_start_index..];
            let result_in_slice =
                search_slice.binary_search_by(|pos| pos.cmp(&target_pos));
            let breakpoint_index_in_slice = match result_in_slice {
                Ok(idx) => idx,
                Err(idx) => idx, /* Index where element *would* be inserted, i.e.,
                                  * first element >= target */
            };

            // Convert the index in the slice back to an index in the original
            // 'positions' vector.
            let breakpoint_index =
                current_search_start_index + breakpoint_index_in_slice;
            let final_index = std::cmp::min(breakpoint_index, size); // Ensure index is within bounds [0, size]

            expected_partition_breakpoints.push(final_index);
            current_search_start_index = final_index; // Start next search from
                                                      // here
        }
        debug_assert_eq!(expected_partition_breakpoints.len(), n_fragments_eq_len - 1); // This is what discretise passes to partition


        let mut expected_rel_pos: Vec<f64> = Vec::with_capacity(n_fragments_eq_len);
        let mut partition_internal_bps = expected_partition_breakpoints.clone();
        // This check should match partition's internal logic: if breakpoints.last() !=
        // &size OR breakpoints is empty, push size.
        if partition_internal_bps
            .last()
            .map(|p| p != &size)
            .unwrap_or(true)
        {
            partition_internal_bps.push(size);
        }
        // Now partition_internal_bps is [b_0, ..., b_{n-2}, size] (length n) (where n =
        // n_fragments_eq_len) partition calculates rel pos using indices 0..n-1
        // of THIS vector. The indices used are the *values* in
        // partition_internal_bps at indices 0..n-1.
        for &idx_value in partition_internal_bps.iter().take(n_fragments_eq_len - 1) {
            // iterates over values b_0 .. b_{n-2}
            // idx_value is an index *into* the original positions vector.
            let pos_value = positions[idx_value] as f64;
            expected_rel_pos.push((pos_value - start as f64) / genomic_length);
        }
        expected_rel_pos.push(1.0); // Always add 1.0 at the end


        let mut segmentation_indices = vec![0];
        segmentation_indices.extend_from_slice(&expected_partition_breakpoints);
        // partition's logic: if breakpoints.last() != &size OR breakpoints is empty,
        // push size. Here expected_partition_breakpoints is NOT empty (size > 1
        // implies n_fragments > 1, so n_fragments-1 >= 1). We need to check if
        // its last element is != size.
        if segmentation_indices
            .last()
            .map(|p| p != &size)
            .unwrap_or(true)
        {
            // Handles empty case too, though not relevant here
            segmentation_indices.push(size);
        }
        // Now segmentation_indices is [0, b_0, ..., b_{n-2}, size] (length n+1)

        let mut expected_densities: Vec<f64> = Vec::with_capacity(n_fragments_eq_len);
        for window in segmentation_indices.windows(2) {
            let start_i = window[0]; // start index of the segment (inclusive)
            let end_i = window[1]; // end index of the segment (exclusive)
                                   // Note: slice indexing is [start..end)
            let segment_densities = &densities_f32[start_i..end_i];
            expected_densities.push(agg_fn(segment_densities));
        }
        debug_assert_eq!(expected_densities.len(), n_fragments_eq_len); // Should have n_fragments densities


        // Call the actual method
        let (rel_pos_eq, densities_eq) =
            batch.discretise(n_fragments_eq_len, agg_method.clone())?;

        // Assert results match expected values
        assert_eq!(rel_pos_eq.len(), n_fragments_eq_len);
        assert_eq!(densities_eq.len(), n_fragments_eq_len);

        // Compare relative positions
        assert_eq!(rel_pos_eq.len(), expected_rel_pos.len());
        for (actual, expected) in rel_pos_eq.iter().zip(expected_rel_pos.iter()) {
            assert_approx_eq!(*actual, *expected, 1e-6);
        }

        // Compare densities
        assert_eq!(densities_eq.len(), expected_densities.len());
        for (actual, expected) in densities_eq.iter().zip(expected_densities.iter()) {
            // Handle potential NAN values (e.g., if empty slices occur and agg_fn
            // returns NAN)
            if actual.is_nan() && expected.is_nan() {
                // Both are NAN, this is ok
            }
            else {
                // If agg_fn panics on empty slices, this check won't be reached for
                // problematic cases.
                assert_approx_eq!(*actual, *expected, 1e-6);
            }
        }

        Ok(())
    }

    #[rstest]
    fn test_discretise_n_fragments_one(
        real_batch_slice_1k: BsxBatch
    ) -> anyhow::Result<()> {
        use assert_approx_eq::assert_approx_eq;

        let agg_mean = AggMethod::Mean;
        let n_fragments_one = 1;
        let (rel_pos_one, densities_one) =
            real_batch_slice_1k.discretise(n_fragments_one, agg_mean.clone())?;
        assert_eq!(rel_pos_one, vec![1.0]);
        assert_eq!(densities_one.len(), 1);
        let expected_density_one = real_batch_slice_1k
            .density()
            .drop_nulls()
            .mean()
            .unwrap_or(f64::NAN) as f64;
        assert_approx_eq!(densities_one[0], expected_density_one, 1e-6);

        Ok(())
    }

    #[rstest]
    fn test_discretise_error_n_fragments_zero(
        real_batch_slice_1k: BsxBatch
    ) -> anyhow::Result<()> {
        let agg_mean = AggMethod::Mean;
        let n_fragments_zero = 0;
        let result_zero =
            real_batch_slice_1k.discretise(n_fragments_zero, agg_mean.clone());
        assert!(result_zero.is_err());
        assert_eq!(
            result_zero.unwrap_err().to_string(),
            "Cannot partition batch into 0 fragments"
        );
        Ok(())
    }

    #[rstest]
    #[case(AggMethod::Mean)]
    #[case(AggMethod::Median)]
    #[case(AggMethod::Max)]
    #[case(AggMethod::Min)]
    #[case(AggMethod::GeometricMean)]
    fn test_discretise_different_agg_methods(
        test_batch: BsxBatch,
        #[case] agg_method: AggMethod,
    ) -> anyhow::Result<()> {
        use assert_approx_eq::assert_approx_eq;
        use statrs::statistics::Data;

        let n_fragments_agg = 2;

        let densities_f32: Vec<f32> = test_batch
            .density()
            .drop_nulls()
            .iter()
            .map(|v| v.unwrap())
            .collect();
        let pos_values: Vec<u32> = test_batch.position().into_no_null_iter().collect();
        let start_pos = test_batch.first_pos().unwrap() as f64;
        let end_pos = test_batch.last_pos().unwrap() as f64;
        let genomic_length = end_pos + 1.0 - start_pos;

        let fragment_genomic_length = genomic_length / n_fragments_agg as f64; // 13 / 2 = 6.5
                                                                               // Target positions based on genomic length division
        let mut target_positions: Vec<PosType> =
            Vec::with_capacity(n_fragments_agg - 1);
        for i in 1..n_fragments_agg {
            let target_pos_f64 = start_pos + i as f64 * fragment_genomic_length; // 3.0 + 1 * 6.5 = 9.5
            target_positions.push(target_pos_f64.round() as PosType); // round(9.5) = 10
        }
        let expected_breakpoints_for_partition: Vec<usize> = vec![3]; // These are the indices passed to partition
        let expected_segment1_indices = 0..3; // indices 0, 1, 2 -> positions 3, 5, 9
        let expected_segment2_indices = 3..test_batch.len(); // indices 3, 4 -> positions 12, 15

        let segment1_densities: Vec<f64> = densities_f32
            [expected_segment1_indices.clone()]
        .iter()
        .map(|&v| v as f64)
        .collect();
        let segment2_densities: Vec<f64> = densities_f32
            [expected_segment2_indices.clone()]
        .iter()
        .map(|&v| v as f64)
        .collect();

        let expected_rel_pos: Vec<f64> = {
            let breakpoint_pos_value =
                pos_values[expected_breakpoints_for_partition[0]] as f64; // pos at index 3 is 12
            vec![(breakpoint_pos_value - start_pos) / genomic_length, 1.0] // [(12 - 3) / 13.0, 1.0] = [9.0/13.0, 1.0]
        };


        let (rel_pos_agg, densities_agg) =
            test_batch.discretise(n_fragments_agg, agg_method.clone())?;

        // discretise returns n_fragments relative positions and n_fragments densities
        assert_eq!(rel_pos_agg.len(), n_fragments_agg); // 2
        assert_eq!(densities_agg.len(), n_fragments_agg); // 2

        // Check relative positions are correct (same for all agg methods)
        assert_approx_eq!(rel_pos_agg[0], expected_rel_pos[0], 1e-6); // 9.0/13.0
        assert_approx_eq!(rel_pos_agg[1], expected_rel_pos[1], 1e-6); // 1.0

        // Check densities based on agg method for segments [0..3) and [3..5)
        match agg_method {
            AggMethod::Mean => {
                let expected1 = segment1_densities.iter().sum::<f64>()
                    / segment1_densities.len() as f64;
                let expected2 = segment2_densities.iter().sum::<f64>()
                    / segment2_densities.len() as f64;
                assert_approx_eq!(densities_agg[0], expected1, 1e-6);
                assert_approx_eq!(densities_agg[1], expected2, 1e-6);
            },
            AggMethod::Sum => {
                let expected1 = segment1_densities.iter().sum::<f64>();
                let expected2 = segment2_densities.iter().sum::<f64>();
                assert_approx_eq!(densities_agg[0], expected1, 1e-6);
                assert_approx_eq!(densities_agg[1], expected2, 1e-6);
            },
            AggMethod::GeometricMean => {
                let expected1 = Data::new(segment1_densities.clone())
                    .iter()
                    .geometric_mean();
                let expected2 = Data::new(segment2_densities.clone())
                    .iter()
                    .geometric_mean();
                assert_approx_eq!(densities_agg[0], expected1, 1e-6);
                assert_approx_eq!(densities_agg[1], expected2, 1e-6);
            },
            AggMethod::Median => {
                let mut segment1_densities_sorted = segment1_densities.clone();
                segment1_densities_sorted
                    .sort_by(|a, b| a.partial_cmp(b).unwrap_or(Ordering::Equal));
                let expected1 = Data::new(segment1_densities_sorted).percentile(50);

                let mut segment2_densities_sorted = segment2_densities.clone();
                segment2_densities_sorted
                    .sort_by(|a, b| a.partial_cmp(b).unwrap_or(Ordering::Equal));
                let expected2 = Data::new(segment2_densities_sorted).percentile(50);
                assert_approx_eq!(densities_agg[0], expected1, 1e-6);
                assert_approx_eq!(densities_agg[1], expected2, 1e-6);
            },
            AggMethod::Max => {
                let expected1 = segment1_densities
                    .iter()
                    .cloned()
                    .max_by(|a, b| a.partial_cmp(b).unwrap_or(Ordering::Equal))
                    .unwrap_or(f64::NAN);
                let expected2 = segment2_densities
                    .iter()
                    .cloned()
                    .max_by(|a, b| a.partial_cmp(b).unwrap_or(Ordering::Equal))
                    .unwrap_or(f64::NAN);
                assert_approx_eq!(densities_agg[0], expected1, 1e-6);
                assert_approx_eq!(densities_agg[1], expected2, 1e-6);
            },
            AggMethod::Min => {
                let expected1 = segment1_densities
                    .iter()
                    .cloned()
                    .min_by(|a, b| a.partial_cmp(b).unwrap_or(Ordering::Equal))
                    .unwrap_or(f64::NAN);
                let expected2 = segment2_densities
                    .iter()
                    .cloned()
                    .min_by(|a, b| a.partial_cmp(b).unwrap_or(Ordering::Equal))
                    .unwrap_or(f64::NAN);
                assert_approx_eq!(densities_agg[0], expected1, 1e-6);
                assert_approx_eq!(densities_agg[1], expected2, 1e-6);
            },
        }

        Ok(())
    }

    // // Case 9: Error case: Density contains nulls
    // This test requires injecting nulls into the density column,
    // which involves modifying the DataFrame, potentially violating
    // the assumptions of other tests using the same fixture.
    // Re-enable and fix carefully if needed.
    // #[rstest]
    // fn test_discretise_error_density_nulls(test_batch: BsxBatch) ->
    // anyhow::Result<()> {     let agg_mean = AggMethod::Mean;
    //     let n_fragments = 2;

    //     let mut test_batch_with_nulls = test_batch.clone();
    //     // Modify the density column to have nulls
    //     let mut density_series = test_batch_with_nulls.density().clone();
    //     density_series.set(0, AnyValue::Null)?; // Set first value to null
    //     test_batch_with_nulls.data.with_column(density_series.into_column())?;

    //     let result_with_nulls = test_batch_with_nulls.discretise(n_fragments,
    // agg_mean.clone());     assert!(result_with_nulls.is_err());
    //     assert_eq!(result_with_nulls.unwrap_err().to_string(), "Density contains
    // nulls");     Ok(())
    // }

    #[cfg(feature = "tools")]
    #[test]
    fn test_shrink() -> anyhow::Result<()> {
        use std::fs::File;
        use std::path::PathBuf;

        use crate::io::bsx::BsxFileReader;

        let reader: BsxFileReader = BsxFileReader::try_new(File::open(
            PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data/report.bsx"),
        )?)?;

        for batch_res in reader.into_iter().take(20) {
            let batch = batch_res?;
            let filtered = batch
                .lazy()
                .filter_context(Context::CG)
                .filter_coverage_gt(5)
                .collect()?;
            if filtered.is_empty() {
                continue;
            }

            let (rel_pos, density) =
                filtered.shrink(SegmentAlgorithm::Pelt(None, 10))?;

            assert_eq!(rel_pos.len(), density.len());
            assert!(rel_pos.len() > 0);
            assert!(rel_pos.iter().all(|v| *v >= 0.0 && *v <= 1.0));
            assert!(density.iter().all(|v| *v >= 0.0 && *v <= 1.0));
        }

        Ok(())
    }

    #[rstest]
    #[case::bismark(ReportType::Bismark,
        df!(
            BsxCol::Chr.as_str() => ["chr1", "chr1"],
            BsxCol::Position.as_str() => [100u32, 150],
            "strand" => ["+", "-"], // Bismark uses "+", "-"
            "count_m" => [5u16, 10],
            "count_um" => [5u16, 10],
            "context" => ["CG", "CHG"], // Bismark uses String context
            "trinuc" => ["CGA", "CAG"] // Bismark requires trinuc column
        )
        .unwrap()
    )]
    #[case::cgmap(ReportType::CgMap,
        df!(
            BsxCol::Chr.as_str() => ["chr2", "chr2"],
            "nuc" => ["C", "G"], // CgMap uses nuc ('C' -> '+', 'G' -> '-')
            BsxCol::Position.as_str() => [200u32, 250],
            "context" => ["CG", "CHH"], // CgMap uses String context
            "dinuc" => ["CG", "CA"], // CgMap requires dinuc, not pattern
            BsxCol::Density.as_str() => [0.8f64, 0.4], // CgMap uses density column
            BsxCol::CountM.as_str() => [8u16, 2],
            BsxCol::CountTotal.as_str() => [10u16, 5]
        )
        .unwrap()
    )]
    #[case::coverage(ReportType::Coverage,
        df!(
            BsxCol::Chr.as_str() => ["chr3", "chr3"],
            "start" => [300u32, 350], // Coverage uses start/end
            "end" => [301u32, 351],
            BsxCol::Density.as_str() => [0.75f64, 0.4], // Coverage uses density column
            "count_m" => [15u16, 20],
            "count_um" => [5u16, 30] // Used to calculate count_total
        )
        .unwrap()
    )]
    #[case::bedgraph(ReportType::BedGraph,
        df!(
            BsxCol::Chr.as_str() => ["chr4", "chr4"],
            "start" => [400u32, 450], // BedGraph uses start/end
            "end" => [401u32, 451],
            BsxCol::Density.as_str() => [0.75f64, 0.8] // Only density provided
        )
        .unwrap()
    )]
    fn test_build_bsx_from_report_type(
        #[case] report_type: ReportType,
        #[case] input_df: DataFrame,
    ) -> anyhow::Result<()> {
        let builder = BsxBatchBuilder::no_checks();
        let batch: BsxBatch = builder.build_from_report(input_df, report_type)?;
        let df = batch.into_inner();

        // Check common columns and types
        assert_eq!(
            df.column(BsxCol::Chr.as_str())?.dtype(),
            &BsxCol::Chr.dtype()
        );
        assert_eq!(
            df.column(BsxCol::Position.as_str())?.dtype(),
            &BsxCol::Position.dtype()
        );
        assert_eq!(
            df.column(BsxCol::Strand.as_str())?.dtype(),
            &BsxCol::Strand.dtype()
        );
        assert_eq!(
            df.column(BsxCol::Context.as_str())?.dtype(),
            &BsxCol::Context.dtype()
        );
        assert_eq!(
            df.column(BsxCol::CountM.as_str())?.dtype(),
            &BsxCol::CountM.dtype()
        );
        assert_eq!(
            df.column(BsxCol::CountTotal.as_str())?.dtype(),
            &BsxCol::CountTotal.dtype()
        );
        assert_eq!(
            df.column(BsxCol::Density.as_str())?.dtype(),
            &BsxCol::Density.dtype()
        );

        // Check specific transformations (example for CgMap strand)
        if report_type == ReportType::CgMap {
            // Input nuc "C", "G" -> should be strand true, false
            let expected_strand =
                Series::new(BsxCol::Strand.as_str().into(), [true, false]);
            assert!(df
                .column(BsxCol::Strand.as_str())?
                .equals(&expected_strand.into_column()));
        }
        // Check specific transformations (example for Bismark/Coverage
        // count_total calculation from m + um)
        if report_type == ReportType::Bismark || report_type == ReportType::Coverage {
            let count_m = df.column(BsxCol::CountM.as_str())?.u16()?;
            let count_total = df.column(BsxCol::CountTotal.as_str())?.u16()?;
            assert!(count_m
                .into_iter()
                .zip(count_total.into_iter())
                .all(|(m, t)| m.unwrap_or(0) <= t.unwrap_or(0)));
        }
        // Check specific transformations (example for BedGraph nulls)
        if report_type == ReportType::BedGraph {
            assert!(df.column(BsxCol::CountM.as_str())?.null_count() > 0);
            assert!(df.column(BsxCol::CountTotal.as_str())?.null_count() > 0);
            assert!(df.column(BsxCol::Context.as_str())?.null_count() > 0);
        }

        Ok(())
    }

    #[rstest]
    #[case::bedgraph(ReportType::BedGraph,
        df!(
            "chr" => &["chr1", "chr1", "chr1", "chr1", "chr1"],
            "start" => &[3u32, 5, 9, 12, 15],
            "end" => &[4u32, 6, 10, 13, 16],
            "density" => &[0.5f32, 0.33333334f32, 0.75f32, 1.0f32, 0.5f32]
        )
        .unwrap()
    )]
    #[case::coverage(ReportType::Coverage,
        df!(
            "chr" => &["chr1", "chr1", "chr1", "chr1", "chr1"],
            "start" => &[3u32, 5, 9, 12, 15],
            "end" => &[4u32, 6, 10, 13, 16],
            "density" => &[0.5f32, 0.33333334f32, 0.75f32, 1.0f32, 0.5f32],
            "count_m" => &[5u16, 10, 15, 10, 5],
            "count_um" => &[5u16, 20, 5, 0, 5],
        )
        .unwrap()
    )]
    #[case::bismark(ReportType::Bismark,
        df!(
            "chr" => &["chr1", "chr1", "chr1", "chr1", "chr1"],
            "position" => &[3u32, 5, 9, 12, 15],
            "strand" => &["+", "-", "+", "+", "-"],
            "count_m" => &[5u16, 10, 15, 10, 5],
            "count_um" => &[5u16, 20, 5, 0, 5],
            "context" => &["CG", "CHG", "CG", "CHH", "CHH"],
            "trinuc" => &["CG", "CHG", "CG", "CHH", "CHH"]
        )
        .unwrap()
    )]
    #[case::cgmap(ReportType::CgMap,
        df!(
            "chr" => &["chr1", "chr1", "chr1", "chr1", "chr1"],
            "nuc" => &["C", "G", "C", "C", "G"],
            "position" => &[3u32, 5, 9, 12, 15],
            "context" => &["CG", "CHG", "CG", "CHH", "CHH"],
            "dinuc" => &["CG", "CH", "CG", "CH", "CH"],
            "density" => &[0.5f32, 0.33333334f32, 0.75f32, 1.0f32, 0.5f32],
            "count_m" => &[5u16, 10, 15, 10, 5],
            "count_total" => &[10u16, 30, 20, 10, 10]
        )
        .unwrap()
    )]
    fn test_into_report(
        test_batch: BsxBatch,
        #[case] report_type: ReportType,
        #[case] expected_df: DataFrame,
    ) -> anyhow::Result<()> {
        let result_df = test_batch.clone().into_report(report_type)?;
        // Polars assert_eq! checks content and schema
        assert_eq!(result_df, expected_df);

        Ok(())
    }

    #[test]
    fn test_add_context_data() -> anyhow::Result<()> {
        // Original BsxBatch data
        let original_positions = vec![10, 20, 30];
        let original_strands = vec![false, false, false]; // F, R, F
        let original_contexts = vec![None, None, None]; // CG, CHG, CHH
        let original_count_m = vec![5, 10, 15];
        let original_count_total = vec![10, 20, 30];

        let original_chr = "chr_orig";
        let chr_dtype = Some(create_caregorical_dtype(vec![original_chr.into()]));

        let original_batch = BsxBatch::try_from_columns(
            original_chr,
            chr_dtype.clone(),
            original_positions,
            original_strands,
            original_contexts,
            original_count_m,
            original_count_total,
        )?;

        let context_positions = vec![10, 20, 25, 30];
        let context_strands = vec![
            Strand::Reverse,
            Strand::Reverse,
            Strand::Forward,
            Strand::Forward,
        ];
        let context_contexts =
            vec![Context::CG, Context::CHG, Context::CHH, Context::CG];

        let context_data = ContextData::new(
            context_positions.clone(),
            context_strands.clone(),
            context_contexts.clone(),
        );

        // Add context data to the original batch
        let result_batch = original_batch.add_context_data(context_data)?;
        let result_df = result_batch.into_inner();

        let expected_positions = vec![10u32, 20, 25, 30];
        // Strand and Context are from ContextData (converted to bool)
        let expected_strands = vec![Some(false), Some(false), Some(true), Some(true)]; // R -> false, F -> true
        let expected_contexts = vec![Some(true), Some(false), None, Some(true)]; // CHG -> false, CHH -> None, CG -> true
                                                                                 // Chr is filled from original batch (all rows get the same chr)
        let expected_chr_values =
            vec![original_chr, original_chr, original_chr, original_chr];
        let expected_count_m = vec![5u16, 10, 0, 15];
        let expected_count_total = vec![10u16, 20, 0, 30];
        let expected_density = vec![0.5f32, 0.5, DensityType::NAN, 0.5f32]; // Use DensityType::NAN

        // Create expected DataFrame using the actual BsxCol types and casting
        let expected_df = df!(
            BsxCol::Chr.as_str() => Series::new(BsxCol::Chr.as_str().into(), expected_chr_values).cast(&BsxCol::Chr.dtype())?,
            BsxCol::Position.as_str() => Series::new(BsxCol::Position.as_str().into(), expected_positions).cast(&BsxCol::Position.dtype())?,
            BsxCol::Strand.as_str() => Series::new(BsxCol::Strand.as_str().into(), expected_strands).cast(&BsxCol::Strand.dtype())?,
            BsxCol::Context.as_str() => Series::new(BsxCol::Context.as_str().into(), expected_contexts).cast(&BsxCol::Context.dtype())?,
            BsxCol::CountM.as_str() => Series::new(BsxCol::CountM.as_str().into(), expected_count_m).cast(&BsxCol::CountM.dtype())?,
            BsxCol::CountTotal.as_str() => Series::new(BsxCol::CountTotal.as_str().into(), expected_count_total).cast(&BsxCol::CountTotal.dtype())?,
            BsxCol::Density.as_str() => Series::new(BsxCol::Density.as_str().into(), expected_density).cast(&BsxCol::Density.dtype())?,
        )?;

        // Check equality (content and schema)
        assert_eq!(result_df, expected_df);

        Ok(())
    }

    #[rstest]
    fn test_column_getters(test_batch: BsxBatch) {
        assert_eq!(test_batch.len(), 5);

        // Test chr column
        let chr = test_batch.chr();
        assert!(chr.iter_str().all(|c| c.unwrap() == "chr1"));

        // Test position column
        let positions = test_batch.position();
        assert_eq!(positions.into_iter().collect::<Vec<_>>(), vec![
            Some(3),
            Some(5),
            Some(9),
            Some(12),
            Some(15)
        ]);

        // Test strand column
        let strands = test_batch.strand();
        assert_eq!(strands.into_iter().collect::<Vec<_>>(), vec![
            Some(true),
            Some(false),
            Some(true),
            Some(true),
            Some(false)
        ]);

        // Test context column
        let contexts = test_batch.context();
        assert_eq!(contexts.into_iter().collect::<Vec<_>>(), vec![
            Some(true),
            Some(false),
            Some(true),
            None,
            None
        ]);

        // Test count_m column
        let count_m = test_batch.count_m();
        assert_eq!(count_m.into_iter().collect::<Vec<_>>(), vec![
            Some(5),
            Some(10),
            Some(15),
            Some(10),
            Some(5)
        ]);

        // Test count_total column
        let count_total = test_batch.count_total();
        assert_eq!(count_total.into_iter().collect::<Vec<_>>(), vec![
            Some(10),
            Some(30),
            Some(20),
            Some(10),
            Some(10)
        ]);

        // Test density column
        let density = test_batch.density();
        assert_eq!(
            density
                .into_iter()
                .map(|v| v.map(|f| (f * 100.0).round() / 100.0))
                .collect::<Vec<_>>(),
            vec![Some(0.5), Some(0.33), Some(0.75), Some(1.0), Some(0.5)]
        );
    }

    #[rstest]
    fn test_split_at(test_batch: BsxBatch) {
        let (left, right) = test_batch.split_at(2);

        assert_eq!(left.len(), 2);
        assert_eq!(right.len(), 3);

        assert_eq!(left.position().into_iter().collect::<Vec<_>>(), vec![
            Some(3),
            Some(5)
        ]);
        assert_eq!(right.position().into_iter().collect::<Vec<_>>(), vec![
            Some(9),
            Some(12),
            Some(15)
        ]);
    }

    #[rstest]
    fn test_slice(test_batch: BsxBatch) {
        let sliced = test_batch.slice(1, 3);

        assert_eq!(sliced.len(), 3);
        assert_eq!(sliced.position().into_iter().collect::<Vec<_>>(), vec![
            Some(5),
            Some(9),
            Some(12)
        ]);
    }

    #[rstest]
    fn test_position_methods(test_batch: BsxBatch) {
        assert_eq!(test_batch.first_pos(), Some(3));
        assert_eq!(test_batch.last_pos(), Some(15));
        assert_eq!(test_batch.seqname(), Some("chr1"));

        let first_genomic_pos = test_batch.first_genomic_pos();
        assert_eq!(
            first_genomic_pos,
            Some(GenomicPosition::new("chr1".into(), 3))
        );

        let last_genomic_pos = test_batch.last_genomic_pos();
        assert_eq!(
            last_genomic_pos,
            Some(GenomicPosition::new("chr1".into(), 15))
        );

        let contig = test_batch.as_contig();
        assert_eq!(
            contig,
            Some(Contig::new("chr1".into(), 3, 15, Strand::None))
        );
    }

    #[rstest]
    fn test_get_methylation_stats(test_batch: BsxBatch) {
        use assert_approx_eq::assert_approx_eq;

        let meth_stats = test_batch.get_methylation_stats();

        // Calculate expected values based on test_batch data (densities: [0.5,
        // 0.333..., 0.75, 1.0, 0.5])
        let expected_mean =
            (0.5 + 10.0 / 30.0 + 15.0 / 20.0 + 10.0 / 10.0 + 5.0 / 10.0) / 5.0; // 0.616666668
                                                                                // Variance calculation from test_batch density values (using f64 for
                                                                                // calculation accuracy)
        let density_f64: Vec<f64> = test_batch
            .density()
            .drop_nulls()
            .into_no_null_iter()
            .map(|v| v as f64)
            .collect();
        let expected_var = statrs::statistics::Data::new(density_f64).iter().variance(); // Sample variance

        assert_approx_eq!(meth_stats.mean_methylation() as f64, expected_mean, 1e-6);
        // Check variance with a reasonable tolerance
        assert_approx_eq!(meth_stats.methylation_var() as f64, expected_var, 1e-6);
        // Note: coverage_dist, context_stats, and strand_stats are tested in
        // their own functions below.
    }

    #[rstest]
    fn test_get_coverage_dist(test_batch: BsxBatch) {
        let cov_dist = test_batch.get_coverage_dist();

        // Expected distribution from test_batch (CountM, CountTotal):
        // [(5, 10), (10, 30), (15, 20), (10, 10), (5, 10)]
        // Grouped: (5, 10) -> 2, (10, 30) -> 1, (15, 20) -> 1, (10, 10) -> 1
        let mut expected_cov_dist: HashMap<u16, u32> = HashMap::new();
        expected_cov_dist.insert(10, 3);
        expected_cov_dist.insert(20, 1);
        expected_cov_dist.insert(30, 1);

        assert_eq!(cov_dist.len(), 3);
        assert_eq!(cov_dist, expected_cov_dist);
    }

    #[rstest]
    fn test_get_context_stats(test_batch: BsxBatch) {
        use assert_approx_eq::assert_approx_eq;

        let context_stats = test_batch.get_context_stats();

        // Expected stats from test_batch (filtered for non-null density):
        // (Context, Density): [(CG, 0.5), (CHG, 0.333...), (CG, 0.75), (CHH, 1.0),
        // (CHH, 0.5)] CG: sum=0.5+0.75=1.25, count=2
        // CHG: sum=0.333..., count=1
        // CHH: sum=1.0+0.5=1.5, count=2
        let mut expected_context_stats: HashMap<Context, (f64, u32)> = HashMap::new();
        expected_context_stats.insert(Context::CG, (1.25, 2));
        expected_context_stats.insert(Context::CHG, (10.0 / 30.0 as f64, 1));
        expected_context_stats.insert(Context::CHH, (1.5, 2));


        assert_eq!(context_stats.len(), 3);
        assert_approx_eq!(
            context_stats
                .get(&Context::CG)
                .copied()
                .unwrap_or_default()
                .0 as f64,
            expected_context_stats
                .get(&Context::CG)
                .copied()
                .unwrap_or_default()
                .0,
            1e-6
        );
        assert_eq!(
            context_stats
                .get(&Context::CG)
                .copied()
                .unwrap_or_default()
                .1,
            expected_context_stats
                .get(&Context::CG)
                .copied()
                .unwrap_or_default()
                .1
        );

        assert_approx_eq!(
            context_stats
                .get(&Context::CHG)
                .copied()
                .unwrap_or_default()
                .0 as f64,
            expected_context_stats
                .get(&Context::CHG)
                .copied()
                .unwrap_or_default()
                .0,
            1e-6
        );
        assert_eq!(
            context_stats
                .get(&Context::CHG)
                .copied()
                .unwrap_or_default()
                .1,
            expected_context_stats
                .get(&Context::CHG)
                .copied()
                .unwrap_or_default()
                .1
        );

        assert_approx_eq!(
            context_stats
                .get(&Context::CHH)
                .copied()
                .unwrap_or_default()
                .0 as f64,
            expected_context_stats
                .get(&Context::CHH)
                .copied()
                .unwrap_or_default()
                .0,
            1e-6
        );
        assert_eq!(
            context_stats
                .get(&Context::CHH)
                .copied()
                .unwrap_or_default()
                .1,
            expected_context_stats
                .get(&Context::CHH)
                .copied()
                .unwrap_or_default()
                .1
        );
    }

    #[rstest]
    fn test_get_strand_stats(test_batch: BsxBatch) {
        use assert_approx_eq::assert_approx_eq;

        let strand_stats = test_batch.get_strand_stats();

        // Expected stats from test_batch (filtered for non-null density):
        // (Strand, Density): [(Forward, 0.5), (Reverse, 0.333...), (Forward, 0.75),
        // (Forward, 1.0), (Reverse, 0.5)] Forward: sum=0.5+0.75+1.0=2.25,
        // count=3 Reverse: sum=0.333...+0.5=0.833..., count=2
        let mut expected_strand_stats: HashMap<Strand, (f64, u32)> = HashMap::new();
        expected_strand_stats.insert(Strand::Forward, (2.25, 3));
        expected_strand_stats.insert(Strand::Reverse, (10.0 / 30.0 as f64 + 0.5, 2));

        assert_eq!(strand_stats.len(), 2);
        assert_approx_eq!(
            strand_stats
                .get(&Strand::Forward)
                .copied()
                .unwrap_or_default()
                .0 as f64,
            expected_strand_stats
                .get(&Strand::Forward)
                .copied()
                .unwrap_or_default()
                .0,
            1e-6
        );
        assert_eq!(
            strand_stats
                .get(&Strand::Forward)
                .copied()
                .unwrap_or_default()
                .1,
            expected_strand_stats
                .get(&Strand::Forward)
                .copied()
                .unwrap_or_default()
                .1
        );

        assert_approx_eq!(
            strand_stats
                .get(&Strand::Reverse)
                .copied()
                .unwrap_or_default()
                .0 as f64,
            expected_strand_stats
                .get(&Strand::Reverse)
                .copied()
                .unwrap_or_default()
                .0,
            1e-6
        );
        assert_eq!(
            strand_stats
                .get(&Strand::Reverse)
                .copied()
                .unwrap_or_default()
                .1,
            expected_strand_stats
                .get(&Strand::Reverse)
                .copied()
                .unwrap_or_default()
                .1
        );
    }

    #[rstest]
    fn test_as_binom(test_batch: BsxBatch) {
        let binom_batch = test_batch.as_binom(0.5, 0.05).unwrap();
        assert_eq!(binom_batch.len(), 5);

        let count_m_series = binom_batch.count_m();
        let count_total_series = binom_batch.count_total();
        let density_series = binom_batch.density();

        // Check types remain the same
        assert_eq!(count_m_series.dtype(), &DataType::UInt16);
        assert_eq!(count_total_series.dtype(), &DataType::UInt16);
        assert_eq!(density_series.dtype(), &DataType::Float32);

        // Define expected values based on binomial test results
        // (5/10, p=0.5) -> p-value approx 1.0 > 0.05 -> Not significant -> (0, 1, 0.0)
        // (10/30, p=0.5) -> p-value approx 0.050174 > 0.05 -> Not significant -> (1, 1,
        // 1.0) (15/20, p=0.5) -> p-value approx 0.01478577 < 0.05 ->
        // Significant -> (1, 1, 1.0) (10/10, p=0.5) -> p-value approx 0.0 <
        // 0.05 -> Significant -> (1, 1, 1.0) (5/10, p=0.5) -> p-value approx
        // 0.24609375 > 0.05 -> Not significant -> (0, 1, 0.0)
        let expected_count_m = vec![0u16, 0, 1, 1, 0];
        let expected_count_total = vec![1u16, 1, 1, 1, 1];
        let expected_density = vec![0.0f32, 0.0, 1.0, 1.0, 0.0];

        // Assert that the resulting series match the expected ones
        assert_eq!(
            count_m_series.to_vec_null_aware().left().unwrap(),
            expected_count_m
        );
        assert_eq!(
            count_total_series.to_vec_null_aware().left().unwrap(),
            expected_count_total
        );
        assert_eq!(
            density_series.to_vec_null_aware().left().unwrap(),
            expected_density
        );
    }

    #[rstest]
    fn test_is_empty(test_batch: BsxBatch) {
        assert!(!test_batch.is_empty());
        assert!(BsxBatch::empty(None).is_empty());
    }

    #[rstest]
    #[case::no_chr(None, None)]
    #[case::both_chr(Some(get_categorical_dtype(vec!["chr1".into()])), Some(get_categorical_dtype(vec!["chr1".into()])))]
    #[should_panic]
    #[case::different_types(None, Some(get_categorical_dtype(vec!["chr1".into()])))]
    fn test_can_extend(
        #[case] first_dtype: Option<DataType>,
        #[case] second_dtype: Option<DataType>,
    ) {
        let batch1 = BsxBatch::empty(first_dtype.as_ref());
        let batch2 = BsxBatch::try_from_columns(
            "chr1",
            second_dtype,
            vec![1, 2, 3],
            vec![true, false, true],
            vec![Some(true), Some(false), None],
            vec![1, 2, 3],
            vec![3, 6, 9],
        )
        .unwrap();

        assert!(matches!(
            batch1.column(BsxCol::Chr).dtype(),
            DataType::Categorical(_, _) | DataType::Enum(_, _)
        ));
        assert!(matches!(
            batch2.column(BsxCol::Chr).dtype(),
            DataType::Categorical(_, _) | DataType::Enum(_, _)
        ));

        let vstack = batch1.data().vstack(&batch2.data()).unwrap();
        assert!(matches!(
            vstack.column(BsxCol::Chr.as_str()).unwrap().dtype(),
            DataType::Categorical(_, _) | DataType::Enum(_, _)
        ));
        assert!(vstack
            .column(BsxCol::Chr.as_str())
            .unwrap()
            .categorical()
            .unwrap()
            .get_rev_map()
            .find("chr1")
            .is_some());
    }
}

mod builder_tests {
    use polars::series::IsSorted;

    use super::builder::build::*;
    use super::*;

    #[fixture]
    fn test_df() -> DataFrame {
        df!(
            BsxCol::Chr.as_str() => ["chr1", "chr1", "chr1"],
            BsxCol::Position.as_str() => [100u32, 150, 200],
            BsxCol::Strand.as_str() => [true, false, true],
            BsxCol::Context.as_str() => [true, false, true],
            BsxCol::CountM.as_str() => [5u16, 10, 15],
            BsxCol::CountTotal.as_str() => [10u16, 20, 30],
            BsxCol::Density.as_str() => [0.5f32, 0.5, 0.5]
        )
        .unwrap()
    }

    #[fixture]
    fn no_cols_df() -> DataFrame {
        DataFrame::empty()
    }

    #[rstest]
    #[case::pos_duplicates({
        df!(
            BsxCol::Chr.as_str() => ["chr1", "chr1", "chr1"],
            BsxCol::Position.as_str() => [100u32, 100, 200], // Duplicate at position 100
        ).unwrap()
    }, check_no_duplicates)]
    #[case::unsorted({
        df!(
            BsxCol::Position.as_str() => [200u32, 100, 150], // Unsorted positions
        ).unwrap()
    }, check_sorted)]
    #[case::nulls_chr({
        df!(
            BsxCol::Chr.as_str() => [Some("chr1"), Some("chr1"), None],
            BsxCol::Position.as_str() => [100u32, 150, 200],
        ).unwrap()
    }, check_no_nulls)]
    #[case::nulls_pos({
        df!(
            BsxCol::Chr.as_str() => ["chr1", "chr1", "chr1"],
            BsxCol::Position.as_str() => [Some(100), Some(150), None],
        ).unwrap()
    }, check_no_nulls)]
    #[case::multi_chr({
        df!(
            BsxCol::Chr.as_str() => ["chr1", "chr2", "chr3"],
            BsxCol::Position.as_str() => [100u32, 150, 200],
        ).unwrap()
    }, check_single_chr)]
    fn test_checks(
        test_df: DataFrame,
        no_cols_df: DataFrame,
        #[case] should_not_pass_df: DataFrame,
        #[case] check_fn: fn(&DataFrame) -> PolarsResult<bool>,
    ) {
        assert_eq!(check_fn(&test_df).unwrap(), true);
        assert_eq!(check_fn(&should_not_pass_df).unwrap(), false);
        assert!(check_fn(&no_cols_df).is_err());
    }

    #[rstest]
    fn test_set_flags(mut test_df: DataFrame) {
        assert!(set_flags(&mut test_df).is_ok());

        // Verify the sorting flag was set
        let pos_col = test_df.column(BsxCol::Position.as_str()).unwrap();
        assert_eq!(pos_col.is_sorted_flag(), IsSorted::Ascending);
    }

    #[test]
    fn test_conversion_expressions() {
        // Test nuc_to_bool_expr
        let df = df!("nuc" => ["C", "G", "A"]).unwrap();
        let res = df
            .lazy()
            .with_column(nuc_to_bool_expr().alias("result"))
            .collect()
            .unwrap();
        let results = res.column("result").unwrap().bool().unwrap();
        assert_eq!(results.get(0), Some(true));
        assert_eq!(results.get(1), Some(false));
        assert_eq!(results.get(2), None);

        // Test strand_to_bool_expr
        let df = df!("strand" => ["+", "-", "."]).unwrap();
        let res = df
            .lazy()
            .with_column(strand_to_bool_expr().alias("result"))
            .collect()
            .unwrap();
        let results = res.column("result").unwrap().bool().unwrap();
        assert_eq!(results.get(0), Some(true));
        assert_eq!(results.get(1), Some(false));
        assert_eq!(results.get(2), None);

        // Test context_to_bool_expr
        let df = df!("context" => ["CG", "CHG", "CHH"]).unwrap();
        let res = df
            .lazy()
            .with_column(context_to_bool_expr().alias("result"))
            .collect()
            .unwrap();
        let results = res.column("result").unwrap().bool().unwrap();
        assert_eq!(results.get(0), Some(true));
        assert_eq!(results.get(1), Some(false));
        assert_eq!(results.get(2), None);
    }

    #[test]
    fn test_count_total_and_density() {
        let df = df!(
            "count_m" => [5i64, 10],
            "count_um" => [5i64, 10]
        )
        .unwrap();

        // Test count_total_col_expr
        let res = df
            .clone()
            .lazy()
            .with_column(count_total_col_expr().alias("count_total"))
            .collect()
            .unwrap();
        let totals = res.column("count_total").unwrap().i64().unwrap();
        assert_eq!(totals.get(0), Some(10));
        assert_eq!(totals.get(1), Some(20));

        // Test density_col_expr
        let res = df
            .lazy()
            .with_column(count_total_col_expr().alias("count_total"))
            .with_column(density_col_expr().alias("density"))
            .collect()
            .unwrap();
        let densities = res.column("density").unwrap().f64().unwrap();
        assert_eq!(densities.get(0), Some(0.5));
        assert_eq!(densities.get(1), Some(0.5));
    }

    fn create_bismark_df() -> DataFrame {
        df!(
            "chr" => ["chr1", "chr1"],
            "position" => [100u32, 150],
            "strand" => ["+", "-"], // Decoded uses "+", "-"
            "count_m" => [5u32, 10],
            "count_um" => [5u32, 10],
            "context" => ["CG", "CHG"],
            "trinuc" => ["CGA", "CAG"] // Bismark requires trinuc column
        )
        .unwrap()
    }

    fn create_cgmap_df() -> DataFrame {
        df!(
            "chr" => ["chr2", "chr2"],
            "nuc" => ["C", "G"], // Determines strand ('C' -> '+', 'G' -> '-')
            "position" => [200u32, 250],
            "context" => ["CG", "CHH"],
            "dinuc" => ["CG", "CA"], // CgMap requires dinuc, not pattern
            "density" => [0.8f64, 0.4], // CgMap has density column
            "count_m" => [8u32, 2],
            "count_total" => [10u32, 5]
        )
        .unwrap()
    }

    fn create_coverage_df() -> DataFrame {
        df!(
            "chr" => ["chr3", "chr3"],
            "start" => [300u32, 350], // Becomes position
            "end" => [301u32, 351],
            "density" => [0.75f64, 0.4], // Coverage has density column
            "count_m" => [15u32, 20],
            "count_um" => [5u32, 30] // Used to calculate count_total
        )
        .unwrap()
    }

    fn create_bedgraph_df() -> DataFrame {
        df!(
            "chr" => ["chr4", "chr4"],
            "start" => [400u32, 450], // Becomes position
            "end" => [401u32, 451],
            "density" => [0.75f64, 0.8] // Only density provided
        )
        .unwrap()
    }

    #[rstest]
    #[case::bismark(ReportType::Bismark, create_bismark_df())]
    #[case::cgmap(ReportType::CgMap, create_cgmap_df())]
    #[case::coverage(ReportType::Coverage, create_coverage_df())]
    #[case::bedgraph(ReportType::BedGraph, create_bedgraph_df())]
    fn test_build_decoded_from_report_type(
        #[case] report_type: ReportType,
        #[case] input_df: DataFrame,
    ) -> anyhow::Result<()> {
        let batch =
            BsxBatchBuilder::all_checks().build_from_report(input_df, report_type)?;
        let _data = batch.data();
        Ok(())
    }
}

mod lazybatch_tests {
    use rstest::{
        fixture,
        rstest,
    };

    use super::*;

    #[fixture]
    fn test_lazybatch() -> LazyBsxBatch {
        batch_tests::test_batch().into()
    }

    #[rstest]
    #[case::filter_pos_lt(|lb: LazyBsxBatch| lb.filter_pos_lt(10), 3, vec![3, 5, 9])]
    #[case::filter_pos_gt(|lb: LazyBsxBatch| lb.filter_pos_gt(10), 2, vec![12, 15])]
    #[case::filter_coverage_gt(|lb: LazyBsxBatch| lb.filter_coverage_gt(10), 2, vec![5, 9])]
    #[case::filter_strand_true(|lb: LazyBsxBatch| lb.filter_strand(Strand::Forward), 3, vec![3, 9, 12])]
    #[case::filter_strand_false(|lb: LazyBsxBatch| lb.filter_strand(Strand::Reverse), 2, vec![5, 15])]
    #[case::filter_context_true(|lb: LazyBsxBatch| lb.filter_context(Context::CG), 2, vec![3, 9])]
    #[case::filter_context_false(|lb: LazyBsxBatch| lb.filter_context(Context::CHG), 1, vec![5])]
    fn test_ops(
        test_lazybatch: LazyBsxBatch,
        #[case] test_fn: fn(LazyBsxBatch) -> LazyBsxBatch,
        #[case] expected_len: usize,
        #[case] expected_pos: Vec<u32>,
    ) -> PolarsResult<()> {
        let filtered_batch = test_fn(test_lazybatch).collect()?;
        assert_eq!(filtered_batch.len(), expected_len);
        let pos = filtered_batch
            .data()
            .column(BsxCol::Position.as_str())?
            .u32()?;
        assert_eq!(pos.iter().flatten().collect::<Vec<_>>(), expected_pos);

        Ok(())
    }
}

mod other_tests {
    use polars::df;
    use polars::prelude::*;

    use super::*;

    // --- Helper Functions for Aggregation ---

    fn sum_agg(cols: Vec<&Column>) -> Column {
        if cols.is_empty() {
            panic!("sum_agg received empty input");
        }

        let series_vec: Vec<Series> = cols
            .into_iter()
            .map(|c| c.to_owned().as_materialized_series().clone())
            .collect();

        let initial_series = series_vec[0].clone();

        let sum_series = series_vec
            .into_iter()
            .skip(1) // Skip the first one as it's the initial value
            .fold(initial_series, |acc, s| {
                (&acc + &s).expect("Failed to sum Series in sum_agg") // Add error handling if needed
            });

        sum_series.into() // Convert the final Series back to a Column
    }

    fn mean_agg(cols: Vec<&Column>) -> Column {
        if cols.is_empty() {
            return Series::new_empty("mean_agg_empty".into(), &DataType::Float64)
                .into();
        }

        let series_vec: PolarsResult<Vec<Series>> = cols
            .into_iter()
            .map(|c| {
                c.to_owned()
                    .as_materialized_series()
                    .cast(&DataType::Float64)
            })
            .collect();

        let series_vec = match series_vec {
            Ok(vec) => vec,
            Err(e) => {
                panic!("Failed to cast column to Float64 in mean_agg: {}", e);
            },
        };

        let n = series_vec.len() as f64; // Use f64 for division

        let initial_series = series_vec[0].clone();

        let total_sum = series_vec
                .into_iter()
                .skip(1) // Skip the first one as it's the initial value
                .fold(initial_series, |acc, s| {
                    // Perform addition between Series
                    // This should work now as all series are Float64
                    (&acc + &s).expect("Failed to sum Series in mean_agg")
                });

        // Divide the sum by the number of columns
        let mean_series = &total_sum / n;
        mean_series.into() // Convert the final Series back to a Column
    }

    // --- Helper Functions to Create Test Batches ---

    #[allow(unsafe_code)]
    fn create_decoded_batch_1() -> BsxBatch {
        let df = df!(
            BsxColumns::Chr.as_str() => &["chr1", "chr1"],
            BsxColumns::Position.as_str() => &[10u64, 20],
            BsxColumns::Strand.as_str() => &["+", "-"],
            BsxColumns::Context.as_str() => &["CG", "CHG"],
            BsxColumns::CountM.as_str() => &[2u32, 4],
            BsxColumns::CountTotal.as_str() => &[10u32, 8],
            BsxColumns::Density.as_str() => &[0.2f64, 0.5],
        )
        .unwrap();
        unsafe { BsxBatch::new_unchecked(df) }
    }

    #[allow(unsafe_code)]
    fn create_decoded_batch_2() -> BsxBatch {
        let df = df!(
            BsxColumns::Chr.as_str() => &["chr1", "chr1"],
            BsxColumns::Position.as_str() => &[10u64, 20],
            BsxColumns::Strand.as_str() => &["+", "-"],
            BsxColumns::Context.as_str() => &["CG", "CHG"],
            BsxColumns::CountM.as_str() => &[3u32, 6],
            BsxColumns::CountTotal.as_str() => &[10u32, 12],
            BsxColumns::Density.as_str() => &[0.3f64, 0.5],
        )
        .unwrap();
        unsafe { BsxBatch::new_unchecked(df) }
    }

    // --- Test Cases ---

    #[test]
    fn test_merge_replicates_decoded() -> anyhow::Result<()> {
        let batch1 = create_decoded_batch_1();
        let batch2 = create_decoded_batch_2();
        let batches = vec![batch1, batch2];

        let merged_batch: BsxBatch =
            merge_replicates(batches, Box::new(sum_agg), Box::new(mean_agg))?;

        // Expected results
        let expected_df = df!(
            BsxColumns::Chr.as_str() => &["chr1", "chr1"],
            BsxColumns::Position.as_str() => &[10u64, 20],
            BsxColumns::Strand.as_str() => &["+", "-"],
            BsxColumns::Context.as_str() => &["CG", "CHG"],
            BsxColumns::CountM.as_str() => &[5u32, 10],
            BsxColumns::CountTotal.as_str() => &[20u32, 20],
            BsxColumns::Density.as_str() => &[0.25f64, 0.5],
        )?;

        // Compare DataFrames
        assert_eq!(merged_batch.data(), &expected_df);
        // Ensure the type is correct
        assert_eq!(merged_batch.count_m().dtype(), &BsxColumns::CountM.dtype());
        assert_eq!(merged_batch.density().dtype(), &BsxColumns::Density.dtype());

        Ok(())
    }

    #[test]
    fn test_merge_replicates_empty() {
        let batches: Vec<BsxBatch> = vec![];
        let result = merge_replicates(batches, Box::new(sum_agg), Box::new(mean_agg));
        assert!(result.is_err());
        assert!(result
            .unwrap_err()
            .to_string()
            .contains("batches cannot be empty"));
    }

    #[test]
    fn test_merge_replicates_single() -> anyhow::Result<()> {
        let batch1 = create_decoded_batch_1();
        let batches = vec![batch1.clone()]; // Clone batch1 for comparison
        let result_batch: BsxBatch =
            merge_replicates(batches, Box::new(sum_agg), Box::new(mean_agg))?;
        assert_eq!(result_batch, batch1); // Should return the original batch
        Ok(())
    }
}
