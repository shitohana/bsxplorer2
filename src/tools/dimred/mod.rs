use rand::{thread_rng, Rng};
use rand::distributions::{Distribution, Uniform};
use rand_xoshiro::Xoshiro256Plus;
use rand_xoshiro::rand_core::{RngCore, SeedableRng};
use rayon::prelude::*;
use std::cmp;
use rand_distr::Normal;

/// Represents a random convolutional kernel for the ROCKET algorithm
#[derive(Clone, Debug)]
pub struct RandomKernel {
    /// Kernel weights
    pub weights: Vec<f32>,
    /// Dilation factor (spacing between kernel elements)
    pub dilation: usize,
    /// Whether to use padding
    pub padding: bool,
}

impl RandomKernel {
    /// Calculate the effective length of the kernel considering dilation
    pub fn effective_length(&self) -> usize {
        if self.weights.len() <= 1 {
            return self.weights.len();
        }
        (self.weights.len() - 1) * self.dilation + 1
    }

    /// Calculate the padding size if padding is enabled
    pub fn padding_size(&self) -> usize {
        if self.padding {
            (self.effective_length() - 1) / 2
        } else {
            0
        }
    }
}

/// Generates a set of random kernels with varying parameters
///
/// # Arguments
/// * `count` - Number of kernels to generate
/// * `max_input_length` - Maximum length of the input data_structs (used for dilation scaling)
/// * `seed` - Optional random seed for reproducibility
///
/// # Returns
/// A vector of RandomKernel structs
pub fn generate_random_kernels(
    count: usize,
    max_input_length: usize,
    seed: Option<u64>
) -> Vec<RandomKernel> {
    // Initialize random number generator with seed if provided
    // TODO change to seedable rng
    let mut rng = rand::thread_rng();
    // Distributions for sampling kernel parameters
    let length_options = [7, 9, 11];
    let length_dist = Uniform::from(0..length_options.len());
    let mut normal_dist = Normal::new(0.0, 1.0).expect("Failed to generate normal distribution");
    let bias_dist = Uniform::from(-1.0..1.0);

    (0..count)
        .map(|_| {
            // Randomly select kernel length
            let length_idx = length_dist.sample(&mut rng);
            let length = length_options[length_idx];

            // Generate weights from normal distribution
            let mut weights: Vec<f32> = (0..length)
                .map(|_| normal_dist.sample(&mut rng) as f32)
                .collect();

            // Mean center the weights
            let mean: f32 = weights.iter().sum::<f32>() / weights.len() as f32;
            for w in weights.iter_mut() {
                *w -= mean;
            }

            // Add bias term (as the last element in weights)
            let bias = bias_dist.sample(&mut rng) as f32;
            weights.push(bias);

            // Generate dilation - exponential scale
            // A = log_2((input_length - 1) / (kernel_length - 1))
            let a = if length > 1 && max_input_length > 1 {
                (max_input_length as f64 - 1.0).log2() - (length as f64 - 1.0).log2()
            } else {
                0.0
            };
            let x = rng.gen_range(0.0..=a);
            let dilation = (2.0f64.powf(x)).floor() as usize;
            dilation.max(1); // Ensure minimum dilation of 1

            // Randomly decide whether to use padding
            let padding = rng.gen_bool(0.5);

            RandomKernel {
                weights,
                dilation,
                padding,
            }
        })
        .collect()
}

/// Applies a kernel to input values using convolution
///
/// # Arguments
/// * `values` - Input methylation values (0.0 to 1.0)
/// * `kernel` - The RandomKernel to apply
///
/// # Returns
/// Vector of convolution results
pub fn apply_kernel(values: &[f32], kernel: &RandomKernel) -> Vec<f32> {
    if values.is_empty() || kernel.weights.is_empty() {
        return Vec::new();
    }

    let kernel_len = kernel.weights.len() - 1; // Last element is bias
    let dilation = kernel.dilation;
    let bias = *kernel.weights.last().unwrap();
    let padding_size = kernel.padding_size();

    let output_len = if kernel.padding {
        values.len()
    } else {
        values.len() - (kernel.effective_length() - 1)
    };

    if output_len <= 0 {
        return Vec::new();
    }

    let mut result = Vec::with_capacity(output_len);

    for i in 0..output_len {
        let mut sum = bias;

        for (j, &weight) in kernel.weights.iter().take(kernel_len).enumerate() {
            let input_idx = i + j * dilation;

            // Handle padding
            if kernel.padding {
                let adjusted_idx = input_idx as isize - padding_size as isize;
                if adjusted_idx >= 0 && adjusted_idx < values.len() as isize {
                    sum += values[adjusted_idx as usize] * weight;
                }
            } else if input_idx < values.len() {
                sum += values[input_idx] * weight;
            }
        }

        result.push(sum);
    }

    result
}

/// Extracts features from the convolution result (max value and PPV)
///
/// # Arguments
/// * `convolution_result` - Output from apply_kernel
///
/// # Returns
/// (max_value, proportion_of_positive_values)
pub fn extract_features_from_convolution(convolution_result: &[f32]) -> (f32, f32) {
    if convolution_result.is_empty() {
        return (0.0, 0.0);
    }

    // Find maximum value
    let max_value = convolution_result.iter()
        .fold(std::f32::NEG_INFINITY, |a, &b| a.max(b));

    // Calculate proportion of positive values (PPV)
    let positive_count = convolution_result.iter()
        .filter(|&&v| v > 0.0)
        .count();

    let ppv = positive_count as f32 / convolution_result.len() as f32;

    (max_value, ppv)
}

/// Extracts features from input values using multiple kernels
///
/// # Arguments
/// * `values` - Input methylation values (0.0 to 1.0)
/// * `kernels` - Vector of RandomKernel to apply
///
/// # Returns
/// Vector of features (alternating max value and PPV for each kernel)
pub fn extract_features(values: &[f32], kernels: &[RandomKernel]) -> Vec<f32> {
    let mut features = Vec::with_capacity(kernels.len() * 2);

    for kernel in kernels {
        let convolution = apply_kernel(values, kernel);
        let (max_val, ppv) = extract_features_from_convolution(&convolution);
        features.push(max_val);
        features.push(ppv);
    }

    features
}

/// Parallel version of feature extraction using Rayon
///
/// # Arguments
/// * `values` - Input methylation values (0.0 to 1.0)
/// * `kernels` - Vector of RandomKernel to apply
///
/// # Returns
/// Vector of features (alternating max value and PPV for each kernel)
pub fn extract_features_parallel(values: &[f32], kernels: &[RandomKernel]) -> Vec<f32> {
    // Process kernels in parallel, collect results
    let results: Vec<(f32, f32)> = kernels.par_iter()
        .map(|kernel| {
            let convolution = apply_kernel(values, kernel);
            extract_features_from_convolution(&convolution)
        })
        .collect();

    // Flatten results into a single feature vector
    let mut features = Vec::with_capacity(results.len() * 2);
    for &(max_val, ppv) in &results {
        features.push(max_val);
        features.push(ppv);
    }

    features
}

/// Processes multiple samples in parallel
///
/// # Arguments
/// * `samples` - Vector of samples, each a vector of methylation values
/// * `kernels` - Vector of RandomKernel to apply
///
/// # Returns
/// Vector of feature vectors, one per sample
pub fn process_samples_parallel(
    samples: &[Vec<f32>],
    kernels: &[RandomKernel]
) -> Vec<Vec<f32>> {
    samples.par_iter()
        .map(|sample| extract_features_parallel(sample, kernels))
        .collect()
}

/// Processes chunks of methylation data_structs in a streaming fashion
///
/// # Arguments
/// * `chunks_iterator` - Iterator yielding chunks of methylation data_structs
/// * `kernels` - Vector of RandomKernel to apply
///
/// # Returns
/// Combined feature vector
pub fn process_streaming<I, T>(
    chunks_iterator: I,
    kernels: &[RandomKernel]
) -> Vec<f32>
where
    I: Iterator<Item = T>,
    T: AsRef<[f32]>,
{
    // Initialize with zeros
    let feature_count = kernels.len() * 2;
    let mut max_features = vec![std::f32::NEG_INFINITY; kernels.len()];
    let mut positive_counts = vec![0usize; kernels.len()];
    let mut total_counts = vec![0usize; kernels.len()];

    // Process each chunk
    for chunk in chunks_iterator {
        let chunk_data = chunk.as_ref();

        for (i, kernel) in kernels.iter().enumerate() {
            let convolution = apply_kernel(chunk_data, kernel);

            // Update max
            if let Some(&chunk_max) = convolution.iter().max_by(|a, b| a.partial_cmp(b).unwrap_or(cmp::Ordering::Equal)) {
                max_features[i] = max_features[i].max(chunk_max);
            }

            // Update PPV counts
            let positive = convolution.iter().filter(|&&v| v > 0.0).count();
            positive_counts[i] += positive;
            total_counts[i] += convolution.len();
        }
    }

    // Combine results
    let mut features = Vec::with_capacity(feature_count);
    for i in 0..kernels.len() {
        features.push(max_features[i]);
        features.push(if total_counts[i] > 0 {
            positive_counts[i] as f32 / total_counts[i] as f32
        } else {
            0.0
        });
    }

    features
}

/// Converts methylation counts to binary values using binomial test
///
/// # Arguments
/// * `methylated_counts` - Vector of methylated read counts
/// * `total_counts` - Vector of total read counts
/// * `expected_rate` - Expected methylation rate (e.g., 0.5)
/// * `alpha` - Significance level (e.g., 0.05)
///
/// # Returns
/// Vector of methylation values (0.0 or 1.0 for binary, or ratio for non-significant)
pub fn apply_binomial_test(
    methylated_counts: &[u32],
    total_counts: &[u32],
    expected_rate: f64,
    alpha: f64
) -> Vec<f32> {
    assert_eq!(methylated_counts.len(), total_counts.len(),
               "Methylated and total count vectors must have same length");

    use statrs::distribution::{Binomial, DiscreteCDF};

    methylated_counts.iter()
        .zip(total_counts.iter())
        .map(|(&meth, &total)| {
            if total == 0 {
                return 0.0;
            }

            let ratio = meth as f32 / total as f32;

            // For small counts, just return the ratio
            if total < 10 {
                return ratio;
            }

            // Apply binomial test
            let dist = Binomial::new(expected_rate, total as u64).unwrap();

            // Two-tailed test
            let lower_tail = dist.cdf(meth as u64);
            let upper_tail = 1.0 - dist.cdf((meth as u64).saturating_sub(1));
            let p_value = lower_tail.min(upper_tail) * 2.0;

            if p_value <= alpha {
                // Significantly different from expected - return binary value
                if (meth as f64 / total as f64) > expected_rate {
                    1.0
                } else {
                    0.0
                }
            } else {
                // Not significantly different - return original ratio
                ratio
            }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn full_test() {
        // Example methylation data (could be streaming in real application)
        let methylated_counts: Vec<u32> = vec![5, 0, 10, 2, 8, 20, 0, 15];
        let total_counts: Vec<u32> = vec![10, 10, 10, 10, 10, 20, 10, 15];

        // Convert to values (binary representation optional)
        let use_binary = true;
        let methylation_values = if use_binary {
            apply_binomial_test(&methylated_counts, &total_counts, 0.5, 0.05)
        } else {
            methylated_counts.iter()
                .zip(total_counts.iter())
                .map(|(&m, &t)| if t > 0 { m as f32 / t as f32 } else { 0.0 })
                .collect()
        };

        // Generate random kernels
        let num_kernels = 1000;
        let max_length = 100; // Adjust based on your data window size
        let kernels = generate_random_kernels(num_kernels, max_length, Some(42));

        // Extract features using parallel processing
        let features = extract_features_parallel(&methylation_values, &kernels);

        println!("Extracted {} features", features.len());

        // In a real application, you'd process multiple samples and compare them
        // For example:
        let sample1 = vec![0.1, 0.5, 0.8, 0.2, 0.9, 0.1, 0.3, 0.7];
        let sample2 = vec![0.2, 0.6, 0.7, 0.3, 0.8, 0.2, 0.4, 0.6];

        let samples = vec![sample1, sample2];
        let all_features = process_samples_parallel(&samples, &kernels);

        println!("Processed {} samples", all_features.len());
    }
}