macro_rules! define_strings {
    (
        $($name:ident = $value:literal);*$(;)?
    ) => {
        $(
            pub const $name: &str = $value;
        )*
    };
}


pub mod dmr {
    define_strings! {
        PADJ =
            "Adjusted P-value threshold for DMR identification using \
            2D-Kolmogorov-Smirnov test. Segments with a p-value smaller \
            than specified will be reported as DMRs.";
        MERGE_P =
            "Mann-Whitney U-test P-value threshold for merging adjacent \
            segments during recursive segmentation. Smaller p-values \
            result in more iterations and fewer falsely merged \
            segments.";
        TOLERANCE =
            "Tolerance for merging adjacent segments after the Total Variation \
            denoising step (Condat's algorithm).  Smaller values result in more \
            segments being merged. Should be very small to avoid \
            over-segmentation after denoising.";
        L_COEF =
            "Coefficient by which `initial_l` is divided in each iteration of the \
            segmentation algorithm. Smaller values perform more segmentation \
            iterations.";
        L_MIN =
            "Minimum value for the regularization parameter.  The regularization \
            parameter is decreased during segmentation until it is smaller than \
            this value.";
        INITIAL_L =
            "Initial regularization parameter for the Condat algorithm.  Larger \
            values result in stronger smoothing.";
        MAX_DIST =
            "Maximum distance between adjacent cytosines in a segment.  Cytosines \
            further apart than this distance will be in separate segments.";
        DIFF_THRESHOLD =
            "Set minimum difference threshold. DMRs with an absolute difference in \
            methylation proportion between the two groups smaller than this value \
            will be discarded.";
        MIN_CYTOSINES =
            "Set minimum number of cytosines threshold. DMRs with cytosine count \
            below this value will be discarded.";
        MIN_COVERAGE =
            "Set coverage threshold. Cytosines with coverage below this value in \
            any of the samples will be discarded.";
        N_MISSING =
            "Set missing values threshold. Cytosines with no data_structs in more \
            than specified number of samples will be discarded.";
        CONTEXT =
            "Select cytosine methylation context. Only cytosines in this context will \
            be used for DMR calling. CG/CHG/CHH.";
        FORCE =
            "Automatically confirm selected paths.";
        OUTPUT =
            "Prefix for the generated output files.";
        GROUP_A =
            "Paths to BSX files of the first sample group.";
        GROUP_B =
            "Paths to BSX files of the second sample group.";
    }
}

pub mod dimred {
    define_strings!{
        FILES =
            "Paths to BSX files";
        OUTPUT =
            "Path to output file";
        CONTEXT =
            "Methylation context";
        COVERAGE =
            "Minimal coverage";
        BETA =
            "Beta parameter value for PELT algorithm. If none - BIC (ln(length)) \
            is selected.";
        MIN_SIZE =
            "Minimal segment size (number of cytosines)";
        CHUNK =
            "Chunk size for PELT algorithm";
        INTERSECTION =
            "Batch intersection size for parallel processing";
        JOINT =
            "Changepoint joint distance for parallel processing";
    }
}
