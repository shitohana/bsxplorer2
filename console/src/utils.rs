use indicatif::{ProgressBar, ProgressStyle};

pub fn init_pbar(total: usize) -> _lib::exports::anyhow::Result<ProgressBar> {
    let progress_bar = ProgressBar::new(total as u64);
    progress_bar.set_style(
        ProgressStyle::default_bar()
            .template("{spinner:.green} [{elapsed_precise}, ETA: {eta}] [{bar:40.cyan/blue}] {pos:>5.green}/{len:5} {msg}")?
            .progress_chars("#>-"),
    );
    progress_bar.set_message("Processing...");
    Ok(progress_bar)
}
