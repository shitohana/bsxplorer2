use clap::{Parser, Subcommand}; // For parsing command-line arguments.
use console::style; // For colored text output.
use dialoguer::Confirm; // For interactive user confirmations.
use glob::glob;
use std::path::PathBuf;

#[derive(Parser, Debug)]
#[command(
    author = env!("CARGO_PKG_AUTHORS"),
    version = env!("CARGO_PKG_VERSION"),
    about = env!("CARGO_PKG_DESCRIPTION"),
    long_about = None
)]
struct Cli {
    /// Select a subcommand to execute.
    #[command(subcommand)]
    command: Commands,
}

/// Define the available subcommands.
#[derive(Subcommand, Debug)]
enum Commands {
    /// Install a package or file(s) matching a pattern.
    ///
    /// The `package` argument accepts a file pattern (e.g. "*.pkg") and uses glob
    /// matching to find files. Use `--force` to bypass confirmation.
    Install {
        /// The package name or file pattern to install.
        #[arg(value_parser)]
        package: String,

        /// Optional flag to force installation without confirmation.
        #[arg(short, long)]
        force: bool,
    },

    /// A dummy subcommand for demonstration purposes.
    Update {
        /// The target to update.
        #[arg(value_parser)]
        target: String,
    },
}

/// The main entry point.
///
/// This function first expands the command–line arguments using `wild` (helpful for Windows
/// wildcard expansion), then parses them with `clap` and dispatches to the correct subcommand.
fn main() {
    // Expand command line arguments for wildcard support.
    // This is particularly useful on Windows, where the shell may not expand glob patterns.
    let args: Vec<String> = wild::args().collect();

    // Parse arguments using the `Cli` struct.
    let cli = Cli::parse_from(args);

    // Display a colorful welcome message.
    println!(
        "{}",
        style("Welcome to DummyPkg Manager!")
            .bold()
            .underlined()
            .green()
    );

    // Dispatch the command based on the provided subcommand.
    match cli.command {
        Commands::Install { package, force } => {
            handle_install(package, force);
        }
        Commands::Update { target } => {
            handle_update(target);
        }
    }
}

/// Handles the `install` subcommand.
///
/// This function demonstrates pattern matching using the `glob` crate, confirmation with `dialoguer`,
/// and colored output using `console`.
fn handle_install(package: String, force: bool) {
    println!(
        "{}",
        style("Starting installation process...").bold().blue()
    );

    // Check if the package argument is a pattern (e.g., contains '*' or '?').
    let is_pattern = package.contains('*') || package.contains('?');

    // Collect matching files.
    let files: Vec<PathBuf> = if is_pattern {
        // Use `glob` to expand the provided pattern into matching file paths.
        // In a production scenario, ensure proper error handling if the pattern is invalid.
        match glob(&package) {
            Ok(paths) => paths.filter_map(Result::ok).collect(),
            Err(e) => {
                eprintln!("{} {}", style("Error parsing pattern:").red(), e);
                return;
            }
        }
    } else {
        // If no pattern is detected, assume the input is a direct file or package name.
        vec![PathBuf::from(package)]
    };

    // If no files were found, report an error.
    if files.is_empty() {
        eprintln!("{}", style("No files matched the provided pattern.").red());
        return;
    }

    // Display the list of files that will be installed.
    println!("{}", style("Files to install:").bold());
    for file in &files {
        println!(" - {}", file.display());
    }

    // If the user did not specify the force flag, confirm before proceeding.
    if !force {
        let prompt = format!(
            "Do you want to proceed with installing {} item(s)?",
            files.len()
        );
        let confirmed = Confirm::new()
            .with_prompt(prompt)
            .default(true)
            .interact()
            .unwrap_or(false);

        if !confirmed {
            println!("{}", style("Installation aborted by the user.").red());
            return;
        }
    }

    // Simulate the installation process.
    for file in files {
        println!("{} {}", style("Installing:").green(), file.display());
        // Insert actual installation logic here.
    }

    println!(
        "{}",
        style("Installation completed successfully!").bold().green()
    );
}

/// Handles the `update` subcommand.
///
/// This is a dummy implementation for demonstration purposes.
fn handle_update(target: String) {
    println!("{} {}", style("Updating target:").yellow().bold(), target);
    // Insert actual update logic here.
    println!("{}", style("Update completed!").green());
}
