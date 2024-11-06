use bsx_rs::genome;
use bsx_rs::io::bsx::write::ConvertReportOptions;
use bsx_rs::io::report::types::ReportType;

fn main() {
    let res = ConvertReportOptions::default()
        .finish(
            ReportType::BEDGRAPH,
            "/Users/shitohana/Documents/CX_reports/A_thaliana.bedGraph",
            "/Users/shitohana/Desktop/RustProjects/bsxplorer2_dev/bedgraph_test.ipc",
            "/Users/shitohana/Documents/CX_reports/old/arabidopsis.fa",
            "/Users/shitohana/Documents/CX_reports/old/arabidopsis.fa.fai",
        );
    res.unwrap();
}