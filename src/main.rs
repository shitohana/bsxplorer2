use bsx_rs::genome;

fn main() {
    let annotation = genome::AnnotationBuilder::from_gff(
        "/Users/shitohana/Documents/CX_reports/old/A_thaliana_genomic.gff"
    ).finish().unwrap();
    println!("{}", annotation);
}