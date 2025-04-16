use bio::io::fasta::Record;
use rand_chacha::rand_core::TryRngCore;

fn generate_sequence<R: TryRngCore>(
    rng: &mut R,
    length: usize,
) -> String {
    let mut seq = String::new();
    let chars = ['A', 'C', 'G', 'T'];
    for _ in 0..length {
        seq.push(chars[(rng.try_next_u32().unwrap() % 4) as usize]);
    }
    seq
}

fn generate_genome<R: TryRngCore>(
    rng: &mut R,
    chr_number: usize,
    length: usize,
) -> Vec<Record> {
    let mut records = Vec::new();
    for i in 0..chr_number {
        let seq = generate_sequence(rng, length);
        let name = format!("chr{}", i);
        records.push(Record::with_attrs(&name, None, &seq.as_bytes()));
    }
    records
}

struct DemoReportBuilder {
    chr_lengths: Vec<usize>,
    seed: Option<u64>,
    max_coverage: Option<usize>,
    min_coverage: Option<usize>,
}
