mod fasta_parsing;
mod fastq_parsing;
mod fm_index;
mod nucleotide_stratified;
mod read_mapping_index;
mod run_length_encoding;

use core::iter::Iterator;
use std::env;

use fasta_parsing::read_fasta;
use fastq_parsing::read_fastq;

use read_mapping_index::{MapReadResult, ReadMappingIndex};

fn main() {
    env::set_var("RUST_BACKTRACE", "1");

    let genome =
        read_fasta("data/ncbi_dataset/data/GCF_000011505.1/GCF_000011505.1_ASM1150v1_genomic.fna");

    let read_mapping_index = ReadMappingIndex::new(&genome);

    // read_mapping_index.to_file("MRSAIndex.bin");
    // let read_mapping_index = ReadMappingIndex::from_file("MRSAIndex.bin");

    println!("Mapping Reads");

    let num_of_reads = 1000000;

    let reads = read_fastq("data/SRR11998244.fastq");

    let read_results = reads
        .take(num_of_reads)
        .map(|read| read_mapping_index.map_read(&read, 25, 3));

    let mut seed_attempt_map = [0, 0, 0];
    for read_result in read_results {
        if let Some(MapReadResult { seed_attempt, .. }) = read_result {
            seed_attempt_map[seed_attempt] += 1;
        }
    }

    let successful_maps_count: usize = seed_attempt_map.iter().sum();

    for (seed_attempt, &attempt_count) in seed_attempt_map.iter().enumerate() {
        let percentage = f64::round(100.0 * attempt_count as f64 / successful_maps_count as f64);
        println!("Seed attempt {seed_attempt} - {attempt_count} ({percentage}%)");
    }

    println!("{successful_maps_count} reads successfully mapped out of {num_of_reads} reads");
}
