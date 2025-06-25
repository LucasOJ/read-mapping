use std::fs::File;
use std::io::{BufReader, BufRead};

pub fn read_fastq(file_path: &str) -> impl Iterator<Item = String> {
    let file = File::open(file_path).expect("FASTQ failed to open");
    let reader = BufReader::new(file);

    // TODO: Make more robust - assumes that the 1st, 5th, 9th, 13th lines are the reads 
    let reads = reader
        .lines() // Get an iterator over lines
        .filter_map(Result::ok) // Ignore any errors
        .skip(1)
        .step_by(4)
        .filter(|line| line.chars().all(|c| matches!(c, 'A' | 'T' | 'C' | 'G')));

    return reads;
}