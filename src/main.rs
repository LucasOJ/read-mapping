mod fm_index;
mod nucleotide_stratified;
mod run_length_encoding;

use core::iter::Iterator;
use std::fs::File;
use std::io::{self, BufRead};

use fm_index::FMIndex;

fn main() {
    let file_path = "data/ncbi_dataset/data/GCF_000011505.1/GCF_000011505.1_ASM1150v1_genomic.fna";
    
    let file = File::open(file_path).expect("open file failed"); // Open the file
    let reader = io::BufReader::new(file);

    let lines: Vec<String> = reader
        .lines() // Get an iterator over lines
        .skip(1) // Skip the first line
        .filter_map(Result::ok) // Ignore any errors
        .collect();
    
    let real_text = lines.join("") + "$";
    // let text = "ACGCGCTTCGCCTT$";

    let fm_index = FMIndex::new(&real_text, 1, 128);

    // let target = "CGA";
    let target = "GCAAG";

    let start_indicies = fm_index.lookup(&target);

    // Debug
    println!("NUMBER OF MATCHES {}",start_indicies.len());
    for (index, &sa_index) in start_indicies.iter().enumerate() {
        
        // DEBUG
        let target_match = &real_text[sa_index..sa_index+target.len()];

        if target_match != target {
            println!("INDEX: {}, SA ENTRY: {}, ENTRY: {}",index,sa_index,target_match);
        }
    }

    let mut start = 0;
    let mut count = 0;
    while let Some(index) = real_text[start..].find(target) {
        let position = start + index;
        count += 1;
        start = position + 1;  // Move one character forward to allow overlapping matches
    }

    println!("RUST COUNT {}", count)

}
