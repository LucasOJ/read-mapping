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

    // println!("{:?}",fm_index);

    // MATCHES RUST COUNT FOR THIS TARGET
    // let target = "CGATTGTT";

    // TODO: GETS A WRONG ANSWER WITH THIS TARGET
    // INCLUDES 1 NON-MATCH -> BOTTOM SA INTERVAL INDEX IS TOO LOW 
    // let target = "CCGCTTGTTGA";

    // INCLUDES 2 NON-MATCHES when thinning factor = 128: 1 too high, 1 too low
    // Includes 1 NON-MATCH when thinning factor = 1: 1 too low
    let target = "CGATTTT";

    let start_indicies = fm_index.lookup(&target);
    println!("NUMBER OF MATCHES {}",start_indicies.len());
    for i in start_indicies {
        println!("INDEX: {}, ENTRY: {}",i,&real_text[i..i+target.len()]);
    }

    println!("RUST MATCH COUNT {}", real_text.match_indices(target).count());

}
