mod fm_index;
mod nucleotide_stratified;
mod run_length_encoding;

use core::iter::Iterator;
use core::option::Option;
use std::collections::HashMap;
use std::fmt::Debug;
use std::fs::File;
use std::io::{self, BufRead};
use std::iter;
use std::cmp::{max,min};

use fm_index::FMIndex;
use run_length_encoding::RunLengthEncodedString;



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

    let rle = RunLengthEncodedString::new(&real_text, 128);

    let text = "ACGCGCTTCGCCTT$";

    let fm_index = FMIndex::new(text, 1, 1);

    println!("{:?}",fm_index);

    fm_index.lookup("CGC");
}
