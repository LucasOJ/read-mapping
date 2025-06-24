use std::fs::File;
use std::io::{self, BufRead};

pub fn read_fasta(file_path: &str) -> String {
    print!("Loading FASTA in {} to memory", file_path);

    let file = File::open(file_path).expect("FASTA failed to open");
    let reader = io::BufReader::new(file);

    let lines: Vec<String> = reader
        .lines() // Get an iterator over lines
        .filter_map(Result::ok) // Ignore any errors
        .map(|line| line.to_ascii_uppercase())
        .filter(|line| line.chars().next().unwrap() != '>')
        .collect();
    
    let genome = lines.join("") + "$";

    println!("[DONE]");
    println!("Genome: {}...", &genome[..30]);

    return genome;
}