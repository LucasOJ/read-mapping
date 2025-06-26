# read-mapping

## Overview

This project is a read aligner inspired by [Bowtie](https://bowtie-bio.sourceforge.net/index.shtml). 

Given a genome, `ReadMappingIndex` constructs a compressed index which can map reads from experiments (eg RNA-seq) to their location within the genome. 

## ReadMappingIndex

`ReadMappingIndex` consists of 2 [FM-Indexes](https://en.wikipedia.org/wiki/FM-index), which use the [Burrows-Wheeler transform (BWT)](https://en.wikipedia.org/wiki/Burrows%E2%80%93Wheeler_transform) and [Run-length encoding](https://en.wikipedia.org/wiki/Run-length_encoding) to compress the genome. 

On my Intel i5-6200U laptop, I have benchmarked ~36 million reads / hour matching 47% of reads ([SRR11998244](https://www.ncbi.nlm.nih.gov/sra/SRR11998244)) for the [MRSA252](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id=282458) genome.

## Usage

```rust
// Construct ReadMappingIndex
let read_mapping_index = ReadMappingIndex::new(&genome);

// Map a read to the genome using up to 3 seeds of size 25 
let read_map_result = read_mapping_index.map_read(&read, 25, 3); 

// Save index to file
read_mapping_index.to_file("genome_index.bin");

// Load index from file
let read_mapping_index = ReadMappingIndex::from_file("genome_index.bin");
```

A FASTA parser is included to load genomes from file. A FASTQ parser is also included to stream reads from a file (see `main.rs` for an example).

## Mapping Strategy

`ReadMappingIndex` tolerates mismatches by using a seed-and-extend strategy. First, it attempts to find an exact match between a short section of the read (the seed) and the genome. The seed acts like a fingerprint for the rest of the read, so needs a reasonable length (we choose 25 nucleotides) to reduce the likelihood of matching irrelevant parts of the genome by chance. Then, the index extends the seed as far as it can in both directions along the genome until it reaches a mismatch or the end of the read.

We assume that the start of the read is more accurate, so begin by using seeds at the start of the read and take seeds from later in the read only when necessary. On the MRSA252 example we found that for the first 1 million reads, 62% of matches came from the first-choice seed at the start of the read, 36% from the second-choice seed and 2% from the third seed. Adding more seeds increases the runtime since each read which fails to match attempts an FMIndex lookup for each seed.

To extend the seeds, we use a property of the FMIndex: given a specific nucleotide, we can find the nucleotide preceding it in the genome, and so on. This property allows us to extend the read backwards from the seed. To extend the seed in both directions `ReadMappingIndex` uses two FMIndexes: one for the original genome and another for the same genome reversed, since extending backwards in the reversed genome corresponds to extending forwards in the original genome.

## Profiling

`cargo flamegraph` was used for profiling, which revealed the majority of CPU time is spent traversing the run-length encoded BWT.

## TODO
- Increase matching rate by allowing mismatches during extension
- Use quality data from FASTQ to improve seed selection
- Perform banded local alignments instead of naive-extension 
- Allow gapped alignments over multiple exons
- Parallel read-mapping using Rayon
