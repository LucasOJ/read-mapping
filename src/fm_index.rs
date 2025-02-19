use std::fmt::Debug;

use crate::nucleotide_stratified::NucStratified;
use crate::run_length_encoding::RunLengthEncodedString;

#[derive(Debug)]
pub struct FMIndex {
    compressed_bwt: RunLengthEncodedString,
    partial_suffix_array: Vec<usize>,
    first_bwt_column_index: NucStratified<usize>,
    ranks: NucStratified<Vec<usize>>,
    rank_lookup_thinning_factor: usize,
    seq_len: usize
}

impl FMIndex {
    // Assumes sentinel terminated string
    pub fn new(str: &str, suffix_array_thinning_factor: usize, rank_lookup_thinning_factor: usize) -> Self {
        let suffix_array = construct_suffix_array(str);
        let bwt = construct_bwt(str, &suffix_array);
        let partial_suffix_array = suffix_array.into_iter()
            .step_by(suffix_array_thinning_factor)
            .collect();

        // TODO better way than this?
        let mut a_count = 0;
        let mut c_count = 0;
        let mut g_count = 0;
        let mut t_count = 0;

        let mut a_rank = Vec::new();
        let mut c_rank = Vec::new();
        let mut g_rank = Vec::new();
        let mut t_rank = Vec::new();

        for (i, nuc) in bwt.chars().enumerate() {
            match nuc {
                'A' => a_count += 1,
                'C' => c_count += 1,
                'G' => g_count += 1,
                'T' => t_count += 1,
                x => println!("Got a {}", x)
            }

            // TODO: think about sentinel char!
            if i % rank_lookup_thinning_factor == 0 {
                a_rank.push(a_count);
                c_rank.push(c_count);
                g_rank.push(g_count);
                t_rank.push(t_count);
            }
        }

        let compressed_bwt = RunLengthEncodedString::new(&bwt, rank_lookup_thinning_factor);

        FMIndex { 
            compressed_bwt,
            partial_suffix_array,

            // TODO: rename and find better way of doing this
            first_bwt_column_index: NucStratified {
                a: 1,
                c: 1 + a_count,
                g: 1 + a_count + c_count, 
                t: 1 + a_count + c_count + g_count,
            },

            ranks: NucStratified {
                a: a_rank,
                c: c_rank,
                g: g_rank,
                t: t_rank
            },
            rank_lookup_thinning_factor,
            seq_len: str.len()
        }
    }

    pub fn lookup(&self, target_str: &str) -> Vec<usize> {
        println!("STARTED LOOKUP FOR {}", target_str);
        // Inclusive
        let mut bottom = 0;
        
        // Exclusive
        let mut top = self.seq_len;

        // Given a top and bottom, each iteration should extend by one charecter and find new top and bottom in first column of bwt table
        // Once no charecters left, lookup range in SA
        for nucleotide in target_str.chars().rev() {
            /*
             * The rank of the first `nucleotide`` in the [bottom, top) range is one larger
             * than the largest rank seen before the range.
             * If bottom = 0, the first rank in any range would be 1.
             */
            let first_bwt_rank = match bottom {
                0 => 1,
                _ => self.get_rank_for_index(nucleotide, bottom - 1) + 1
            };
            
            let last_bwt_rank = self.get_rank_for_index(nucleotide, top - 1);

            // No occourences of `nucleotide` in [bottom, top) so no matches
            if first_bwt_rank == last_bwt_rank {
                return Vec::new();
            }

            let first_column_bottom_index = self.first_bwt_column_index.get(nucleotide);

            // Last-to-First mapping - chars of the same rank in the first and last column are the same

            // -1 for rank to index offset
            bottom = first_column_bottom_index + first_bwt_rank - 1;
            
            // rank-to-offset -1 here cancelled by +1 since top is exclusive of the interval
            top = first_column_bottom_index + last_bwt_rank;
        }

        // println!("MATCH INTERVAL [{},{})",bottom,top);

        let substr_start_indicies = self.partial_suffix_array[bottom..top].to_vec();
        return substr_start_indicies;

    }

    // Should panic if gets wrong index
    fn get_rank_for_index(&self, nucleotide: char, index: usize) -> usize {
        let checkpoint_index = index / self.rank_lookup_thinning_factor;

        let rank_checkpoint = self.ranks.get(nucleotide)
                                        .get(checkpoint_index)
                                        .expect("coarse_rank_start_index not in range");

        // TODO: Fix multipliciation in argument list
        let missing_nuc_instances = self.compressed_bwt.count_matches_from_checkpoint(nucleotide, checkpoint_index * self.rank_lookup_thinning_factor, index);

        return rank_checkpoint + missing_nuc_instances;
    }
}

fn construct_suffix_array(str: &str) -> Vec<usize> {
    let mut suffix_array: Vec<usize> = (0..str.len()).collect();

    // Sort suffix indices based on their corresponding suffix slices
    suffix_array.sort_by_key(|&i| &str[i..]);
    
    // println!("Suffix Array {:?}", suffix_array);

    return suffix_array;
}

fn construct_bwt(str: &str, suffix_array: &Vec<usize>) -> String{
    let str_bytes = str.as_bytes();
    let str_len = str.len();
    let bwt_bytes = suffix_array.iter()
                 .map(|&lex_pos| str_bytes[(lex_pos + str_len - 1) % str_len]) // Add str_len so we mod +ve numbers to [0,str_len) range
                 .collect();
    return String::from_utf8(bwt_bytes).expect("not valid UTF8");
}