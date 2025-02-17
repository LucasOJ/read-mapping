use std::cmp::max;
use std::fmt::Debug;

use crate::nucleotide_stratified::NucStratified;

#[derive(Debug)]
pub struct FMIndex {
    bwt: String,
    partial_suffix_array: Vec<usize>,
    first_bwt_column_index: NucStratified<usize>,
    ranks: NucStratified<Vec<usize>>,
    rank_lookup_thinning_factor: usize
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

        FMIndex { 
            bwt,
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
            rank_lookup_thinning_factor
        }
    }

    pub fn lookup(&self, target_str: &str) {
        // Inclusive
        let mut bottom = 0;
        
        // Exclusive
        let mut top = self.bwt.len();

        // Given a top and bottom, each iteration should extend by one charecter and find new top and bottom in first column of bwt table
        // Once no charecters left, lookup range in SA
        for nucleotide in target_str.chars().rev() {
            if (top <= bottom) {
                println!("NO MATCH");
                return
            }

            let first_bwt_rank = self.get_rank_for_index(nucleotide, bottom);
            let last_bwt_rank = self.get_rank_for_index(nucleotide, top - 1);

            let first_column_bottom_index = self.first_bwt_column_index.get(nucleotide);

            // Last-to-First mapping - chars of the same rank in the first and last column are the same

            /*
             * Why the max?
             * 
             * Note max(first_bwt_rank,1) - 1 is the same as max(first_bwt_rank - 1,0) wihtout underflow.
             * 
             * If `first_bwt_rank` >= 1, we have seen at least one occurence of `nucleotide` by the 
             * `bottom`th nucleotide in the bwt. In which case we -1 to convert rank to index offset.
             * 
             * If `first_bwt_rank` = 0, no seen any occurences yet, so we add no index offset.
             */ 
            bottom = first_column_bottom_index + max(first_bwt_rank,1) - 1;
            
            // rank-to-offset -1 here cancelled by +1 since top is exclusive of the interval
            top = first_column_bottom_index + last_bwt_rank;
        }

        println!("MATCH INTERVAL [{},{})",bottom,top);

        println!("SUBSTRING START INDICIES {:?}", &self.partial_suffix_array[bottom..top]);
        
    }

    // Should panic if gets wrong index
    fn get_rank_for_index(&self, nucleotide: char, index: usize) -> usize {
        let rank_checkpoint_index = index / self.rank_lookup_thinning_factor;

        let rank_checkpoint = self.ranks.get(nucleotide)
                                        .get(rank_checkpoint_index)
                                        .expect("coarse_rank_start_index not in range");

        let bwt = &self.bwt;

        // todo might be inefficient - probably better to for loop and count
        // let missing_nuc_instances = &bwt[rank_checkpoint_index..index].matches(nucleotide).count();
        // TODO: add in later when using rank lookup thinning
        let missing_nuc_instances = 0;

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