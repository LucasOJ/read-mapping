use std::collections::HashMap;
use std::fmt::Debug;

use suffix_array::SuffixArray;

use crate::nucleotide_stratified::NucStratified;
use crate::run_length_encoding::RunLengthEncodedString;

type SuffixArrayIndex = usize;

#[derive(Debug)]
pub struct FMIndex {
    compressed_bwt: RunLengthEncodedString,
    sampled_suffix_array: HashMap<usize, usize>,
    first_bwt_column_index: NucStratified<usize>,
    ranks: NucStratified<Vec<usize>>,
    rank_sampling_step_size: usize,
    suffix_array_sampling_step_size: usize,
    seq_len: usize,
}

impl FMIndex {
    // Assumes sentinel terminated string `str`
    pub fn new(
        str: &str,
        suffix_array_sampling_step_size: usize,
        rank_sampling_step_size: usize,
    ) -> Self {
        let suffix_array = construct_suffix_array(str);

        let bwt = construct_bwt(str, &suffix_array);
        let compressed_bwt = RunLengthEncodedString::new(&bwt, rank_sampling_step_size);

        /*
         * Constructing the BWT requires the full suffix array
         * Since 'sample_suffix_array' takes ownership of the suffix array we can only take a sample
         * after the BWT has been contructed
         */
        let sampled_suffix_array =
            sample_suffix_array(suffix_array, suffix_array_sampling_step_size);

        let mut counts: NucStratified<usize> = NucStratified::default();
        let mut ranks: NucStratified<Vec<usize>> = NucStratified::default();

        for (i, nuc) in bwt.chars().enumerate() {
            if nuc != '$' {
                // += operator does not do automatic dereferencing so have to do it ourselves
                *counts.get_mut(nuc) += 1;
            }

            if i % rank_sampling_step_size == 0 {
                for nuc_key in ['A', 'C', 'G', 'T'] {
                    let current_nuc_count = *counts.get(nuc_key);
                    ranks.get_mut(nuc_key).push(current_nuc_count);
                }
            }
        }

        /*
         * Stores the index at which each run begins in the first column of the BWT matrix
         *
         * Note we're using the final nuc counts from the BWT rather than the original string,
         * but the overall counts will match since the BWT is a permutation of the original string
         */
        let first_bwt_column_index = get_first_bwt_column_index(&counts);

        FMIndex {
            compressed_bwt,
            sampled_suffix_array,
            first_bwt_column_index,
            ranks,
            rank_sampling_step_size,
            suffix_array_sampling_step_size,
            seq_len: str.len(),
        }
    }

    pub fn lookup(&self, target_str: &str) -> impl Iterator<Item = SuffixArrayIndex> {
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
            let bottom_rank = match bottom {
                0 => 1,
                _ => self.get_rank_for_index(nucleotide, bottom - 1) + 1,
            };
            let top_rank = self.get_rank_for_index(nucleotide, top - 1);

            // No occourences of `nucleotide` in [bottom, top) so no matches
            if bottom_rank > top_rank {
                return 0..0; // Empty range
            }

            bottom = self.last_to_first_mapping(nucleotide, bottom_rank);

            // +1 since top is exclusive of the interval
            top = self.last_to_first_mapping(nucleotide, top_rank) + 1;
        }

        return bottom..top;
    }

    // Should panic if gets wrong index
    fn get_rank_for_index(&self, nucleotide: char, index: usize) -> usize {
        let checkpoint_index = index / self.rank_sampling_step_size;

        let rank_checkpoint = self
            .ranks
            .get(nucleotide)
            .get(checkpoint_index)
            .expect("coarse_rank_start_index not in range");

        // TODO: Fix multipliciation in argument list
        let missing_nuc_instances = self.compressed_bwt.count_matches_from_checkpoint(
            nucleotide,
            checkpoint_index * self.rank_sampling_step_size,
            index,
        );

        return rank_checkpoint + missing_nuc_instances;
    }

    pub fn get_genome_position(&self, target_index: SuffixArrayIndex) -> usize {
        let mut current_index = target_index;

        /*
         * General idea is to `walk back` through the string by repeatedly applying the LF mapping until we
         * find a suffix with a sampled SA entry.
         *
         * `self.sampled_suffix_array` samples the SA every `suffix_array_sampling_step_size` entries
         * by entry index - not the index of the suffix in the original string! This ensures that it takes
         * at most `suffix_array_sampling_step_size` steps to find a sampled SA entry.
         */
        for steps in 0..self.suffix_array_sampling_step_size {
            /*
             * Case 1: Hit a sampled SA entry.
             *
             * Note:
             * 1. The SA maps each char in the BWT to its position in the original string (by construction)
             * 2. In each step we walk back one char in the original string following the LF mapping
             *
             * so the original index of the suffix [target_index:] is
             * sampled entry + the numbers of chars we walked back.
             */
            if let Some(suffix_array_entry) = self.sampled_suffix_array.get(&current_index) {
                return suffix_array_entry + steps;
            }

            /*
             * Case 2: Current entry not sampled in SA array.
             *
             * 1. Lookup the charecter preceding the current suffix in the original string.
             * This is just the char at index `current_index` in the BWT (by construction).
             *
             * 2. Find the index of the row in the BW matrix that starts with the preceding
             * charecter we just found using the LF mapping
             */
            let nucleotide = self.compressed_bwt.get_char_from_position(current_index);
            let rank = self.get_rank_for_index(nucleotide, current_index);

            current_index = self.last_to_first_mapping(nucleotide, rank);
        }
        panic!("DIDN'T FIND SA ENTRY WHEN SHOULD HAVE")
    }

    // TODO: Could also return genome_pos of extension start
    // Match index is the row of the bwt where the match begins
    pub fn count_extension_matches(&self, match_index: SuffixArrayIndex, extension: &str) -> usize {
        let mut current_index = match_index;

        for (match_count, read_nuc) in extension.chars().enumerate() {
            // Read the previous character in the string
            let nucleotide = self.compressed_bwt.get_char_from_position(current_index);

            if nucleotide != read_nuc {
                return match_count;
            }

            // Find the BWT row assoicated with the preceding character
            let rank = self.get_rank_for_index(nucleotide, current_index);
            current_index = self.last_to_first_mapping(nucleotide, rank);
        }

        // Matched whole extension
        return extension.len();
    }

    /*
     * By construction of the BWT, charecters of the same rank in the first and last column refer to the same
     * charecter.
     *
     * To prove, imagine taking all the rows of the BW matrix ending with char C, and moving the C to the
     * front.
     * Note that:
     * 1. The strings generated are rotations of the original string, so are the set of all C-prefixed rows
     * in the BW Matrix
     * 2. Appending the C does not change the relative order of these rows
     * 3. All the rows beginning with C are contiguous so can be defined by an interval of rows (due to
     * the lexicographical ordering)
     *
     * This means that given a char C with rank i (C_i) in the BWT, we can find the corresponding
     * row in the BW matrix beginning with C_i.
     *
     * This function takes a char and its rank in the BWT, and returns the index of the row in the BW matrix
     * which starts with this char
     */
    fn last_to_first_mapping(&self, nucleotide: char, rank: usize) -> usize {
        let first_occurence_index = self.first_bwt_column_index.get(nucleotide);

        // -1 to convert rank to index
        return first_occurence_index + rank - 1;
    }
}

fn construct_suffix_array(str: &str) -> Vec<usize> {
    /*
     * suffix_array crate uses https://github.com/y-256/libdivsufsort under the hood.
     * 
     * Stripping the final character as the crate adds its own sentinel
     */ 
    let suffix_array = SuffixArray::new(&str[..str.len() - 1].as_bytes());

    return suffix_array
        .into_parts()
        .1
        .into_iter()
        .map(|x| x as usize)
        .collect();
}

// Note: This function takes ownership of the suffix_array since the sampled structure should have ownership of all the elements
fn sample_suffix_array(
    suffix_array: Vec<usize>,
    suffix_array_sampling_step_size: usize,
) -> HashMap<usize, usize> {
    /*
     * To ensure O(1) lookup of the SA entries which have not been sampled, sample every
     * `suffix_array_sampling_step_size`th entry by value.
     *
     * In `get_genome_position` we traverse back one suffix at a time in their lexicographical ordering,
     * so using this hashmap we will find a sampled SA entry in at most `suffix_array_sampling_step_size`
     * steps.
     */
    return suffix_array
        .into_iter()
        .enumerate()
        .filter(|(_, entry)| entry % suffix_array_sampling_step_size == 0)
        .collect(); // Maps index in original SA to entry in SA
}

fn construct_bwt(str: &str, suffix_array: &Vec<usize>) -> String {
    let str_bytes = str.as_bytes();
    let str_len = str.len();
    /*
     * The lexicographical ordering of suffixes given by the suffix array yeilds the BWT matrix
     * To the the BWT itself, we just take the charecter in the final column for each row.
     * This is the character before the start of each suffix in the original string by the cyclic
     * nature of the permuations making up the rows in the BWT matrix.
     *
     * For example, "GCCAT$" yields:
     * SA BWT Matrix
     * 5  $|G C C A [T]
     * 3  A T $|G C [C]
     * 2  C A T $|G [C]
     * 1  C C A T $|[G]
     * 0  G C C A T [$]
     * 4  T $|G C C [A]
     *
     */
    let bwt_bytes = suffix_array
        .iter()
        // Add str_len to (lex_pos - 1) since without it we underflow when lex_pos = 0
        .map(|&lex_pos| str_bytes[(lex_pos + str_len - 1) % str_len])
        .collect();
    return String::from_utf8(bwt_bytes).expect("not valid UTF8");
}

fn get_first_bwt_column_index(nuc_counts: &NucStratified<usize>) -> NucStratified<usize> {
    return NucStratified {
        // Since $ is less than all other chars the run of A's starts from position 1
        a: 1,
        c: 1 + nuc_counts.get('A'),
        g: 1 + nuc_counts.get('A') + nuc_counts.get('C'),
        t: 1 + nuc_counts.get('A') + nuc_counts.get('C') + nuc_counts.get('G'),
    };
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn lookup() {
        let fm_index = FMIndex::new("ACGCGCTTCGCCTT$", 3, 4);
        assert!(fm_index.lookup("ATG").eq(0..0));
        assert!(fm_index.lookup("C").eq(2..8));
        assert!(fm_index.lookup("CGC").eq(3..6));
    }

    #[test]
    fn count_extension_matches() {
        let fm_index = FMIndex::new("ACGCGCTTCGCCTT$", 3, 4);
        // Partial match against chars preceding CCTT$
        assert_eq!(fm_index.count_extension_matches(2, "GCTAAA"), 3);

        assert_eq!(fm_index.count_extension_matches(7, ""), 0);

        // Full match against chars preceding GCTTCGCCTT$
        assert_eq!(fm_index.count_extension_matches(10, "CGCA"), 4);
    }

    #[test]
    fn get_genome_position() {
        /*
         * Suffix array for ACGCGCTTCGCCTT$:
         *  0 14: $
         *  1  0: ACGCGCTTCGCCTT$
         *  2 10: CCTT$
         *  3  8: CGCCTT$
         *  4  1: CGCGCTTCGCCTT$
         *  5  3: CGCTTCGCCTT$
         *  6 11: CTT$
         *  7  5: CTTCGCCTT$
         *  8  9: GCCTT$
         *  9  2: GCGCTTCGCCTT$
         * 10  4: GCTTCGCCTT$
         * 11 13: T$
         * 12  7: TCGCCTT$
         * 13 12: TT$
         * 14  6: TTCGCCTT$
         */
        let fm_index = FMIndex::new("ACGCGCTTCGCCTT$", 3, 4);

        let actual_sa_values: Vec<_> = (0..15).map(|i| fm_index.get_genome_position(i)).collect();

        let expected_sa_values: Vec<usize> =
            Vec::from([14, 0, 10, 8, 1, 3, 11, 5, 9, 2, 4, 13, 7, 12, 6]);

        assert_eq!(actual_sa_values, expected_sa_values);
    }
}
