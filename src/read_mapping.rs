use std::{cmp::min, collections::HashMap};

use crate::fm_index::FMIndex;

pub struct ReadMapper {
    // forwards_fm_index: FMIndex,
    reverse_fm_index: FMIndex,
    genome_length: usize,
}

impl ReadMapper {
    pub fn new(genome: &str) -> Self {
        // let forwards_fm_index = FMIndex::new(&genome, 128, 128);

        let mut reversed_genome: String = genome.chars().rev().collect();
        reversed_genome.push('$');

        let reverse_fm_index = FMIndex::new(&reversed_genome, 128, 128);

        return ReadMapper {
            // forwards_fm_index,
            reverse_fm_index,
            genome_length: genome.len() + 1, // +1 for sentinel
        };
    }

    pub fn to_file(self, filename: &str) {
        todo!()
    }

    pub fn from_file(self, filename: &str) {
        todo!()
    }

    // TODO: Return type should include seed number (and location?)
    pub fn map_read(
        &self,
        read: &str,
        seed_length: usize,
        max_seeds: usize,
    ) -> Option<(usize, usize)> {
        assert!(read.len() >= seed_length);
        let reversed_read: String = read.chars().rev().collect();

        for i in 0..min(read.len() / seed_length, max_seeds) {
            /*
             * Seed bounds in the forward direction
             * We start by choosing seeds from the beginning of the reads since
             * reads are typically more accurate at the beginning
             */
            let seed_start_index = i * seed_length;
            let seed_end_index = seed_start_index + seed_length;

            // Corresponding bounds for the same seed in the reverse direction
            let reverse_seed_start_index = read.len() - seed_end_index;
            let reverse_seed_end_index = read.len() - seed_start_index;
            let reverse_seed = &reversed_read[reverse_seed_start_index..reverse_seed_end_index];

            /*
             * We check for the seed in the reverse genome as often if there is a match for the 
             * first seed (which happens often), then we can extend to the rest of the read without
             * consulting the forwards_fm_index
             */
            let mut reverse_lookup_results = self.reverse_fm_index.lookup(&reverse_seed).peekable();

            // Try next seed if current seed not found in genome
            if reverse_lookup_results.peek().is_none() {
                continue;
            }

            // Get the read extension from where the seed ends
            let reverse_extension = &read[seed_end_index..];

            let mut extension_length_map: HashMap<_, _> = reverse_lookup_results
                .map(|lookup_result| {
                    let reverse_genome_pos =
                        self.reverse_fm_index.get_genome_position(lookup_result);
                    let extension_length = self
                        .reverse_fm_index
                        .count_extension_matches(lookup_result, reverse_extension);

                    /*
                     * Converts to a genome_pos in the forwards genome.
                     * 
                     * `self.genome_length - reverse_genome_pos` is classic index reversal
                     * 
                     * `-seed_length` since the beginning of the seed in the reverse genome 
                     * is the end of the seed in the forwards genome
                     * 
                     * `-1` since the sentinel is at the end of the genomes in both cases
                     */
                    let genome_pos = self.genome_length - reverse_genome_pos - seed_length - 1;
                    let match_length = seed_length + extension_length;
                    return (genome_pos, match_length);
                })
                .collect();

            // TODO: If length of forwards continuation != 0, try adding forwards

            return extension_length_map
                .into_iter()
                .max_by_key(|(_, length_ref)| *length_ref);
        }

        // None of the seeds could be found
        return None;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn map_read() {
        let read_mapper = ReadMapper::new("ATACTTTATCAAATGTAAAAGTATCTCCTTCGTTTACGTCTAATTTTT");
        assert_eq!(read_mapper.map_read("ATAC", 4, 1), Some((0, 4)));
        assert_eq!(
            read_mapper.map_read("ATACTTTATCAAATGTAA", 5, 2),
            Some((0, 18))
        );
        assert_eq!(read_mapper.map_read("ATCAAATGTAAAAG", 7, 2), Some((7, 14)));
        assert_eq!(read_mapper.map_read("ATCAATTGTAAAA", 7, 2), None);
        assert_eq!(read_mapper.map_read("ATCCAATGTAAAAG", 4, 1), None);
    }

    #[test]
    #[should_panic]
    fn map_read_panics_on_invalid_params() {
        let read_mapper = ReadMapper::new("ATACTTTATCAAATGTAAAAGTATCTCCTTCGTTTACGTCTAATTTTT");
        read_mapper.map_read("ATAC", 5, 4);
    }
}
