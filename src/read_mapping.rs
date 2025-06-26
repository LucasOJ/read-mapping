use bincode::decode_from_std_read;
use bincode::{Encode, Decode,config::standard, encode_into_std_write, error::EncodeError};

use std::{cmp::min, collections::HashMap};
use std::fs::File;

use crate::fm_index::FMIndex;

#[derive(Encode, Decode)]
pub struct ReadMapper {
    forwards_fm_index: FMIndex,
    reverse_fm_index: FMIndex,
    genome_length: usize,
}

impl ReadMapper {
    pub fn new(genome: &str) -> Self {

        // TODO: Not nice having to reallocate `str` just to add the sentinel
        let mut forwards_genome = String::from(genome);
        forwards_genome.push('$');

        let forwards_fm_index = FMIndex::new(&forwards_genome, 128, 128);

        let mut reversed_genome: String = genome.chars().rev().collect();
        reversed_genome.push('$');

        let reverse_fm_index = FMIndex::new(&reversed_genome, 128, 128);

        return ReadMapper {
            forwards_fm_index,
            reverse_fm_index,
            genome_length: genome.len() + 1, // +1 for sentinel
        };
    }

    pub fn to_file(&self, filename: &str) -> Result<usize, EncodeError> {
        let mut file = File::create(filename).unwrap();
        encode_into_std_write(self, &mut file, standard())
    }

    pub fn from_file(filename: &str) -> Self {
        let mut file = File::open(filename).unwrap();
        decode_from_std_read(&mut file, standard()).unwrap()
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

            let mut match_length_map: HashMap<_, _> = reverse_lookup_results
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

            let forwards_extension = &reversed_read[read.len() - seed_start_index..read.len()];

            /*
             * Only consider forwards_fm_index if we have something to gain.
             * Note that the forwards_extension_length can be at most the seed_length,
             * as if it was longer then the previous seed would have matched in the 
             * reverse_fm_index lookup.
             */ 
            if forwards_extension.len() > 0 {
                let forwards_seed = &read[seed_start_index..seed_end_index];

                let forwards_lookup_results = self.forwards_fm_index.lookup(forwards_seed);

                for lookup_result in forwards_lookup_results {
                    let forwards_genome_pos =
                        self.forwards_fm_index.get_genome_position(lookup_result);
                    let forwards_extension_length = self
                        .forwards_fm_index
                        .count_extension_matches(lookup_result, forwards_extension);

                    // The position where the left-extension in the forwards genome begins
                    let genome_pos = forwards_genome_pos - forwards_extension_length;

                    /*
                     * If we already have the length of a reverse match from the same seed extract it
                     * 
                     *       <--forwards_extension_length--|-----reverse_match_length------>
                     *       <-----forwards_extension------|--seed--|--reverse_extension--->
                     * ------|-----------------------------|--------|----------------------|------
                     *       ^                             ^
                     *    genome_pos              forwards_genome_pos
                     * 
                     * Should always be something here since the reverse lookup should populate it 
                     * with at least the length of the seed. Hence, we do `unwrap` to panic if 
                     * something has gone wrong
                     */ 
                    let reverse_match_length =
                        match_length_map.remove(&forwards_genome_pos).unwrap();

                    // TODO: Can there already be something at genome_pos?
                    match_length_map
                        .insert(genome_pos, forwards_extension_length + reverse_match_length);
                }
            }

            // Return longest match
            return match_length_map
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
        
        // |------| is used to highlight the matches between the reads and the genome

        //                               |--|
        assert_eq!(read_mapper.map_read("ATAC", 4, 1), Some((0, 4)));

        assert_eq!(
            //                    |----------------|
            read_mapper.map_read("ATACTTTATCAAATGTAA", 5, 2),
            Some((0, 18))
        );

        //                               |------------|
        assert_eq!(read_mapper.map_read("ATCAAATGTAAAAG", 7, 2), Some((7, 14)));

        assert_eq!(read_mapper.map_read("ATCAATTGTAAAA", 7, 2), None);

        assert_eq!(read_mapper.map_read("ATCCAATGTAAAAG", 4, 1), None);

        assert_eq!(
            //                     |---------------|
            read_mapper.map_read("TTACTTTATCAAATGTAA", 5, 2),
            Some((1, 17))
        );
        
        assert_eq!(
            //                      |--------------|
            read_mapper.map_read("AAACTTTATCAAATGTAA", 5, 2),
            Some((2, 16))
        );

        assert_eq!(
            //                              |--------------|
            read_mapper.map_read("GTATCTTCTACGTTTACGTCTAATTT", 7, 3),
            Some((30, 16))
        );

        assert_eq!(
            //                              |-----------|
            read_mapper.map_read("GTATCTTCTACGTTTACGTCTAAATT", 7, 3),
            Some((30, 13))
        );
    }

    #[test]
    #[should_panic]
    fn map_read_panics_on_invalid_params() {
        let read_mapper = ReadMapper::new("ATACTTTATCAAATGTAAAAGTATCTCCTTCGTTTACGTCTAATTTTT");
        read_mapper.map_read("ATAC", 5, 4);
    }
}
