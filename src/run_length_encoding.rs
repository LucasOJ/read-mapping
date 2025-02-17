use std::cmp::min;

#[derive(Debug)]
struct RunLengthEncodingEntry {
    char: char,
    count: usize
}

#[derive(Debug)]
pub struct RunLengthEncodedString {
    seq: Vec<RunLengthEncodingEntry>,
    index_to_seq_checkpoints: Vec<usize>, 
    block_size: usize
}

impl RunLengthEncodedString {
    pub fn new(str: &str, block_size: usize) -> Self {
        let mut maybe_prev_char: Option<char> = None;
        let mut count = 0;
        let mut rle_sequence = Vec::new();

        let mut index_to_seq_checkpoints: Vec<usize> = Vec::new();

        for (index, char) in str.chars().enumerate() {
            let is_block_start = index % block_size == 0;

            if let Some(prev_char) = maybe_prev_char { // true on every iteration except the first

                // Terminate the previous run if the charecter changes or have reached a checkpoint
                if prev_char != char || is_block_start {
                    rle_sequence.push(RunLengthEncodingEntry { char: prev_char, count });
                    count = 0;
                }
            }

            if is_block_start {
                index_to_seq_checkpoints.push(rle_sequence.len());
            }

            count += 1;
            maybe_prev_char = Some(char);
        }

        // TODO: think about this - one character we haven't pushed yet?
        // Could be the start of a k-block?
        if let Some(prev_char) = maybe_prev_char {
            rle_sequence.push(RunLengthEncodingEntry { char: prev_char, count });
        }

        return RunLengthEncodedString {
            seq: rle_sequence,
            index_to_seq_checkpoints,
            block_size
        };
    }

    pub fn count_matches_from_checkpoint(&self, target_char: char, checkpoint_index: usize, target_index: usize) -> usize {        
        if checkpoint_index % self.block_size != 0 {
            panic!("Index {} is not a checkpoint for block size {}", checkpoint_index, self.block_size);
        }

        let mut match_count = 0;

        // Current index in the original string
        let mut current_str_index = checkpoint_index;

        // Finds the RLE seq entry from the checkpoints array
        let mut rle_seq_index = self.index_to_seq_checkpoints[checkpoint_index / self.block_size];


        while current_str_index < target_index {
            let rle_entry = &self.seq[rle_seq_index];

            if rle_entry.char == target_char {
                /*
                 * Use min in case target index in middle of RLE entry eg when matching on 'A':
                 * 012345
                 * AAAAAA -> ('A',6)
                 *   ^ 
                 * when target_index = 3, we don't want the count from the entry (6)
                 * +1 needed to convert from indicies to counts
                 */
                match_count += min(rle_entry.count, target_index - current_str_index + 1);
            }

            current_str_index += rle_entry.count;

            // Consider next RLE entry on next iteration
            rle_seq_index += 1;
        }

        return match_count;
    }
}