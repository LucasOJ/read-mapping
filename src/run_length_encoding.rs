use std::cmp::min;

#[derive(Debug)]
struct RunLengthEncodingEntry {
    char: char,
    count: usize
}

#[derive(Debug)]
struct RunLengthEncodingCheckpoint {
    entry_index: usize,
    offset: usize
}

#[derive(Debug)]
pub struct RunLengthEncodedString {
    seq: Vec<RunLengthEncodingEntry>,
    index_to_seq_checkpoints: Vec<RunLengthEncodingCheckpoint>, 
    block_size: usize
}

impl RunLengthEncodedString {
    pub fn new(str: &str, block_size: usize) -> Self {
        let mut maybe_prev_char: Option<char> = None;
        let mut count = 0;
        let mut rle_sequence = Vec::new();

        let mut index_to_seq_checkpoints: Vec<RunLengthEncodingCheckpoint> = Vec::new();

        for (index, char) in str.chars().enumerate() {

            if let Some(prev_char) = maybe_prev_char { // true on every iteration except the first

                // Terminate the previous run if the character changes
                if prev_char != char {
                    rle_sequence.push(RunLengthEncodingEntry { char: prev_char, count });
                    count = 0;
                }
            }

            if index % block_size == 0 {
                let checkpoint = RunLengthEncodingCheckpoint {
                    entry_index: rle_sequence.len(),
                    offset: count
                };
                index_to_seq_checkpoints.push(checkpoint);
            }

            count += 1;
            maybe_prev_char = Some(char);
        }

        // TODO: think about this - one character we haven't pushed yet?
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

        if checkpoint_index == target_index {
            return 0;
        }

        let mut match_count = 0;

        let checkpoint = &self.index_to_seq_checkpoints[checkpoint_index / self.block_size];    

        // First block needs special care since the checkpoint index may be somewhere in the 
        // run - not neccessarily at the front
        let first_rle_entry = &self.seq[checkpoint.entry_index];

        if first_rle_entry.char == target_char {
            match_count += min(
                /*
                 * Number of chars left in the run after the checkpoint is
                 * first_rle_entry.count - checkpoint.offset
                 *   
                 * -1 since we don't want to count the first char
                 */
                first_rle_entry.count - checkpoint.offset - 1, 
                target_index - checkpoint_index + 1
            );
        }

        // index in the string at start of the next run
        let mut current_str_index = checkpoint_index + first_rle_entry.count - checkpoint.offset;
        let mut rle_seq_index = checkpoint.entry_index + 1;

        while current_str_index <= target_index {
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