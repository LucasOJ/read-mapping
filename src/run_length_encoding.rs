#[derive(Debug)]
struct RunLengthEncodingEntry {
    char: char,
    count: usize
}


pub struct RunLengthEncodedString {
    seq: Vec<RunLengthEncodingEntry>,
    checkpoint_index: Vec<usize>

}

impl RunLengthEncodedString {
    pub fn new(str: &str, block_size: usize) -> Self {
        let mut maybe_prev_char: Option<char> = None;
        let mut count = 0;
        let mut rle_sequence = Vec::new();

        let mut checkpoint_index: Vec<usize> = Vec::new();

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
                checkpoint_index.push(rle_sequence.len());
            }

            count += 1;
            maybe_prev_char = Some(char);
        }

        // TODO: think about this - one character we haven't pushed yet?
        // Could be the start of a k-block?
        if let Some(prev_char) = maybe_prev_char {
            rle_sequence.push(RunLengthEncodingEntry { char: prev_char, count });
        }

        // TODO: return RLE checkpoint index
        return RunLengthEncodedString {
            seq: rle_sequence,
            checkpoint_index
        };
    }
}