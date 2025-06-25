use std::fmt::Debug;

#[derive(Debug, Default)]
pub struct NucStratified<T: Debug> {
    pub a: T,
    pub c: T,
    pub g: T,
    pub t: T
}

impl<T: Debug> NucStratified<T> {
    pub fn get(&self, nucleotide: char) -> &T {
        match nucleotide {
            'A' => &self.a,
            'C' => &self.c,
            'G' => &self.g,
            'T' => &self.t,
             _  => panic!("{nucleotide} IS NOT A NUCLEOTIDE!") 
        }
    }

    pub fn get_mut(&mut self, nucleotide: char) -> &mut T {
        match nucleotide {
            'A' => &mut self.a,
            'C' => &mut self.c,
            'G' => &mut self.g,
            'T' => &mut self.t,
             _  => panic!("{nucleotide} IS NOT A NUCLEOTIDE!") 
        }
    }
}