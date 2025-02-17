use std::fmt::Debug;

// TODO: Should handle all nucleotide switching but still some manual in FMIndex
#[derive(Debug)]
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
             _  => panic!("NOT A NUCLEOTIDE!") 
        }
    }
}