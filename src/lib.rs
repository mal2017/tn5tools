extern crate bio;
extern crate indicatif;
extern crate rust_htslib;
extern crate csv;
extern crate rayon;
extern crate serde;
#[macro_use(array)]
extern crate ndarray;

// for pulling regions from bed
use rust_htslib::bam::IndexedReader;
pub mod regions;

pub mod count;

// pub mod sizes