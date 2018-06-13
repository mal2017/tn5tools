extern crate bio;
extern crate indicatif;
extern crate rust_htslib;
extern crate csv;


// for pulling regions from bed
use rust_htslib::bam::IndexedReader;
pub mod regions;

pub mod count;

// pub mod sizes