extern crate bio;
extern crate indicatif;
extern crate rust_htslib;



// for pulling regions from bam
use rust_htslib::bam::IndexedReader;
pub mod regions;


// for shifting reads in iterator
//pub mod shift;

// pub mod count

// pub mod sizes