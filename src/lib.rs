extern crate bio;
extern crate indicatif;
extern crate rust_htslib;
extern crate csv;
extern crate rayon;
extern crate serde;
#[macro_use(array)]
#[macro_use(s)]
extern crate ndarray;
extern crate linear_map;
extern crate itertools;
extern crate gnuplot;
#[macro_use] extern crate quick_error;

pub mod regions;

pub mod count;

pub mod profile;

pub mod tn5shift;

pub mod atacbam;

pub mod fragsizes;