use rust_htslib::bam;
use rust_htslib::prelude::*;
use std::iter;


// 1. expand cigar
// 2. consume cigar while following rules
// 3. as this happens, build read specific rule for sequence


fn expand_cigarstring(cs: &bam::record::CigarString) -> Vec<char> {
	cs.iter()
	  .flat_map(|a| vec![a.char(); a.len() as usize])
	  .collect::<Vec<char>>()
}


pub fn trim_cigar_string_tn5(cs: &bam::record::CigarString, is_rev: &bool) {

	let mut cs_vec = expand_cigarstring(cs);
	let rng = if *is_rev {
		(0..(cs_vec.len() - 1)).rev()
	} else {
		(0..(cs_vec.len() - 1))
	};

	println!("{:?}", is_rev);
	for i in rng {
		println!("{:?}", i);
	}
}