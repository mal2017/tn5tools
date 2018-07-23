use rust_htslib::bam;
use rust_htslib::prelude::*;

pub fn trim_cigar_pos(c: bam::record::CigarString, readlen: u32) {
	//let new_cig = bam::record::Cigar::Match(0);

	let len_c = c.len();
	let mut consumed = 0;
	let mut idx = 0;

	// TODO init vec here
	let mut new_cig_vec = Vec::new();

	while consumed < 5 {
		let current_cigar_len  = c[idx].len();
		let current_cigar_char = c[idx].char();

		if  current_cigar_len >= readlen {
			new_cig_vec.push(format!("{}{}",current_cigar_len-5,current_cigar_char));
		}	
		consumed +=5;
	}

}