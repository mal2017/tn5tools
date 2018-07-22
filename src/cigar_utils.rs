use rust_htslib::bam;
use rust_htslib::prelude::*;

pub fn trim_cigar_pos(c: bam::record::CigarString, readlen: u32) {
	let new_cig = bam::record::Cigar::Match(0);

	let len_c = c.len();
	let mut consumed = 0;
	let mut idx = 0;


	while consumed < 5 {
		if c[idx].len() < readlen {
			println!("{:?} - {:?}", c[idx].len(), c[idx].char());
		}	
		consumed +=5;
	}

}