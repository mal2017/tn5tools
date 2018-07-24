use rust_htslib::bam;
use rust_htslib::prelude::*;

pub fn trim_cigar_string_tn5(cigar: bam::record::CigarString, is_rev: bool) {


	let c_owned: bam::record::CigarString = if is_rev {
		bam::record::CigarString(cigar.iter().rev().map(|a|a.clone()).collect())
	} else {
		cigar
	};


	let c_new: Vec<bam::record::Cigar> = Vec::new();

	for c in c_owned.iter() {
		trim_single_cigar(c, 5);
	}

}

fn trim_single_cigar(c: &bam::record::Cigar, left_to_trim: u32) {
	use rust_htslib::bam::record::Cigar;
	match c {
		Cigar::Match(u32) | Cigar::Equal(u32) | 
			Cigar::Diff(u32) | Cigar::Ins(u32) |
			Cigar::SoftClip(u32) => println!("consume"),
		_ => return(None),
	}

}

//pub consume_cigar()