use rust_htslib::bam;
use rust_htslib::prelude::*;

pub fn trim_cigar_string_tn5(cigar: &bam::record::CigarString, is_rev: &bool) {


	for mut i in &mut cigar.into_iter() {
		trim_single_cigar(i,0);
	}



}

fn trim_single_cigar(c: &bam::record::Cigar, left_to_trim: u32) {
	use rust_htslib::bam::record::Cigar;
	match c {
		Cigar::Match(u32) | Cigar::Equal(u32) | 
			Cigar::Diff(u32) | Cigar::Ins(u32) |
			Cigar::SoftClip(u32) => (),
		_ => (),
	};
}

//pub consume_cigar()