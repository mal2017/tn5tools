use rust_htslib::bam;
use rust_htslib::sam;
use rust_htslib::bam::header::Header;
use rust_htslib::prelude::*;
use std::str;
use ::cigar_utils;
use ::header_utils;
use linear_map::LinearMap;

// shifts all reads in a bam

pub fn shift_bam(ib: &str, ob: &str) {
	use std::ops::Deref;
	use rust_htslib::bam::record::{Cigar, CigarString, Aux};

	let mut bam = bam::Reader::from_path(ib).unwrap();
	let mut header = header_utils::edit_hdr_unsrt(bam.header());

	//let h3 = header_utils::from_hashmap(&h2);

	let mut obam = bam::Writer::from_path(ob, &header).unwrap();
	//let mut osam = sam::Writer::from_stdout(&header).unwrap();

	// only initialize once + use while loop
	let mut bam_rec = bam::Record::new();
	let mut new_cigarstr = CigarString(vec![Cigar::Match(0)]);
	let mut new_insize: i32;
	let mut is_rev: bool;
	let mut seq: Vec<u8>;
	let mut slen: usize;
	let mut qual: Vec<u8>;
	let mut qname: Vec<u8>;
	let mut new_qual: Vec<u8>;
	let mut new_seq: Vec<u8>;
	let mut pos: i32;
	let mut idx: usize;
	let mut mpos: i32;
	let mut bin: u16;

	while let Ok(_r) = bam.read(&mut bam_rec) {
		// TODO: if is not mapped, write record and skip
	 	is_rev = bam_rec.is_reverse();
	 	qual   = bam_rec.qual().to_owned();
	 	qname  = bam_rec.qname().to_owned();
	 	seq    = bam_rec.seq().as_bytes();
	 	slen   = bam_rec.seq().len();

		// 	//cigar_utils::trim_cigar_string_tn5(&cigar, &is_rev);
		// 	// https://www.biostars.org/p/76892/

	 	if is_rev {
	 		idx = slen - 4;
	 		new_seq = seq[..idx].to_vec();
	 		pos = bam_rec.pos() as i32 - 5;
	 		new_qual = qual[..idx].to_vec();
	 		new_insize = bam_rec.insert_size() + 9;
	 		mpos = bam_rec.mpos() as i32 + 4;
	 	} else {
	 		pos = bam_rec.pos() as i32 + 4;
	 		new_seq = seq[5..].to_vec();
	 		new_qual = qual[5..].to_vec();
	 		new_insize = bam_rec.insert_size() - 9;
	 		mpos = bam_rec.mpos() as i32 - 5;
	 	}

	 	bin = reg2bin(pos, pos + new_seq.len() as i32);

	 	if bin != bam_rec.bin() {
	 		println!("{:?} : {:?} : {:?}",pos, bin, bam_rec.bin());
	 	}

	 	bam_rec.set_bin(bin);

	 	new_cigarstr = CigarString(vec![Cigar::Match(new_seq.len() as u32)]);

	 	// account for offsets with new insert size
	 	bam_rec.set_insert_size(new_insize);
	
	 	bam_rec.set_pos(pos);
	 	bam_rec.set_mpos(mpos);
	
	 	bam_rec.set(&qname,
	 		&new_cigarstr,
	 		&new_seq,
	 		&new_qual);

		obam.write(&bam_rec).unwrap();
	}	

}

// For details, see SAM V1 format spec.
// Calculates bin for bai style bin scheme.
fn reg2bin(beg: i32, e: i32) -> u16 {
	let end = e - 1;
	if beg>>14 == end>>14 {
		return(((1<<15)-1)/7 + (beg>>14) as u16)
	};
	if beg>>17 == end>>17 {
		return(((1<<12)-1)/7 + (beg>>17) as u16)
	};
	if beg>>20 == end>>20 {
		return(((1<<9)-1)/7 + (beg>>20) as u16)
	};
	if beg>>23 == end>>23 {
		return(((1<<6)-1)/7 + (beg>>23) as u16)
	};
	if beg>>26 == end>>26 {
		return(((1<<3)-1)/7 + (beg>>26) as u16)
	};
	return(0u16)
}





