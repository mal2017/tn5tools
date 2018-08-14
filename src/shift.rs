use rust_htslib::bam;
use rust_htslib::sam;
use rust_htslib::prelude::*;
use std::str;
use ::cigar_utils;

// shifts all reads in a bam

pub fn shift_bam(ib: &str) {
	use std::ops::Deref;
	use rust_htslib::bam::record::{Cigar, CigarString, Aux};

	let mut bam = bam::Reader::from_path(ib).unwrap();
	let header = bam::Header::from_template(bam.header());

	//let mut obam = bam::Writer::from_path(ob, &header).unwrap();
	let mut osam = sam::Writer::from_stdout(&header).unwrap();

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
	let mut pos: u32;
	let mut idx: usize;

	while let Ok(_r) = bam.read(&mut bam_rec) {

		is_rev = bam_rec.is_reverse();
		qual   = bam_rec.qual().to_owned();
		qname  = bam_rec.qname().to_owned();
		seq    = bam_rec.seq().as_bytes();
		slen   = bam_rec.seq().len();

		//cigar_utils::trim_cigar_string_tn5(&cigar, &is_rev);
		// https://www.biostars.org/p/76892/

		if is_rev {
			idx = slen - 4;
			new_seq = seq[..idx].to_vec();
			pos = bam_rec.pos() as u32 - 5;
			new_qual = qual[..idx].to_vec();
			new_insize = bam_rec.insert_size() + 9;
		} else {
			pos = bam_rec.pos() as u32 + 4;
			new_seq = seq[5..].to_vec();
			new_qual = qual[5..].to_vec();
			new_insize = bam_rec.insert_size() - 9;
		}

		new_cigarstr = CigarString(vec![Cigar::Match(new_seq.len() as u32)]);

		// account for offsets with new insert size
		bam_rec.set_insert_size(new_insize);
		
		bam_rec.set_pos(pos as i32);

		bam_rec.set(&qname,
			&new_cigarstr,
			&new_seq,
			&new_qual);

		osam.write(&bam_rec).unwrap();
	}	

}
