use rust_htslib::bam;
use rust_htslib::prelude::*;
use std::str;
use ::cigar_utils;

// shifts all reads in a bam

pub fn shift_bam(ib: &str, ob: &str, p: usize) {
	use std::ops::Deref;
	use rust_htslib::bam::record::{Cigar, CigarString};

	let mut bam = bam::Reader::from_path(ib).unwrap();
	let header = bam::Header::from_template(bam.header());
	let mut obam = bam::Writer::from_path(ob, &header).unwrap();

	if p > 1 {
		obam.set_threads(p);
	}

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

	while let Ok(_r) = bam.read(&mut bam_rec) {

		is_rev = bam_rec.is_reverse();
		if is_rev {
			new_insize = bam_rec.insert_size() + 9
		} else {
			new_insize = bam_rec.insert_size() - 9
		};

		qual   = bam_rec.qual().to_owned();
		qname  = bam_rec.qname().to_owned();
		seq    = bam_rec.seq().as_bytes();
		slen    = bam_rec.seq().len();

		//cigar_utils::trim_cigar_string_tn5(&cigar, &is_rev);

		// https://www.biostars.org/p/76892/
		let new_seq: &[u8] = if is_rev {
			let idx = slen - 4;
			&seq[..idx]
		} else {
			&seq[5..]
		};

		let new_qual: &[u8] = if is_rev {
			let idx = slen - 4;
			&qual[..idx]
		} else {
			&qual[5..]
		};

		let pos = if is_rev {
			bam_rec.pos() - 5
		} else {
			bam_rec.pos() + 4
		} as u32;

		new_cigarstr = CigarString(vec![Cigar::Match(new_seq.len() as u32)]);

		bam_rec.set(&qname,
					&new_cigarstr,
					&new_seq,
					&new_qual);

		// account for offsets with new insert size
		bam_rec.set_insert_size(new_insize);
		
		bam_rec.set_pos(pos as i32);

		//bam_rec.set(&qname, &cigar, new_seq, &qual);
		obam.write(&bam_rec).unwrap();
	}	

}


pub fn shift_read(r: bam::Record) {
	// TODO: impl for record
}