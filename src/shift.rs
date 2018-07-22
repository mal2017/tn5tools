use rust_htslib::bam;
use rust_htslib::prelude::*;
use std::str;
use ::cigar_utils;

// shifts all reads in a bam

pub fn shift_bam(ib: &str, ob: &str, p: usize) {
	use std::ops::Deref;
	let mut bam = bam::Reader::from_path(ib).unwrap();
	let header = bam::Header::from_template(bam.header());
	let mut obam = bam::Writer::from_path(ob, &header).unwrap();

	// only initialize once + use while loop
	let mut bam_rec = bam::Record::new();

	while let Ok(_r) = bam.read(&mut bam_rec) {

		let insize = bam_rec.insert_size();

		let qname = bam_rec.qname().to_owned();
		let cigar = bam_rec.cigar().deref().clone();
		let qual  = bam_rec.qual().to_owned();
		let seq   = bam_rec.seq().as_bytes();

		cigar_utils::trim_cigar_pos(cigar, 31);

		// https://www.biostars.org/p/76892/
		let new_seq: &[u8] = if bam_rec.is_reverse() {
			let idx = seq.len() - 4;
			&seq[..idx]
		} else {
			&seq[5..]
		};

		let pos = if bam_rec.is_reverse() {
			bam_rec.pos() - 5
		} else {
			bam_rec.pos() + 4
		} as u32;


		// account for offsets with new insert size
		bam_rec.set_insert_size(insize - 9);
		
		bam_rec.set_pos(pos as i32);

		//bam_rec.set(&qname, &cigar, new_seq, &qual);
		obam.write(&bam_rec).unwrap();
	}	

}


pub fn shift_read(r: bam::Record) {
	// TODO: impl for record
}