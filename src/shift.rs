use rust_htslib::bam;
use rust_htslib::sam;
use rust_htslib::bam::header::Header;
use rust_htslib::prelude::*;
use std::str;
use ::cigar_utils;
use ::header_utils;
use linear_map::LinearMap;
use rayon::iter::ParallelBridge;
use rayon::prelude::*;
use rayon::ThreadPoolBuilder;
use std::sync::{Arc, Mutex};


// shifts all reads in a bam

pub fn shift_bam(ib: &str, ob: &str, p: usize) {
 
	use std::ops::Deref;
	use rust_htslib::bam::record::{Cigar, CigarString, Aux};
	let mut bam = bam::Reader::from_path(ib).unwrap();

	let header = header_utils::edit_hdr_unsrt(bam.header());

	let mut obam = bam::Writer::from_path(ob, &header).unwrap();

	let mut bam_rec = bam::Record::new();

	if p >= 2 {
		let p2 = if ((p % 2) == 0 ){
			p / 2
		} else {
			(p - 1) / 2
		};
		bam.set_threads(p2);
		obam.set_threads(p-p2);
	} else {
		bam.set_threads(1);
		obam.set_threads(1);
	}

	let mut is_rev: bool;
	let mut m_new_cigarstr;
	let mut new_cigarstr;
	let mut qual;
	let mut qname;
	let mut seq;
	let mut slen;
	let mut new_insize;
	let mut new_qual;
	let mut new_seq;
	let mut pos;
	let mut idx;
	let mut mpos;
	let mut bin;
	
	while let Ok(_r) = bam.read(&mut bam_rec) {
		// TODO: throw out unmapped
		// TODO: throw out improper pairs

		is_rev = bam_rec.is_reverse();

		m_new_cigarstr = cigar_utils::trim_cigar_adjust_shift_tn5(&CigarString::from_bytes(bam_rec.aux(b"MC")
	 																							 .unwrap()
	 																							 .string()).unwrap(), &!is_rev);

	 	new_cigarstr = cigar_utils::trim_cigar_adjust_shift_tn5(&bam_rec.cigar(), &is_rev);
	 	//print!("{:?} -- ", new_cigarstr.cigar);

	 	pos    = bam_rec.pos();  
	 	qual   = bam_rec.qual().to_owned();
	 	qname  = bam_rec.qname().to_owned();
	 	seq    = bam_rec.seq().as_bytes();
	 	slen   = bam_rec.seq().len();
	 	mpos   = bam_rec.mpos();
	 	new_insize = bam_rec.insert_size() - 
	 		(m_new_cigarstr.seq_shift + 
	 			new_cigarstr.seq_shift) as i32;

		// https://www.biostars.org/p/76892/
		// https://jef.works/blog/2017/03/28/CIGAR-strings-for-dummies/

	 	if is_rev {
	 		idx = slen - new_cigarstr.seq_shift as usize;
	 		new_seq = seq[..idx].to_vec();
	 		new_qual = qual[..idx].to_vec();
	 		bam_rec.set_mpos(mpos + m_new_cigarstr.pos_shift as i32);
	 	} else {
	 		pos = bam_rec.pos() + new_cigarstr.pos_shift as i32;
	 		new_seq = seq[(new_cigarstr.pos_shift as usize)..].to_vec();
	 		new_qual = qual[(new_cigarstr.pos_shift as usize)..].to_vec();
	 		bam_rec.set_pos(pos);
	 	}

	 	//println!("{:?} -- {:?}",is_rev, new_seq.len());

	 	bin = reg2bin(pos, pos + new_seq.len() as i32);

	 	bam_rec.set_bin(bin);
		bam_rec.set_insert_size(new_insize);
	 	
	 	bam_rec.set(&qname,
	 		&new_cigarstr.cigar,
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





