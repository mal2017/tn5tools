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

	ThreadPoolBuilder::new().num_threads(p).build_global().unwrap();

	use std::ops::Deref;
	use rust_htslib::bam::record::{Cigar, CigarString, Aux};
	let mut bam = bam::Reader::from_path(ib).unwrap();

	let header = header_utils::edit_hdr_unsrt(bam.header());

	let mut obam = bam::Writer::from_path(ob, &header).unwrap();
	
	// only initialize once + use while loop
	let mut bam_rec = bam::Record::new();

	let mut new_cigarstr = cigar_utils::CigarTrimParams {
		cigar: CigarString(vec![Cigar::Match(0)]),
		seq_shift: 0,
		pos_shift: 0,
	};
	let mut m_new_cigarstr = cigar_utils::CigarTrimParams {
		cigar: CigarString(vec![Cigar::Match(0)]),
		seq_shift: 0,
		pos_shift: 0,
	};
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

	let mut shift_var: u32;


	let a = bam.records().into_iter().par_bridge();

	a.map(|k| () ).collect::<()>();

	
	/*while let Ok(_r) = bam.read(&mut bam_rec) {
		// todo parallize with rayon
		// TODO: avoid error if unpaired/untagged
		if bam_rec.is_unmapped() {
			obam.write(&bam_rec).unwrap();
			continue;
		}

	 	is_rev = bam_rec.is_reverse();
	 	qual   = bam_rec.qual().to_owned();
	 	qname  = bam_rec.qname().to_owned();
	 	seq    = bam_rec.seq().as_bytes();
	 	slen   = bam_rec.seq().len();
	 	new_insize = bam_rec.insert_size() - (m_new_cigarstr.seq_shift + new_cigarstr.seq_shift) as i32;

	 	m_new_cigarstr = cigar_utils::trim_cigar_adjust_shift_tn5(&CigarString::from_bytes(bam_rec.aux(b"MC")
	 																							 .unwrap()
	 																							 .string()).unwrap(),
	 		&!is_rev);


		//TODO delete this
	 	//new_cigarstr = CigarString(vec![Cigar::Match(new_seq.len() as u32)]); // this will safely call all matches
	 	new_cigarstr = cigar_utils::trim_cigar_adjust_shift_tn5(&bam_rec.cigar(), &is_rev);
	 	

		// https://www.biostars.org/p/76892/
		// https://jef.works/blog/2017/03/28/CIGAR-strings-for-dummies/
	 	if is_rev {
	 		idx = slen - new_cigarstr.seq_shift as usize;
	 		new_seq = seq[..idx].to_vec();
	 		pos = bam_rec.pos() - new_cigarstr.pos_shift as i32;
	 		new_qual = qual[..idx].to_vec();
	 		mpos = bam_rec.mpos() + new_cigarstr.pos_shift as i32;
	 	} else {
	 		pos = bam_rec.pos() + new_cigarstr.pos_shift as i32;
	 		new_seq = seq[(new_cigarstr.pos_shift as usize +1 )..].to_vec();
	 		new_qual = qual[(new_cigarstr.pos_shift as usize +1 )..].to_vec();
	 		mpos = bam_rec.mpos() - m_new_cigarstr.pos_shift as i32;
	 	}

	 	bin = reg2bin(pos, pos + new_seq.len() as i32);
	 	bam_rec.set_bin(bin);


	 	bam_rec.set_insert_size(new_insize);
	 	bam_rec.set_pos(pos);
	 	bam_rec.set_mpos(mpos);
	 	bam_rec.set(&qname,
	 		&new_cigarstr.cigar,
	 		&new_seq,
	 		&new_qual);
		obam.write(&bam_rec).unwrap();
	}*/
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





