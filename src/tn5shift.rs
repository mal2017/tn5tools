use rust_htslib::bam;
use rust_htslib::prelude::*;
use std::str;
use ::atacbam::header_utils;
use ::atacbam::tn5record;

// Shift all mapped reads in bam and output to new bam.
pub fn tn5shift_bam(ib: &str, ob: &str, p: usize) {
	// todo make sure has mate info in tags	
	let mut bam = bam::Reader::from_path(ib).unwrap();
	let header = header_utils::edit_hdr_unsrt(bam.header());
	let mut obam = bam::Writer::from_path(ob, &header).unwrap();
	
	if p >= 2 {
		let p2 = if (p % 2) == 0 {
			p / 2
		} else {
			(p - 1) / 2
		};
		bam.set_threads(p2).unwrap();
		obam.set_threads(p-p2).unwrap();
	} else {
		bam.set_threads(1).unwrap();
		obam.set_threads(1).unwrap();
	}

	let mut tn5_rec: tn5record::Tn5Record; 
	while let Some(x) = bam.records().next() {
		tn5_rec = tn5record::Tn5Record::from_record(x.unwrap()).unwrap();
		tn5_rec.tn5shift();
		obam.write(&tn5_rec.inner).unwrap();
	}	
}


