use bio;
use bio::io::bed;
use rust_htslib::bam::IndexedReader;
use rust_htslib::bam::Records;
use rust_htslib::bam::Record;
use std::sync::{Arc, Mutex};
use std;
use rust_htslib::bam;
use rust_htslib::sam;
use rust_htslib::prelude::*;

pub fn get_count_in_region(idxr: &Arc<Mutex<IndexedReader>>, rec: &bed::Record) -> u32 {

	let chrom_as_bytes = rec.chrom().as_bytes();

	let mut idxr = idxr.lock().unwrap();
	
	let tid = idxr.header().tid(chrom_as_bytes).unwrap();
	idxr.fetch(tid, rec.start() as u32, rec.end() as u32);
	
	let mut bam_rec = Record::new();
	
	let mut count = 0;

	while let Ok(r) = idxr.read(&mut bam_rec) {
		let pos = if bam_rec.is_reverse() {
			bam_rec.pos() - 4
		} else {
			bam_rec.pos() + 4 // bam is 0based beds, then plus 4 for shift
			//TODO: check on left-rightness of bam records
			//TODO: check if bam/bed zero base assumptions are correct
		} as u32;

		if pos >= rec.start() as u32 && pos <= rec.end() as u32 {
			count += 1;
		}
	}
	count as u32
}

// -> bam::Records<'_,IndexedReader>