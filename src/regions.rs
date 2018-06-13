use csv;
use std::convert::AsRef;
use std::fs;
use std::io;
use std::path::Path;
use bio;


// shifts the coords of each region
// TODO: safety checks to prevent going out of bounds for a given chrom
// 		This would be higher priority if we were to release, but unlikely that
//		any region we would look at would be within 5nt of the start or end of
//		a contig
pub fn expand_region(mut r: bio::io::bed::Record,
	start_shift: i64, end_shift: i64) -> bio::io::bed::Record {
	let new_start = r.start() as i64 + start_shift;
	let new_end = r.start() as i64 + end_shift;
	r.set_start(new_start as u64);
	r.set_end(new_end as u64);
	r
}










pub fn bed_as_strings(bf: &str) -> Vec<(String, u32, u32)> {
	use bio::io::bed;
	use bio::io::bed::Records;
	use bio::io::bed::Reader;
	use std::fs;
	use indicatif::{ProgressBar, HumanDuration};
	use std::time::Instant;

	let mut reader = bed::Reader::from_file(bf).unwrap();

	let mut regions_unwrapped: Vec<(String, u32, u32)> = Vec::new();
	let spin = ProgressBar::new_spinner();
	let begin = Instant::now();
	for r in reader.records() {
		let rec = r.unwrap();
		regions_unwrapped.push((rec.chrom().to_string(), 
			rec.start() as u32, 
			rec.end() as u32));
		spin.tick();
	}

	let time = format!("BED imported in {}",HumanDuration(begin.elapsed()));

	spin.finish_with_message(&time);
	regions_unwrapped
}

pub fn get_reads_in_region(bam: &str, chrom: &str, start: &u32, end: &u32) -> u32 {
	use rust_htslib::bam;
	use rust_htslib::sam;
	use rust_htslib::prelude::*;
	use rust_htslib::bam::IndexedReader;
	use bio::io::bed::Record as BedRecord;

	let mut idxr = IndexedReader::from_path(bam).unwrap();

	let chrom_as_bytes = chrom.as_bytes();
	let tid = idxr.header().tid(chrom_as_bytes).unwrap();

	idxr.fetch(tid, *start -5 , *end + 5 );

	let mut count = 0;
	for s in idxr.records() {
		let mut record = s.unwrap();

		let pos = if record.is_reverse() {
			record.pos() + 1 - 5
		} else {
			record.pos() + 1 + 4 // bam is 0based beds, then plus 4 for shift
		} as u32;

		if pos >= *start && pos <= *end {
			count += 1; // remember to write a condition here, just a test for now.
		}
	}
	count
}