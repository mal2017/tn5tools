use bio;
use rust_htslib::bam::IndexedReader;
use rust_htslib::bam::Records;
use rust_htslib::bam::Record;
use std::sync::{Arc, Mutex};
use std;
use rust_htslib::bam;
use rust_htslib::sam;
use rust_htslib::prelude::*;

pub fn get_reads_in_region(idxr: &mut IndexedReader) {
	let chrom_as_bytes = "chr4".as_bytes();
	let tid = idxr.header().tid(chrom_as_bytes).unwrap();
	idxr.fetch(tid, 45500, 50000);
	let mut rec = Record::new();

	while let Ok(r) = idxr.read(&mut rec) {
		println!("{:?}", rec);
	}
}

// -> bam::Records<'_,IndexedReader>