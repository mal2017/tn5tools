use bio;
use rust_htslib::bam::IndexedReader;
use std::sync::{Arc, Mutex};
use std;
use rust_htslib::bam;
use rust_htslib::sam;
use rust_htslib::prelude::*;

pub fn get_reads_in_region(idxr: &mut IndexedReader)  {
	let chrom_as_bytes = "chr4".as_bytes();
	let tid = idxr.header().tid(chrom_as_bytes).unwrap();
	idxr.fetch(tid, 1000, 50000);
	idxr.records();
}

//-> bam::Records<'_,IndexedReader>