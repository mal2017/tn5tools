use rust_htslib::bam;
use rust_htslib::prelude::*;
use std::str;
use ::atacbam::header_utils;
use ::atacbam::tn5record;
use std::collections::HashMap;

pub fn fragsizes(ib: &str, p: usize) {
	// todo make sure has mate info in tags	
	let mut bam = bam::Reader::from_path(ib).unwrap();

	bam.set_threads(p).unwrap();

	let mut bam_rec: bam::Record = bam::Record::new(); 
	//let mut sizes: Vec<i32> = Vec::new();
	let mut tracker: HashMap<i32,i32> = HashMap::with_capacity(4000);
	let mut size: i32;
	let mut ct: i32;

	while let Ok(x) = bam.read(&mut bam_rec) {
		if (!bam_rec.is_proper_pair() | bam_rec.is_reverse()) { continue };
		size = bam_rec.insert_size();
		ct = match tracker.get(&size) {
			Some(i) => i.clone(),
			None => 0,
		};

		tracker.insert(size,ct + 1);

	}	
	let min_isize = tracker.keys().min();
	let max_isize = tracker.keys().max();
	println!("{:?}..{:?}", min_isize, max_isize)
	
}


