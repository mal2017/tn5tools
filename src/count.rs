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

pub fn counts(bed_path: &str, bams: &Vec<&str>, p: usize) {
	use rayon::prelude::*;
	use csv;
	use rust_htslib::bam::IndexedReader;
	use rust_htslib::bam;
	use super::regions::expand_region;
	use super::regions;
	use indicatif::{ProgressBar, HumanDuration};
	use std::sync::{Arc, Mutex};
	use rayon::ThreadPoolBuilder;
	use std::io;
	use bio::io::bed;
	use bio::io::bed::Records;
	use bio::io::bed::Reader;
	use std::fs;
	use std::time::Instant;
	use ndarray;
	use ndarray::prelude::*;
	use serde::ser::{Serialize, Serializer, SerializeStruct};


	let mut reader = bed::Reader::from_file(bed_path).unwrap();

	// vector of regions with expanded coords
	let recs: Vec<bed::Record> = reader.records()
									  .map(|a| a.unwrap())
									  .map(|a| expand_region(a, -5, 5)) // 5 both sides
					                  .collect();

    ThreadPoolBuilder::new().num_threads(p).build_global().unwrap();

    let n_row = recs.len();

    let n_col = bams.len();

    let mut cuts_vec: Vec<u32> = Vec::new();

	for bam in bams {

		//let pb = Arc::new(Mutex::new(ProgressBar::new(n_regions)));

		let idxr = Arc::new(Mutex::new(IndexedReader::from_path(bam).unwrap()));

		let cuts: Vec<u32> = recs.par_iter()
    						  .map(|a| get_count_in_region(&idxr, &a))
    						  .collect();

    	cuts_vec.extend(cuts);

		//pb.lock().unwrap().finish_with_message("Done");   
    }

    let arr = Array::from_shape_vec((n_col,n_row), cuts_vec).unwrap()
    														.reversed_axes();


    let mut csv_header = vec!["region"];

    csv_header.append(&mut bams.clone());

    let mut wtr = csv::Writer::from_writer(io::stdout());
    
    wtr.write_record(&csv_header);


    for i in 0..recs.len() {

    	wtr.write_field(regions::region_as_string(&recs[i]));

    	wtr.serialize(arr.row(i).to_vec());

    }

    wtr.flush();
}