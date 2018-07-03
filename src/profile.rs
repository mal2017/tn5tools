use bio::io::bed;
use rust_htslib::bam::IndexedReader;
use std::sync::{Arc, Mutex};

use rust_htslib::prelude::*;
use rust_htslib::bam::buffer::RecordBuffer;

pub fn get_profile_in_region(idxr: &mut IndexedReader, rec: &bed::Record) -> u32 {

	use rust_htslib::bam;
		
	let chrom_as_bytes = rec.chrom().as_bytes();

	let tid = idxr.header().tid(chrom_as_bytes).unwrap();
	
	idxr.fetch(tid, rec.start() as u32, rec.end() as u32).unwrap();

	let mut count = 0;

	let recs = idxr.records();

	for r in recs {
		r.unwrap().set_pos(1);
		println!("{:?}", r);
	}
	count as u32
	//TODO make sure im not shifting reads twice
}


pub fn profile(bed_path: &str, bam: &str, p: usize) {
	//&Vec<&str>
	use rayon::prelude::*;
	use csv;
	use rust_htslib::bam::IndexedReader;
	use super::regions::expand_region;
	use super::regions;
	use indicatif::{ProgressBar, ProgressStyle};
	use rayon::ThreadPoolBuilder;
	use std::io;
	use bio::io::bed;
	use ndarray::prelude::*;

	let mut reader = bed::Reader::from_file(bed_path).unwrap();

	let mut idxr = IndexedReader::from_path(bam).unwrap();

	// vector of regions with expanded coords
	let recs: Vec<bed::Record> = reader.records()
									   .map(|a| a.unwrap())
									   .map(|a| expand_region(a, -5, 5))
									   .collect();


	let cuts: Vec<u32> = recs.iter()
							.map(|a| { get_profile_in_region(&mut idxr, &a)})
							.collect();


	//println!("{:?}", cuts);
									   //.par_iter()
									   //.unwrap();

	

	//let res: Vec<bed::Record> = recs.iter()
	//								  .map(|a| expand_region(a, -5, 5)) // 5 both sides
	//				                  .collect();


	//let cuts: Vec<u32> = recs.par_iter()
    //						  .map(|a| { pb.inc(1); get_count_in_region(&idxr, &a)})
    //						  .collect();


    //let arr = Array::from_shape_vec((n_col as usize, n_row as usize), cuts_vec).unwrap()
    //														.reversed_axes();

    //let mut csv_header = vec!["region"];

    //csv_header.append(&mut bams.clone());

    //let mut wtr = csv::Writer::from_writer(io::stdout());
    
    //wtr.write_record(&csv_header).unwrap();


    //for (i, c) in recs.iter().enumerate() {

    //	wtr.write_field(regions::region_as_string(c)).unwrap();
    	//&recs[i]
    //	wtr.serialize(arr.row(i).to_vec()).unwrap();

    //}

    //wtr.flush().unwrap();
}