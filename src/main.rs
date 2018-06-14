#[macro_use]
extern crate clap;
extern crate tn5tools;
extern crate rust_htslib;
extern crate rayon;
extern crate indicatif;
extern crate serde;
extern crate bio;
#[macro_use(s)]
extern crate ndarray;

#[macro_use]
extern crate utah;

extern crate csv;

use utah::prelude::*;

fn main() {
	use clap::App;

	let yml = load_yaml!("cli.yml");
	let m = App::from_yaml(yml).get_matches();

	
 	if let Some(m) = m.subcommand_matches("counts") {
        // Safe to use unwrap() because of the required() option
        let bams: Vec<&str> = m.values_of("BAMS").unwrap().collect();
    	let bed: &str = m.value_of("BED").unwrap();
    	let threads: u32 =  m.value_of("threads").unwrap_or("1").to_string().parse().unwrap();
        counts(bed, &bams, threads as usize);
	}

}


fn counts(bed_path: &str, bams: &Vec<&str>, p: usize) {
	use rayon::prelude::*;
	use rust_htslib::bam::IndexedReader;
	use rust_htslib::bam;
	use tn5tools::regions::expand_region;
	use tn5tools::regions;
	use tn5tools::count;
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
    						  .map(|a| count::get_count_in_region(&idxr, &a))
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













	// matrix with sizes from min size found to
	// max size found on rows, samples in cols
	// each output is for a region or set of regions