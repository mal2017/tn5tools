use bio::io::bed;
use rust_htslib::bam::IndexedReader;
use std::sync::{Arc, Mutex};

use rust_htslib::prelude::*;

pub fn get_count_in_region(idxr: &Arc<Mutex<IndexedReader>>, rec: &bed::Record) -> u32 {
	use rust_htslib::bam;
	let chrom_as_bytes = rec.chrom().as_bytes();

	let mut idxr = idxr.lock().unwrap();
	
	let tid = idxr.header().tid(chrom_as_bytes).unwrap();
	
	idxr.fetch(tid, rec.start() as u32, rec.end() as u32).unwrap();
	
	let mut bam_rec = bam::Record::new();
	
	let mut count = 0;

	while let Ok(_r) = idxr.read(&mut bam_rec) {
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
	use super::regions::expand_region;
	use super::regions;
	use indicatif::{ProgressBar, ProgressStyle};
	use rayon::ThreadPoolBuilder;
	use std::io;
	use bio::io::bed;
	use ndarray::prelude::*;

	let mut reader = bed::Reader::from_file(bed_path).unwrap();

	// vector of regions with expanded coords
	let recs: Vec<bed::Record> = reader.records()
									  .map(|a| a.unwrap())
									  .map(|a| expand_region(a, -5, 5)) // 5 both sides
					                  .collect();

    

    let n_row = recs.len() as u64;

    let n_col = bams.len() as u64;

    ThreadPoolBuilder::new().num_threads(p).build_global().unwrap();
    
    //------------------------------------	
    
    // BASIC PARALLEL VERSION
    /*let mut cuts_vec: Vec<u32> = Vec::new();

	for bam in bams {

		let pb = ProgressBar::new(n_row);

		pb.set_style(ProgressStyle::default_bar()
    	  .template("[{eta_precise}] {bar:40.red/blue} {pos:>7}/{len:7} {msg}")
    	  .progress_chars("$$-"));

		let idxr = Arc::new(Mutex::new(IndexedReader::from_path(bam).unwrap()));

		let cuts: Vec<u32> = recs.par_iter()
    						  .map(|a| { pb.inc(1); get_count_in_region(&idxr, &a)})
    						  .collect();

    	cuts_vec.extend(cuts);

		pb.finish_with_message("Cash Money!!");   
    }

    let arr = Array::from_shape_vec((n_col as usize, n_row as usize), cuts_vec).unwrap()
    														.reversed_axes();*/

    //------------------------------------	
    // MEGA PARALLEL	
    let pb = ProgressBar::new(n_row * n_col);

	pb.set_style(ProgressStyle::default_bar()
    	  .template("Counting transposition events... {bar:40.blue/red} {percent}% {msg}")
    	  .progress_chars("##-"));

	// TODO make this a flat matrix and elim the flattening below
    let cuts_vec: Vec<Vec<u32>> = bams.par_iter().map(|bam| {

		let idxr = Arc::new(Mutex::new(IndexedReader::from_path(bam).unwrap()));

		let cuts: Vec<u32> = recs.par_iter()
    						.map(|a| { pb.inc(1); get_count_in_region(&idxr, &a)})
    						.collect();
		
		cuts
    }).collect();
    pb.finish_with_message("Complete.");


    let cuts_vec_flat: Vec<u32> = cuts_vec.iter()
    				  .flat_map(|a| a.iter())
    				  .cloned()
    				  .collect();

    let arr = Array::from_shape_vec((n_col as usize, n_row as usize), cuts_vec_flat).unwrap()
    														.reversed_axes();

    //------------------------------------


    let mut csv_header = vec!["region"];

    csv_header.append(&mut bams.clone());

    let mut wtr = csv::Writer::from_writer(io::stdout());
    
    wtr.write_record(&csv_header).unwrap();


    for i in 0..recs.len() {

    	wtr.write_field(regions::region_as_string(&recs[i])).unwrap();

    	wtr.serialize(arr.row(i).to_vec()).unwrap();

    }

    wtr.flush().unwrap();
}