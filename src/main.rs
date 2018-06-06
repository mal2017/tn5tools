#[macro_use]
extern crate clap;
extern crate tn5tools;
extern crate rust_htslib;
extern crate rayon;
extern crate indicatif;
extern crate csv;

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
	use tn5tools::regions::get_reads_in_region;
	use tn5tools::regions;
	use indicatif::{ProgressBar, HumanDuration};
	use std::sync::{Arc, Mutex};
	use rayon::ThreadPoolBuilder;
	let regions = regions::bed_as_strings(bed_path);

	let n_regions = regions.len() as u64;
	ThreadPoolBuilder::new().num_threads(p).build_global().unwrap();
	for bam in bams {
		let pb = Arc::new(Mutex::new(ProgressBar::new(n_regions)));

		let counter = regions.par_iter()
							 .map(|x| {		 	
							 	get_reads_in_region(bam, &x.0, &x.1, &x.2)
							 });

		let counts: Vec<u32> = counter.collect();
		pb.lock().unwrap().finish_with_message("Done");   
    }
    
}













	// matrix with sizes from min size found to
	// max size found on rows, samples in cols
	// each output is for a region or set of regions



// make iterator of the records we're interested in
// make chain iterator for 
//    shifting and creating bed records/interval tree 
//	  counting these records as specified

// for profile, normalize to a percentage of total reads
// per region or per base?? or just raw counts and worry about post
// processing later

// make test bam - by chrom?