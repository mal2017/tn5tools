#[macro_use]
extern crate clap;
extern crate tn5tools;
extern crate rust_htslib;
extern crate rayon;

fn main() {
	use clap::App;

	let yml = load_yaml!("cli.yml");
	let m = App::from_yaml(yml).get_matches();

	
 	if let Some(m) = m.subcommand_matches("counts") {
        // Safe to use unwrap() because of the required() option
        let bams: Vec<&str> = m.values_of("BAMS").unwrap().collect();
    	let bed: &str = m.value_of("BED").unwrap();
        counts(bed, &bams);
	}

}


fn counts(bed_path: &str, bams: &Vec<&str>) {
	use rayon::prelude::*;
	use rust_htslib::bam::IndexedReader;
	use tn5tools::regions::get_reads_in_region;
	use tn5tools::regions;

	let regions = regions::bed_as_strings(bed_path);

	let n_regions = regions.len() as u32;

	for bam in bams {
		let counter = regions.par_iter()
							 .map(|x| get_reads_in_region(bam, &x.0, &x.1, &x.2));

		let counts: Vec<u32> = counter.collect();     
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