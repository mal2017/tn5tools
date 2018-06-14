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
extern crate csv;


fn main() {
	use clap::App;
	use tn5tools::count;

	let yml = load_yaml!("cli.yml");
	let m = App::from_yaml(yml).get_matches();

	
 	if let Some(m) = m.subcommand_matches("counts") {
        // Safe to use unwrap() because of the required() option
        let bams: Vec<&str> = m.values_of("BAMS").unwrap().collect();
    	let bed: &str = m.value_of("BED").unwrap();
    	let threads: u32 =  m.value_of("threads").unwrap_or("1").to_string().parse().unwrap();
        count::counts(bed, &bams, threads as usize);
	}

}
















	// matrix with sizes from min size found to
	// max size found on rows, samples in cols
	// each output is for a region or set of regions