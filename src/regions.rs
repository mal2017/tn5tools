//use std::path::Path;
use bio;


// shifts the coords of each region
// TODO: safety checks to prevent going out of bounds for a given chrom
// 		This would be higher priority if we were to release, but unlikely that
//		any region we would look at would be within 5nt of the start or end of
//		a contig
pub fn expand_region(mut r: bio::io::bed::Record,
	start_shift: i64, end_shift: i64) -> bio::io::bed::Record {
	let new_start = r.start() as i64 + start_shift;
	let new_end = r.end() as i64 + end_shift;
	r.set_start(new_start as u64);
	r.set_end(new_end as u64);
	r
}

pub fn region_as_string(b: &bio::io::bed::Record) -> String {
	// todo check on 0-1 based stuff
	format!("{}:{}-{}",b.chrom(),b.start(),b.end())
}
