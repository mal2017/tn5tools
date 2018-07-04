use rust_htslib::bam;
use rust_htslib::prelude::*;


// shifts all reads in a bam

pub fn shift_bam(ib: &str, ob: &str, p: usize) {
	let mut bam = bam::Reader::from_path(ib).unwrap();
	let header = bam::Header::from_template(bam.header());
	let mut obam = bam::Writer::from_path(ob, &header).unwrap();

	let mut bam_rec = bam::Record::new();

	while let Ok(_r) = bam.read(&mut bam_rec) {
		let pos = if bam_rec.is_reverse() {
			bam_rec.pos() - 4
		} else {
			bam_rec.pos() + 4 // bam is 0based beds, then plus 4 for shift
			//TODO: check on left-rightness of bam records
			//TODO: check if bam/bed zero base assumptions are correct
		} as u32;
		bam_rec.set_pos(pos as i32);
		obam.write(&bam_rec).unwrap();
	}	

}


pub fn shift_read(r: &bam::Record) {
	// TODO: impl for record
}