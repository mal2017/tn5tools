use rust_htslib::bam::record::{CigarString, Aux};
use rust_htslib;
use ::atacbam::cigar_utils::*;
use std::collections::HashMap;

// Useful info: https://gist.github.com/ChrisWellsWood/84421854794037e760808d5d97d21421


trait Empty<T> {
	fn empty() -> T;
}

#[derive(Debug)]
pub struct Tn5Record {
	inner      : rust_htslib::bam::Record,
	mate_cigar : CigarString,
	is_shifted : bool,
	rev        : bool,
	intended_shift: u32,
	shifted_cig  : Option<CigarString>,
	shifted_mate_cig : Option<CigarString>,
	pos_shift  : Option<u32>,
	seq_shift  : Option<u32>,
	pos_mate_shift : Option<u32>,
	seq_mate_shift : Option<u32>,
	aux_hm : Option<HashMap<String,String>>,
}

impl Tn5Record {
	pub fn from_record(r: rust_htslib::bam::Record) -> Result<Self, Tn5RecordError> {
		let is_rev = r.is_reverse();

		let mate_cig = CigarString::from_bytes(r.aux(b"MC").unwrap().string()).unwrap();

		Ok(Tn5Record{
			inner: r,
			mate_cigar : mate_cig,
			is_shifted: false,
			rev: is_rev,
			intended_shift: if is_rev {4} else {5},
			shifted_cig: None,
			shifted_mate_cig: None,
			pos_shift: None,
			seq_shift: None,
			pos_mate_shift: None,
			seq_mate_shift: None,
			aux_hm : None,

		})
	}

	pub fn record(&self) -> &rust_htslib::bam::record::Record {
		&self.inner
	}

	fn set_cigar_shift_tn5(&mut self) {
		let cigar_config = get_tn5shift_params(&self.inner.cigar(), &self.rev);
		self.shifted_cig = Some(cigar_config.cigar);
		self.pos_shift = Some(cigar_config.pos_shift);
		self.seq_shift = Some(cigar_config.seq_shift);
	}

	fn set_mate_cigar_shift_tn5(&mut self) {
		let cigar_config = get_tn5shift_params(&self.mate_cigar, &!self.rev);
		self.shifted_mate_cig = Some(cigar_config.cigar);
		self.pos_mate_shift = Some(cigar_config.pos_shift);
		self.seq_mate_shift = Some(cigar_config.seq_shift);
	}

	pub fn is_tn5shifted(&self) -> bool {
		self.is_shifted
	}

	pub fn tn5shift(&self) {
		// TODO
	}

	fn tn5_seq_shift(&self) {
		// TODO
	} 

	fn tn5_pos_shift(&self) {
		// TODO
	}

	fn tn5_qual_shift(&self) {
		// TODO
	}

	fn tn5_new_insize(&self) {
		//TODO
	}

	fn set_new_bin(&self) {
		// TODO
	}

	fn update_cigar_in_record(&self) {

	}

	fn update_aux_in_record(&self) {

	}


}

quick_error! {
    #[derive(Debug, Clone)]
    pub enum Tn5RecordError {
        InvalidRecord {
            description("Record is unmapped/unpaired")
        }
    }
}