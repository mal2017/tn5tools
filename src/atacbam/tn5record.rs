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
	pub inner  : rust_htslib::bam::Record,
	is_shifted : bool,
	rev        : bool,
	intended_shift: u32,
	shifted_cigar_config: CigarTrimParams,
	shifted_mate_cigar_config: CigarTrimParams,
	attributes: RecordAttributes,
	aux_hm : Option<HashMap<String,String>>,
}

impl Tn5Record {
	pub fn from_record(r: rust_htslib::bam::Record) -> Result<Self, Tn5RecordError> {
		let is_rev = r.is_reverse();
		let cigar = r.cigar();
		let mate_cigar = CigarString::from_bytes(r.aux(b"MC").unwrap().string()).unwrap();
		let ra = RecordAttributes::from_record(&r);
		Ok(Tn5Record{
			inner: r,
			is_shifted: false,
			rev: is_rev,
			intended_shift: if is_rev {5} else {4},
			shifted_cigar_config: get_tn5shift_params(&cigar, &is_rev),
			shifted_mate_cigar_config: get_tn5shift_params(&mate_cigar, &is_rev),
			attributes: ra,
			aux_hm : None,

		})
	}

	pub fn record(&self) -> &rust_htslib::bam::record::Record {
		&self.inner
	}

	fn is_tn5shifted(&self) -> bool {
		self.is_shifted
	}


	pub fn tn5shift(&self) {
		// TODO
	}

	fn tn5_seq_shift(&mut self) {

		if self.rev {
			
		}
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

#[derive(Debug)]
struct RecordAttributes {
	qual: Vec<u8>,
	pos: i32,
	qname: Vec<u8>,
	seq: Vec<u8>,
	mpos: i32,
}

impl RecordAttributes {
	fn from_record(r: &rust_htslib::bam::Record) -> Self {
		RecordAttributes {
			qual: r.qual().to_owned(),
			pos: r.pos(),
			qname: r.qname().to_owned(),
			seq: r.seq().as_bytes(),
			mpos: r.mpos(),
		}
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