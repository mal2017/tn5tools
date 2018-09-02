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
	rev        : bool,
	shifted_cigar_config: CigarTrimParams,
	shifted_mate_cigar_config: CigarTrimParams,
	orig_attr: RecordAttributes,
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
			rev: is_rev,
			shifted_cigar_config: get_tn5shift_params(&cigar, &is_rev),
			shifted_mate_cigar_config: get_tn5shift_params(&mate_cigar, &!is_rev),
			orig_attr: ra,
			aux_hm : None,

		})
	}

	pub fn record(&self) -> &rust_htslib::bam::record::Record {
		&self.inner
	}

	pub fn tn5shift(&mut self) {
		if self.rev {
			let idx = self.orig_attr.seq.len() - self.shifted_cigar_config
										   				  .seq_shift as usize;

			self.inner.set_mpos(self.orig_attr.mpos + 
				self.shifted_mate_cigar_config.pos_shift as i32);

			self.inner.set_insert_size(self.orig_attr.insize + 
				(self.shifted_mate_cigar_config.seq_shift + 
					self.shifted_cigar_config.seq_shift) as i32);

			self.inner.set(&self.orig_attr.qname,
						   &self.shifted_cigar_config.cigar,
						   &self.orig_attr.seq[..idx].to_vec(),
						   &self.orig_attr.qual[..idx].to_vec());
			
		} else {
			let idx = self.shifted_cigar_config.seq_shift as usize;

			self.inner.set_pos(self.orig_attr.pos + self.shifted_cigar_config.pos_shift as i32);

			self.inner.set_insert_size(self.orig_attr.insize - 
				(self.shifted_mate_cigar_config.seq_shift + 
					self.shifted_cigar_config.seq_shift) as i32);

			self.inner.set(&self.orig_attr.qname,
						   &self.shifted_cigar_config.cigar,
						   &self.orig_attr.seq[idx..].to_vec(),
						   &self.orig_attr.qual[idx..].to_vec());
		}
		let new_pos = self.inner.pos();
		let new_seqlen =  self.inner.seq().len();
		self.inner.set_bin(reg2bin(new_pos, new_pos + new_seqlen as i32));

	}

	fn update_aux_in_record(&self) {

	}
}

#[derive(Debug)]
struct RecordAttributes {
	qual: Vec<u8>,
	insize: i32,
	pos: i32,
	qname: Vec<u8>,
	seq: Vec<u8>,
	mpos: i32,
}

impl RecordAttributes {
	fn from_record(r: &rust_htslib::bam::Record) -> Self {
		RecordAttributes {
			qual: r.qual().to_owned(),
			insize: r.insert_size(),
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

// For details, see SAM V1 format spec.
// Calculates bin for bai style bin scheme.
fn reg2bin(beg: i32, e: i32) -> u16 {
	let end = e - 1;
	if beg>>14 == end>>14 {
		return ((1<<15)-1)/7 + (beg>>14) as u16
	};
	if beg>>17 == end>>17 {
		return ((1<<12)-1)/7 + (beg>>17) as u16
	};
	if beg>>20 == end>>20 {
		return ((1<<9)-1)/7 + (beg>>20) as u16
	};
	if beg>>23 == end>>23 {
		return ((1<<6)-1)/7 + (beg>>23) as u16
	};
	if beg>>26 == end>>26 {
		return ((1<<3)-1)/7 + (beg>>26) as u16
	};
	return 0u16
}