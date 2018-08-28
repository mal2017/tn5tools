use rust_htslib::bam;
use rust_htslib::prelude::*;
use std::iter;
use rust_htslib::bam::record::{Cigar, CigarString, Aux};
use rayon::prelude::*;

// 1. expand cigar
// 2. consume cigar while following rules
// 3. as this happens, build read specific rule for sequence


fn expand_cigarstring(cs: &bam::record::CigarString) -> Vec<char> {
	cs.iter()
	  .flat_map(|a| vec![a.char(); a.len() as usize])
	  .collect::<Vec<char>>()
}


// TODO char_vec_to_cigarstring
fn char_vec_to_cigarstring(cv: Vec<char>) -> CigarString {
	let mut current_char = 'z'; // dummy char
	let mut current_ct = 1;

	let mut cigar_vec: Vec<Cigar> = Vec::new();

	for x in cv {
		if current_char == 'z' {
			current_char = x;
		} else if current_char == x {
			current_ct += 1;
		} else {
			cigar_vec.push(char_to_cigar(&current_char, current_ct).unwrap());
			current_char = x;
			current_ct = 1;
		}
		cigar_vec.push(char_to_cigar(&current_char, current_ct).unwrap()); // final cigar
	}

	CigarString(cigar_vec)
}

fn char_to_cigar(t: &char, n: u32) -> Option<Cigar> {
	match t {
		'M' => Some(Cigar::Match(n)),
        'I' => Some(Cigar::Ins(n)),
        'D' => Some(Cigar::Del(n)),
        'N' => Some(Cigar::RefSkip(n)),
        'H' => Some(Cigar::HardClip(n)),
        'S' => Some(Cigar::SoftClip(n)),
        'P' => Some(Cigar::Pad(n)),
        '=' => Some(Cigar::Equal(n)),
        'X' => Some(Cigar::Diff(n)),
        _ => {
        	println!("{:?}", t);
        	None
        },
	}
}

pub fn trim_cigar_adjust_shift_tn5(cs: &bam::record::CigarString, is_rev: &bool) -> CigarTrimParams {

	// current cigar string as char vec
	let cs_vec: Vec<char> = if *is_rev {
		expand_cigarstring(cs).into_iter().rev().collect()
	} else {
		expand_cigarstring(cs)
	};

	// counter for remaining moves through the cigar string
	let mut moves = if *is_rev { 
		5
	} else {
		4
	};
	// tracker for adjusted shift count for the sequence (returned by this func)
	let mut seq_shift = 0; 
	// tracker for adjusted shift count for the reference (returned by this func)
	let mut pos_shift = 0;

	// index tracker for iterating through cigar string
	let mut idx = 0;
	let mut cgr: char;

	// KEEP FOR DEBUG
	//let mut has_indel: bool = false;

	while (moves > 0) & ((cs_vec[idx] != 'D') & (cs_vec[idx] != 'N') & (cs_vec[idx] != 'H')) {
		cgr = cs_vec[idx];
		match cgr {
			'M' | '=' | 'X' => {
				moves -= 1;
				seq_shift += 1;
				pos_shift += 1;
				//has_indel =false;
			},
			'D' | 'N' | 'H' | 'P' => {
				// moves stay the same, essentially ignoring dels for the purpose of shifting
				// seq_shift stays the same, because we ignored dels
				pos_shift += 1; // position shift moves so we can pick back up on the other side of the del
				//has_indel = true;
			},
			'I' | 'S' => {
				moves -= 1;
				seq_shift +=1;
				// pos shift remains, we're moving through an insertion not in ref
				//has_indel = true
			},
			_ => {
				moves -= 1;
				seq_shift += 1;
				pos_shift += 1;
			},
		}
		idx += 1;
		
	}

	let mut new_cs_vec: Vec<char> = cs_vec[idx..].to_owned();

	if *is_rev {
		new_cs_vec = new_cs_vec.into_iter().rev().collect();
	}

	let new_cigar_string = char_vec_to_cigarstring(new_cs_vec);
	
	// KEEP FOR DEBUG
	//if has_indel {
	//	println!("is rev:{:?} ---- seq shift:{:?} ---- pos shift:{:?} ---- cigarstr:{:?}",is_rev, seq_shift,pos_shift,cs);
	//}
	

	CigarTrimParams {
		cigar: new_cigar_string,
		seq_shift: seq_shift,
		pos_shift: pos_shift,
	}

	//CigarString(vec![Cigar::Match(35 as u32)])
}


pub struct CigarTrimParams {
	pub cigar: CigarString,
	pub seq_shift: u32,
	pub pos_shift: u32,
}