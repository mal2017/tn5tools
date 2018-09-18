use rust_htslib::bam;
use rust_htslib::prelude::*;
use std::str;
use ::atacbam::header_utils;
use ::atacbam::tn5record;
use std::collections::HashMap;
use gnuplot::{Figure, Caption, Color};
use ndarray::{Array};
use rayon::ThreadPoolBuilder;
use rayon::prelude::*;

pub fn fragsizes(ibs: &Vec<&str>, p: usize) {
	ThreadPoolBuilder::new().num_threads(p).build_global().unwrap();

	let mut sizes: Vec<Vec<i32>> = ibs.par_iter()
			 .map(|a| get_fragsizes(a))
			 .collect();

	let max = sizes.par_iter()
				   .map(|a| a.iter().count())
				   .max()
				   .unwrap();
							   
	for i in &mut sizes {
		while i.len() < max {
			i.push(0);
		}
	}

	let sz_flat: Vec<i32> = sizes.into_iter().flat_map(|a| a.into_iter()).collect();

	let arr = Array::from_shape_vec((ibs.len(),max), sz_flat)
									.unwrap()
									.reversed_axes();

	let mut fg = Figure::new();
	fg.axes2d()
	  .lines( 0..max, arr.slice(s![..,1]), &[Caption("Proper pairs"), Color("red")]);
	fg.set_terminal("png", "out.png");
	fg.show();
}


/// Takes a bam path and returns the fragment sizes 
/// for all proper pairs.
fn get_fragsizes(ib: &str) -> Vec<i32> {
	let mut bam = bam::Reader::from_path(ib).unwrap();
	let mut bam_rec: bam::Record = bam::Record::new(); 
	let mut tracker: HashMap<i32,i32> = HashMap::with_capacity(5000);
	let mut size: i32;
	let mut ct: i32;
	while let Ok(x) = bam.read(&mut bam_rec) {
		if (!bam_rec.is_proper_pair() | bam_rec.is_reverse()) { continue };
		size = bam_rec.insert_size();
		ct = match tracker.get(&size) {
			Some(i) => i.clone(),
			None => 0,
		};
		tracker.insert(size,ct + 1);
	}	
	fragsize_hashmap_to_vec(&tracker)
}


/// Takes a reference to hashmap with sizes as keys and 
/// counts as values and returns a vector with n elements 
/// where n is the size of the largest fragment found.
fn fragsize_hashmap_to_vec(h: &HashMap<i32,i32>) -> Vec<i32> {
	let max_isize = h.keys().max().unwrap().clone();
	let mut counts = vec![0; max_isize as usize];
	for i in (0..max_isize) { //sizes.iter() {
		match h.get(&i) {
			Some(k) => counts[i as usize] = *k,
			None => {();},
		};
	}
	counts
}



	//let mut fg = Figure::new();
	//fg.axes2d()
	//  .lines(&fragsize_array.0, &fragsize_array.1, &[Caption("Proper pairs"), Color("red")]);
	//fg.set_terminal("png", "out.png");
	//fg.show();