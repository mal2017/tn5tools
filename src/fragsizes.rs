use rust_htslib::bam;
use rust_htslib::prelude::*;
use std::str;
use std::collections::HashMap;
use ndarray::{Array};
use rayon::ThreadPoolBuilder;
use rayon::prelude::*;
use csv;

/// Given a list of bams write return csv with fragment sizes
pub fn write_fragsizes(ibs: &Vec<&str>, p: usize, of: &str) {
	ThreadPoolBuilder::new().num_threads(p).build_global().unwrap();

	// if only checking one bam, send all threads to i/o
	let mut sizes: Vec<Vec<i32>> = if ibs.len() > 1 {
		ibs.par_iter().map(|a| get_fragsizes(a, 0)).collect()
	} else {
		ibs.iter().map(|a| get_fragsizes(a, p-1)).collect()
	};

	// get max size observed (len of longest sizes vec)
	let max = sizes.par_iter()
				   .map(|a| a.iter().count())
				   .max()
				   .unwrap();			
	// pad all samples to max size observed		   
	for i in &mut sizes {
		while i.len() < max {
			i.push(0);
		}
	}
	let sz_flat: Vec<i32> = sizes.into_iter().flat_map(|a| a.into_iter()).collect();
	let arr = Array::from_shape_vec((ibs.len(),max), sz_flat)
									.unwrap()
									.reversed_axes();
<<<<<<< HEAD
    let mut wtr = csv::Writer::from_path(of).unwrap();
    wtr.write_record(ibs).unwrap();
    for i in 0..max {
    	wtr.serialize(arr.row(i).to_vec()).unwrap();
    }
    wtr.flush().unwrap();
=======

	let mut fg = Figure::new();
	fg.axes2d()
	  .lines( 0..max, arr.slice(s![..,1]), &[Caption("Proper pairs"), Color("red")]);
	fg.set_terminal("png", "out.png");
	fg.show();
>>>>>>> b72e031e5568bf858071efb71ee8b57a32f93a77
}


/// Takes a bam path and returns the fragment sizes 
/// for all proper pairs.
fn get_fragsizes(ib: &str, p: usize) -> Vec<i32> {
	let mut bam = bam::Reader::from_path(ib).unwrap();
	if p > 0 {bam.set_threads(p);}
	let mut bam_rec: bam::Record = bam::Record::new(); 
	let mut tracker: HashMap<i32,i32> = HashMap::with_capacity(5000);
	let mut size: i32;
	let mut ct: i32;
	while let Ok(_x) = bam.read(&mut bam_rec) {
		if !bam_rec.is_proper_pair() | bam_rec.is_reverse() { continue };
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
	for i in 0..max_isize { 
		match h.get(&i) {
			Some(k) => counts[i as usize] = *k,
			None => {();},
		};
	}
	counts
}


