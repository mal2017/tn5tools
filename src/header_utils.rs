use rust_htslib::bam;
use rust_htslib::prelude::*;
use std::collections::HashMap;
use linear_map::LinearMap;
use itertools::Itertools;


/// From a bam header stored as a HashMap storing a vector of LinearMaps, spit out a new header.
pub fn from_hashmap(hm: &HashMap<String, Vec<LinearMap<String, String>>>) ->  bam::header::Header {

	let mut new_header = bam::header::Header::new();

	let tags = vec!["HD","SQ","PG","CL","RG"];

	let bam_header_strings: Vec<Option<Vec<String>>> = tags.iter()
											  .map(|a| header_tag_to_string_vector(a, &hm))
											  .collect();

	println!("{:?}", bam_header_strings);
	new_header
}

/// From a Linear Map generates a string in the format "KEY:VALUE".
/// Used for working with bam headers stored in HashMaps storing vectors of LinearMaps,
/// as is currently an option in rust_htslib.
fn linear_map_to_header_entry_string(lm: &LinearMap<String, String>) -> String {

	let entries: Vec<String> = lm.iter()
								 .map(|a| format!("{}:{}",a.0,a.1))
								 .collect();
	entries.join("\t")
}

/// For a supplied header tag and a ref to a hashmap of header data,
/// returns a vector of string representations.
fn header_tag_to_string_vector(tag: &str, hm: &HashMap<String, Vec<LinearMap<String, String>>>) -> Option<Vec<String>> {
	
	let record = hm.get(tag);

	match record {
		Some(r) => {
			let string_records: Vec<String> = r.iter()
		  									   .map(|a| linear_map_to_header_entry_string(a))
		  									   .collect();
		  	Some(string_records)
		},
		None => None,
	}

}