extern crate bio;
//extern crate bambam;
extern crate rust_htslib;
use std::{str, fs::File, fs::metadata, convert::TryInto};
use bio::bio_types::sequence::SequenceRead;
use rust_htslib::{bam::record::CigarString, htslib::hts_close};
use std::io::Write;
use std::io::Read;
extern crate regex;
use log::debug;
//use regex::Regex;
use regex::bytes::Regex;
use itertools::Itertools;
use rust_htslib::{bam, bam::Read as BamRead};
//use bio::stats::{PHREDProb, Prob};
use std::convert::TryFrom;
use bio::{alphabets, io};
use niffler::{compression};
use rustc_hash::FxHashMap;
use std::path::Path;
use std::io::{BufRead, BufReader};

/// simple first collection of basic sequence
/// information
#[derive(Default,Debug,Clone,PartialEq)]
pub struct SeqInfo {
    pub id: String,
    pub seq_length : u64,
    // gc content
    pub seq_gc: f32,
    // number of gaps in the sequence
    pub seq_n_gaps: u64,
    // sum of all gaps
    pub seq_bp_gaps: u64,
    // phred score quality if available
    pub seq_qual: u32
}


/// The possible black and white-list
/// as a structur. Holds the path to
/// a potentially provided file for 
/// both possibilities.
pub struct SubSet<'a> {
    pub black : Option<&'a str>,
    pub white : Option<&'a str>
}


/// This represents the formats of input and 
/// output which is currently accepted.
/// For the output this only applies for the conversion
/// and subsetting functions, for the stats this
/// obviously wont apply
#[derive(Debug)]
pub enum InOutFormat {
    Fasta,
    Fastq,
    Ubam
}


/// Describes the different possible
/// filtering options which can exits
pub enum FilterWay {
    /// all IDs in a black list will be ignored
    Black,
    /// only IDs in a black list will be kept
    White,
    /// no subsetting chosen
    None
}



#[derive(Debug,Clone,Copy)]
pub struct VersionInfo <'a>{
    /// the used program/sub-program
    pub program  : &'a str,
    /// the version of the program
    pub version  : &'a str,
    /// the author
    pub author : &'a str,
    /// the executed command
    pub command : &'a str,
}

/// This ones returns a u64 median of 
/// a vector of u64 sequences.
/// Important note: in case of a even number the 
/// float floors
/// ```rust
/// use delta_biomed_seq::median_u64;
/// let mut test1 : Vec<u64> = vec![2,2,5,5];
/// let mut test2 : Vec<u64> = vec![2,2,5,5,5];
/// assert_eq!(median_u64(&mut test1),3);
/// assert_eq!(median_u64(&mut test2),5);
/// ```
pub fn median_u64 (
    input: &mut Vec<u64>,
) -> u64 {
    input.sort();
    let len = input.len();
    if len%2 != 0 {
        input[(len/2)]
    } else {
        (input[(len/2)] + input[(len/2)-1] )/ 2
    }
}


/// This ones returns a f32 median of 
/// a vector of f32 sequences.
/// Important note: in case of a even number the 
/// float floors
/// ```rust
/// use delta_biomed_seq::median_f32;
/// let mut test1 : Vec<f32> = vec![2.0,2.0,5.0,5.0];
/// let mut test2 : Vec<f32> = vec![2.0,2.0,5.0,5.0,5.0];
/// assert_eq!(median_f32(&mut test1),3.5);
/// assert_eq!(median_f32(&mut test2),5.0);
/// ```
pub fn median_f32 (
    input: &mut Vec<f32>,
) -> f32 {
    input.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let len = input.len();
    if len%2 != 0 {
        input[(len/2)]
    } else {
        (input[(len/2)] + input[(len/2)-1] )/ 2.0
    }
}


/// adapted from here https://users.rust-lang.org/t/efficient-way-of-checking-if-two-files-have-the-same-content/74735
/// very useful for tests with external files and to verify that the results is identical
/// to a previously manually generated result file
pub fn is_same_file(file1: &str, file2: &str) -> Result<bool, std::io::Error> {
    let f1 = File::open(file1)?;
    let f2 = File::open(file2)?;

    // Check if file sizes are different
    if f1.metadata().unwrap().len() != f2.metadata().unwrap().len() {
        return Ok(false);
    }
    // Use buf readers since they are much faster
    let f1 = BufReader::new(f1);
    let f2 = BufReader::new(f2);

    // Do a byte to byte comparison of the two files
    for (b1, b2) in f1.bytes().zip(f2.bytes()) {
        if b1.unwrap() != b2.unwrap() {
            return Ok(false);
        }
    }

    Ok(true)
}

/// This ones returns a u64 Lx of 
/// a vector of u64 sequences.
/// Lx is the numer of needed sequences with sorted sequences in decreasing order where 
/// ther cumsum reaches x% of the entire sum.
/// The size can be pre-defined as well to allow calculation of LG(x) values
/// ```rust
/// use delta_biomed_seq::lx_u64;
/// let mut test1 : Vec<u64> = vec![2,3,4,5,6,7,8,9,10];
/// let mut test2 : Vec<u64> = vec![1,2,2,4,4,4,5,5,5];
/// assert_eq!(lx_u64(&mut test1,0.5,None),3);
/// assert_eq!(lx_u64(&mut test2,0.5,None),4);
/// assert_eq!(lx_u64(&mut test1,0.5,Some(54)),3);
/// assert_eq!(lx_u64(&mut test1,0.5,Some(60)),4);
/// assert_eq!(lx_u64(&mut test1,0.0,None),1);
/// assert_eq!(lx_u64(&mut test1,1.0,None),9);
/// ```
pub fn lx_u64 (
    input: &mut Vec<u64>,
    fract: f64,
    g_size: Option<u64>
) -> u64 {
    input.sort();
    let total_sum = match g_size {
        Some(x) => x,
        None => cumsum_u64(input),
    };
    let mut sum: u64 = 0;
    let mut counter = input.len();
    while sum < (total_sum as f64 * fract) as u64 {
        if counter == 0 { 
            debug!("Counter is 0, breaking the loop");
            break;
         };
        counter -= 1;
        //eprintln!("counter: {} element {}", counter, input[counter]);
        sum += input[counter];
    };
    // in case we have a fraction of 0 or similar events
    // the counter will exceed by one, so we to correct this
    if counter == input.len() { counter -=1 };
    input.len().abs_diff(counter) as u64

}


/// This ones returns a u64 Nx of 
/// a vector of u64 sequences.
/// Nx is the length of the sequence from sorted sequences in decreasing order where 
/// ther cumsum reaches x% of the entire sum
/// The size can be pre-defined as well to allow calculation of NG(x) values
/// ```rust
/// use delta_biomed_seq::nx_u64;
/// let mut test1 : Vec<u64> = vec![2,3,4,5,6,7,8,9,10];
/// let mut test2 : Vec<u64> = vec![1,2,2,4,4,4,5,5,5];
/// assert_eq!(nx_u64(&mut test1,0.5,None),8);
/// assert_eq!(nx_u64(&mut test2,0.5,None),4);
/// assert_eq!(nx_u64(&mut test1,0.5,Some(54)),8);
/// assert_eq!(nx_u64(&mut test1,0.5,Some(60)),7);
/// assert_eq!(nx_u64(&mut test1,0.0,None),10);
/// assert_eq!(nx_u64(&mut test1,1.0,None),2);
/// ```
pub fn nx_u64 (
    input: &mut Vec<u64>,
    fract: f64,
    g_size: Option<u64>
) -> u64 {
    input.sort();
    let total_sum = match g_size {
        Some(x) => x,
        None => cumsum_u64(input),
    };
    let mut sum: u64 = 0;
    let mut counter = input.len();
    //debug!("counter: {} sum {} total_sum {} fract {} g_size {:?}", counter, &sum, &total_sum,&fract,&g_size);
    while sum < (total_sum as f64 * fract) as u64 {
        if counter == 0 { 
            debug!("Counter is 0, breaking the loop");
            break;
         };
        counter -= 1;
        //debug!("Counter: {}", counter);
        sum += input[counter];
        //debug!("counter: {} sum {} total_sum {} fract {} g_size {:?}", counter, &sum, &total_sum,&fract,&g_size);
    };
    // in case we have a fraction of 0 or similar events
    // the counter will exceed by one, so we to correct this
    if counter == input.len() { counter -=1 };
    debug!("counter: {} sum {} total_sum {} fract {} g_size {:?} result {}", counter, &sum, &total_sum,&fract,&g_size, input[counter]);
    input[counter]

}


/// This ones returns a cumsum of a u64 vector
/// ```rust
/// use delta_biomed_seq::cumsum_u64;
/// let test1 : Vec<u64> = vec![2,2,5,5];
/// let test2 : Vec<u64> = vec![2,2,5,5,5];
/// assert_eq!(cumsum_u64(&test1),14);
/// assert_eq!(cumsum_u64(&test2),19);
/// ```
pub fn cumsum_u64 (
    input: &Vec<u64>,
) -> u64 {
    let mut out : u64 = 0;
    for entry in input {
        out += entry;
    }
    out
}


/// This ones returns a mean of a u64 vector.
/// Note that the returned value is floored
/// ```rust
/// use delta_biomed_seq::mean_u64;
/// let test1 : Vec<u64> = vec![2,2,5,5,5,5];
/// let test2 : Vec<u64> = vec![2,2,5,5,5];
/// assert_eq!(mean_u64(&test1),4);
/// assert_eq!(mean_u64(&test2),3);
/// ```
pub fn mean_u64 (
    input: &Vec<u64>,
) -> u64 {
    let sum = cumsum_u64(input);
    let len = input.len() as u64;
    sum/len
}



pub fn get_ids(
    my_file: &str
) -> FxHashMap<String, i8> {
    // this function takes a file and returns a hashmap
    // of all entries in the form entry=>1
    // this is used later to serve as list of entries
    // which are queried

    let mut tmp_map : FxHashMap<String, i8> = FxHashMap::default();
    debug!("INFO: ID list {} provided, reading entries...", my_file);
    assert!(
        Path::new(my_file).exists(),
        "ERROR: ID list file {:?} does not exist!",
        my_file
    );

    let input  = File::open(&my_file).expect("Unable to open filter file");
    let reader = BufReader::new(input);
        
    for line in reader.lines() {
        let l = line.expect("ERROR: could not read line!");
        // now we split by empty space
        let e: Vec<&str> = l.split('\t').collect();
        if e.len() != 1 {
            panic!("ERROR: your ID list contains space delimited entries!");
        }
        tmp_map.insert(e[0].to_string(), 1);
    }
    tmp_map
}



/// converts fasta to fasta. This is useful to guarantee that
/// the final output is correctly formatted. Provides as well
/// functionality to only look at a subset of sequences
pub fn fa2fa (
    input: &str, 
    output: &str , 
    filter: SubSet,
    force_upper: bool,
    force_simple: bool,
){
    let (unflate_stream, _compression_format) = niffler::from_path(input).expect("ERROR: could not open path for reading!");
    let reader     = io::fasta::Reader::new(unflate_stream);
    
    let mut writer = match input.ends_with(".gz") {
        true  => io::fasta::Writer::new(niffler::to_path(&output, compression::Format::Gzip, niffler::Level::Six).expect("ERROR: writer issue!")),
        false => io::fasta::Writer::new(niffler::to_path(&output, compression::Format::No, niffler::Level::Six).expect("ERROR: writer issue!")),
    };

    //debug!("Reader {:?}",reader);
    // first verify if we have a black or a white list potentially
    let filtering  = match filter {
        SubSet {white: Some(_x), black: None} => FilterWay::White,
        SubSet {white: None, black: Some(_x)} => FilterWay::Black,
        SubSet {white: None, black: None}     => FilterWay::None,
        SubSet {white: Some(_x), black: Some(_y)}     => panic!("ERROR: black and white list define!"),
    };
    // if we have one or the other than we need to get the IDs 
    // next
    let filter_ids = match filtering {
        FilterWay::Black => get_ids(filter.black.unwrap()),
        FilterWay::White => get_ids(filter.white.unwrap()),
        FilterWay::None  => FxHashMap::default(),
    };
    
    for record in reader.records() {
        let mut record = record.expect("ERROR: could not read fa entry!");

        // if no filter, write directly
        // if blacklist and ID in list then skip
        // if whitelist and ID in list then keep
        match filtering {
            FilterWay::Black => {
                if !filter_ids.contains_key(record.id()){
                    match check_fa(&record) {
                        Ok(_) => (),
                        Err(r) => println!("ERROR: {:?} --> QC NOT PASSED! {:?} ",record,r),
                    };
                    
                    // this is NOT sufficient as they do only
                    // check integrity of header and if ascii for sequence
                    // no check if part of nuc oder AA alphabet!
                    // see https://github.com/rust-bio/rust-bio/issues/472
                    if force_simple {
                        record = simple_iupac(&record);
                    }; 
                    if force_upper {
                        eprintln!("INFO: Forcing conversion from non ATCGN IUPAC compatible nucleotides to Ns");
                        let tmp_seq = std::str::from_utf8(record.seq() ).expect("ERROR: could not convert sequence !").to_uppercase();
                        record = bio::io::fasta::Record::with_attrs(record.id(), record.desc(), tmp_seq.as_bytes());
            
                    }

                    writer.write_record(&record).expect("ERROR: could not write out record!")
                }
            },
            FilterWay::White => {
                if filter_ids.contains_key(record.id()){
                    match check_fa(&record) {
                        Ok(_) => (),
                        Err(r) => println!("ERROR: {:?} --> QC NOT PASSED! {:?} ",record,r),
                    };
                    // this is NOT sufficient as they do only
                    // check integrity of header and if ascii for sequence
                    // no check if part of nuc oder AA alphabet!
                    // see https://github.com/rust-bio/rust-bio/issues/472
                    if force_simple {
                        record = simple_iupac(&record);
                    };
                    if force_upper {
                        let tmp_seq = std::str::from_utf8(record.seq() ).expect("ERROR: could not convert sequence !").to_uppercase();
                        record = bio::io::fasta::Record::with_attrs(record.id(), record.desc(), tmp_seq.as_bytes());
            
                    }

                    writer.write_record(&record).expect("ERROR: could not write out record!")
                }
            },
            FilterWay::None => {
                match check_fa(&record) {
                    Ok(_) => (),
                    Err(r) => println!("ERROR: {:?} --> QC NOT PASSED! {:?} ",record,r),
                };
                // this is NOT sufficient as they do only
                // check integrity of header and if ascii for sequence
                // no check if part of nuc oder AA alphabet!
                // see https://github.com/rust-bio/rust-bio/issues/472
                if force_simple {
                    record = simple_iupac(&record);
                };
                if force_upper {
                    let tmp_seq = std::str::from_utf8(record.seq() ).expect("ERROR: could not convert sequence !").to_uppercase();
                    record = bio::io::fasta::Record::with_attrs(record.id(), record.desc(), tmp_seq.as_bytes());
        
                }

                writer.write_record(&record).expect("ERROR: could not write out record!")
            }
        }
        
    }
    writer.flush().expect("ERROR: could not close output file!");
}



/// converts fastq to fastq. Useful guarantee that files are
/// according to format and subsetting allows to keep/remove entries
/// of interest
/// gathers stats for uBAM
/// sequences
/// ```rust
/// use delta_biomed_seq::SeqInfo;
/// use std::fs::File;
/// use std::io::{self, Write};
/// use tempfile::tempdir;
/// use std::convert::TryFrom;
/// use std::fs;
/// use delta_biomed_seq::SubSet;
/// use delta_biomed_seq::fq2fq;
/// use niffler::{compression};
/// use delta_biomed_seq::is_same_file;
/// 
/// let tmp_dir = tempdir().unwrap();
/// let tmp_path = tmp_dir.into_path();
/// let file_path = &tmp_path.join("tmp.fastq");
/// let file_path_out = &tmp_path.join("tmp2.fastq");
/// let file_path_out2 = &tmp_path.join("tmp3.fastq");
/// let mut tmp_file = File::create(file_path).unwrap();
/// /// here, all of them are phred score 20 
/// let record = bio::io::fastq::Record::with_attrs("id_str", Some("desc"), b"TTTTTTTTTTCCCCCCCCCCAAAAAAAAAAGGGGGGGGGG",b"5555555555555555555555555555555555555555");
/// eprintln!("{:?}",record.qual());
/// let mut writer = bio::io::fastq::Writer::new(niffler::to_path(&file_path, compression::Format::No, niffler::Level::Nine).expect("ERROR: writer issue!"));
/// writer.write_record(&record).expect("ERROR: could not write out record!");
/// drop(writer);
/// // Here we calculate now the quality first from rq field
/// fq2fq(file_path.to_str().unwrap(),file_path_out.to_str().unwrap(),SubSet{ black : None, white : None},Some(0.9));
/// let comparison = is_same_file(file_path.to_str().unwrap(),file_path_out.to_str().unwrap());
/// assert_eq!(comparison.unwrap(), true);
/// fq2fq(file_path.to_str().unwrap(),file_path_out2.to_str().unwrap(),SubSet{ black : None, white : None},Some(0.999));
/// let comparison2 = is_same_file(file_path.to_str().unwrap(),file_path_out2.to_str().unwrap());
/// assert_eq!(comparison2.unwrap(), false);
/// ```
pub fn fq2fq (
    input: &str, 
    output: &str , 
    filter: SubSet,
    rq: Option<f32>
){
    let rq_val : Option <u32> = match rq{
        Some(x) => {
             let er = rq.unwrap();
            let phred : f32;
            if er < 0.0 {
                phred = 10.0 ;
            }else if er > 0.999999 {
                phred = 60.0 ;
            }else{
                phred = -10_f32 * f32::log10(1_f32-er);
                //eprintln!("{}",&phred);
            }
            Some(phred.round() as u32)
        },
        None => None
    };
    let (unflate_stream, _compression_format) = niffler::from_path(input).expect("ERROR: could not open path for reading!");
    let reader     = io::fastq::Reader::new(unflate_stream);
    let mut writer = match output.ends_with(".gz") {
        true  => io::fastq::Writer::new(niffler::to_path(&output, compression::Format::Gzip, niffler::Level::Nine).expect("ERROR: writer issue!")),
        false => io::fastq::Writer::new(niffler::to_path(&output, compression::Format::No, niffler::Level::Nine).expect("ERROR: writer issue!")),
    };

    // first verify if we have a black or a white list potentially
    let filtering  = match filter {
        SubSet {white: Some(_x), black: None} => FilterWay::White,
        SubSet {white: None, black: Some(_x)} => FilterWay::Black,
        SubSet {white: None, black: None}     => FilterWay::None,
        SubSet {white: Some(_x), black: Some(_y)}     => panic!("ERROR: black and white list define!"),
    };
    // if we have one or the other than we need to get the IDs 
    // next
    let filter_ids = match filtering {
        FilterWay::Black => get_ids(filter.black.unwrap()),
        FilterWay::White => get_ids(filter.white.unwrap()),
        FilterWay::None  => FxHashMap::default(),
    };
    for record in reader.records() {
        debug!("Record: {:?}",&record);
        let record = record.unwrap();
        if rq_val.is_some() {
            // pacbio is actually cheating and caps the rq 
            // at 1 which is not correct. As this returns otherwise just
            // the max of u32 I will instead cap to phred 60 for the moment,
            // which is  99.9999% quality
            let mut quals : u32 = 0;
            //let mut quals = trans_record.qual().iter().sum::<u8>() as f64;
            //eprintln!("qual: {:?}",trans_record.qual().to_vec());
            for entry in record.qual(){
                let mut tmp: u32  = entry.to_string().parse().unwrap();
                tmp = tmp-33;
                quals += tmp;
                debug!("entry: {}, sum: {}",&tmp,&quals);
            }
            
            let tmp_result =  quals as u64/record.len() as u64 ;
            debug!("INFO: quality of read: {:?}, quality cut-off {:?}",tmp_result,rq_val);
            if tmp_result < rq_val.unwrap().into(){
                continue
            }
        }
        // if no filter, write directly
        // if blacklist and ID in list then skip
        // if whitelist and ID in list then keep
        match filtering {
            FilterWay::Black => {
                if !filter_ids.contains_key(record.id()){
                // this is NOT sufficient as they do only
                // check integrity of header and if ascii for sequence
                // no check if part of nuc oder AA alphabet!
                // see https://github.com/rust-bio/rust-bio/issues/472
                match check_fq(&record) {
                    Ok(_) => (),
                    Err(r) => println!("ERROR: {:?} --> QC NOT PASSED! {:?} ",record,r),
                };
                    writer.write_record(&record).expect("ERROR: could not write out record!")
                }
            },
            FilterWay::White => {
                if filter_ids.contains_key(record.id()){
                    // this is NOT sufficient as they do only
                    // check integrity of header and if ascii for sequence
                    // no check if part of nuc oder AA alphabet!
                    // see https://github.com/rust-bio/rust-bio/issues/472
                    match check_fq(&record) {
                        Ok(_) => (),
                        Err(r) => println!("ERROR: {:?} --> QC NOT PASSED! {:?} ",record,r),
                    };
                    writer.write_record(&record).expect("ERROR: could not write out record!")
                }
            },
            FilterWay::None => {
                // this is NOT sufficient as they do only
                // check integrity of header and if ascii for sequence
                // no check if part of nuc oder AA alphabet!
                // see https://github.com/rust-bio/rust-bio/issues/472
                match check_fq(&record) {
                    Ok(_) => (),
                    Err(r) => println!("ERROR: {:?} --> QC NOT PASSED! {:?} ",record,r),
                };
                writer.write_record(&record).expect("ERROR: could not write out record!")
            }
        }
    }
    writer.flush().expect("ERROR: could not close output file!");
}


/// converts un-mapped BAM (uBAM) to uBAM. Useful in order guarantee that no
/// format problem might be encountered. 
/// Furthermore allows subsetting through black and white-lists and based on
/// quality in AUX field rq.
pub fn ubam2ubam (
    input: &str,
    output: &str,
    filter: SubSet,
    threads: u32 ,
    args: &str,
    id:  &str,
    vn:  &str,
    rq: Option<f32>,
    trash: Option<&str>
){
    let mut reader = bam::Reader::from_path(input).expect("ERROR: could not open input file");
    let mut header = bam::Header::from_template(reader.header());
    // generate a thread pool for the multithreaded reading and writing of BAM files
    let pool = rust_htslib::tpool::ThreadPool::new(threads).unwrap();
    reader.set_thread_pool(&pool).unwrap();
    
    header.push_record(
        bam::header::HeaderRecord::new(b"PG")
            .push_tag(b"ID", &id)
            .push_tag(b"PN", &id)
            .push_tag(b"VN", &vn)
            .push_tag(b"CL",&args)
    );
    // first verify if we have a black or a white list potentially
    let filtering  = match filter {
        SubSet {white: Some(_x), black: None} => FilterWay::White,
        SubSet {white: None, black: Some(_x)} => FilterWay::Black,
        SubSet {white: None, black: None}     => FilterWay::None,
        SubSet {white: Some(_x), black: Some(_y)}     => panic!("ERROR: black and white list define!"),
    };
    // if we have one or the other than we need to get the IDs 
    // next
    let filter_ids = match filtering {
        FilterWay::Black => get_ids(filter.black.unwrap()),
        FilterWay::White => get_ids(filter.white.unwrap()),
        FilterWay::None  => FxHashMap::default(),
    };

    let mut writer = bam::Writer::from_path(output, &header, bam::Format::Bam ).expect("ERROR: could not open output file");
    // currently to simplfiy this I make a dummy writer 
    // hope that this is a good idea
    let mut trash_writer = match trash {
        Some(x) => bam::Writer::from_path(x, &header, bam::Format::Bam ).expect("ERROR: could not open non-selected output file"),
        None         => bam::Writer::from_path("/tmp/tmp.bam", &header, bam::Format::Bam ).expect("ERROR: could not open dummy path /dev/null"),
    };
    writer.set_thread_pool(&pool).unwrap();
    for r in reader.records() {
        let record = r.expect("ERROR: could not read BAM entry!");
        // here we filter now already the record if it is below a asked quality
        // This is predominantly for HIFI settings
        if rq.is_some() && match record.aux(b"rq").unwrap() {
                rust_htslib::bam::record::Aux::Float(x) => x,
                _ => panic!("ERROR: no float in RQ field encountered!"),
            } < rq.unwrap(){
            if trash.is_some() {
                //println!("writing to trash");
                trash_writer.write(&record).expect("ERROR: could not write trash out record!");
            }
            continue
        }
        
        // if no filter, write directly
        // if blacklist and ID in list then skip
        // if whitelist and ID in list then keep
        match filtering {
            FilterWay::Black => {
                if !filter_ids.contains_key(str::from_utf8(record.qname()).unwrap()){
                    writer.write(&record).expect("ERROR: could not write out record!");
                }else if trash.is_some() {
                    trash_writer.write(&record).expect("ERROR: could not write trash out record!");
                }
            },
            FilterWay::White => {
                if filter_ids.contains_key(str::from_utf8(record.qname()).unwrap()){
                    writer.write(&record).expect("ERROR: could not write out record!");
                }else if trash.is_some() {
                    trash_writer.write(&record).expect("ERROR: could not write trash out record!");
                }
            },
            FilterWay::None => {
                writer.write(&record).expect("ERROR: could not write out record!");
            }
        }
         
    }
}



/// Verifies that output is according to format and provides capability to subset 
/// the sequences
pub fn fa2fq (
    input: &str, 
    output: &str , 
    qual: &str, 
    filter: SubSet
) {
    let (unflate_stream, _compression_format) = niffler::from_path(input).expect("ERROR: could not open path for reading!");
    let reader     = io::fasta::Reader::new(unflate_stream);

    let mut writer = match input.ends_with(".gz") {
        true  => io::fastq::Writer::new(niffler::to_path(&output, compression::Format::Gzip, niffler::Level::Nine).expect("ERROR: writer issue!")),
        false => io::fastq::Writer::new(niffler::to_path(&output, compression::Format::No, niffler::Level::Nine).expect("ERROR: writer issue!")),
    };
    // first verify if we have a black or a white list potentially
    let filtering  = match filter {
        SubSet {white: Some(_x), black: None} => FilterWay::White,
        SubSet {white: None, black: Some(_x)} => FilterWay::Black,
        SubSet {white: None, black: None}     => FilterWay::None,
        SubSet {white: Some(_x), black: Some(_y)}     => panic!("ERROR: black and white list define!"),
    };
    // if we have one or the other than we need to get the IDs 
    // next
    let filter_ids = match filtering {
        FilterWay::Black => get_ids(filter.black.unwrap()),
        FilterWay::White => get_ids(filter.white.unwrap()),
        FilterWay::None  => FxHashMap::default(),
    };
    for record in reader.records() {
        let record = record.unwrap();

        // if no filter, write directly
        // if blacklist and ID in list then skip
        // if whitelist and ID in list then keep
        match filtering {
            FilterWay::Black => {
                if !filter_ids.contains_key(record.id()){
                    match check_fa(&record) {
                        Ok(_) => (),
                        Err(r) => println!("ERROR: {:?} --> QC NOT PASSED! {:?} ",record,r),
                    };
                    writer.write_record(&fa2fq_record(record,qual)).expect("ERROR: could not write out record!");
                }
            },
            FilterWay::White => {
                if filter_ids.contains_key(record.id()){
                    match check_fa(&record) {
                        Ok(_) => (),
                        Err(r) => println!("ERROR: {:?} --> QC NOT PASSED! {:?} ",record,r),
                    };
                    writer.write_record(&fa2fq_record(record,qual)).expect("ERROR: could not write out record!");
                }
            },
            FilterWay::None => {
                match check_fa(&record) {
                    Ok(_) => (),
                    Err(r) => println!("ERROR: {:?} --> QC NOT PASSED! {:?} ",record,r),
                };
                writer.write_record(&fa2fq_record(record,qual)).expect("ERROR: could not write out record!");
            }
        }
        
    }
    writer.flush().expect("ERROR: could not close output file!");
}




/// converts fastq to fasta and removes in the process all of the
/// quality values, Verifies format and provides capability to subset the sequences.
pub fn fq2fa (
    input: &str, 
    output: &str , 
    filter: SubSet
){
    let (unflate_stream, _compression_format) = niffler::from_path(input).expect("ERROR: could not open path for reading!");
    let reader     = io::fastq::Reader::new(unflate_stream);
    //debug!("reader");
    let mut writer = match input.ends_with(".gz") {
        true  => io::fasta::Writer::new(niffler::to_path(&output, compression::Format::Gzip, niffler::Level::Nine).expect("ERROR: writer issue!")),
        false => io::fasta::Writer::new(niffler::to_path(&output, compression::Format::No, niffler::Level::Nine).expect("ERROR: writer issue!")),
    };
   
    // first verify if we have a black or a white list potentially
    let filtering  = match filter {
        SubSet {white: Some(_x), black: None} => FilterWay::White,
        SubSet {white: None, black: Some(_x)} => FilterWay::Black,
        SubSet {white: None, black: None}     => FilterWay::None,
        SubSet {white: Some(_x), black: Some(_y)}     => panic!("ERROR: black and white list define!"),
    };
    // if we have one or the other than we need to get the IDs 
    // next
    let filter_ids = match filtering {
        FilterWay::Black => get_ids(filter.black.unwrap()),
        FilterWay::White => get_ids(filter.white.unwrap()),
        FilterWay::None  => FxHashMap::default(),
    };
    for record in reader.records() {
        let record = record.expect("ERROR: could not read fq entry!");
        
        // if no filter, write directly
        // if blacklist and ID in list then skip
        // if whitelist and ID in list then keep
        match filtering {
            FilterWay::Black => {
                if !filter_ids.contains_key(record.id()){
                    match check_fq(&record) {
                        Ok(_) => (),
                        Err(r) => println!("ERROR: {:?} --> QC NOT PASSED! {:?} ",record,r),
                    };
                    writer.write_record(&fq2fa_record(record)).expect("ERROR: could not write out record!");
                }
            },
            FilterWay::White => {
                if filter_ids.contains_key(record.id()){
                    match check_fq(&record) {
                        Ok(_) => (),
                        Err(r) => println!("ERROR: {:?} --> QC NOT PASSED! {:?} ",record,r),
                    };
                    writer.write_record(&fq2fa_record(record)).expect("ERROR: could not write out record!");
                }
            },
            FilterWay::None => {
                match check_fq(&record) {
                    Ok(_) => (),
                    Err(r) => println!("ERROR: {:?} --> QC NOT PASSED! {:?} ",record,r),
                };
                writer.write_record(&fq2fa_record(record)).expect("ERROR: could not write out record!");
            }
        }
        
    }
    writer.flush().expect("ERROR: could not close output file!");
}



/// converts uBAM to fastq which allows to keep id, sequences and quality but will
/// remove all other meta-data (e.g. PacBio kinetics annotation).
/// Checks that output format is valid and provides capability to subset sequences.
pub fn ubam2fq (
    input: &str,
    output: &str ,
    filter: SubSet,
    threads: u32,
    rq:Option<f32>,
) {
    let mut reader = bam::Reader::from_path(input).expect("ERROR: could not open input file");

    let mut writer = match input.ends_with(".gz") {
        true  => io::fastq::Writer::new(niffler::to_path(&output, compression::Format::Gzip, niffler::Level::Nine).expect("ERROR: writer issue!")),
        false => io::fastq::Writer::new(niffler::to_path(&output, compression::Format::No, niffler::Level::Nine).expect("ERROR: writer issue!")),
    };
    // generate a thread pool for the multithreaded reading and writing of BAM files
    let pool = rust_htslib::tpool::ThreadPool::new(threads).unwrap();
    reader.set_thread_pool(&pool).unwrap();

    // first verify if we have a black or a white list potentially
    let filtering  = match filter {
        SubSet {white: Some(_x), black: None} => FilterWay::White,
        SubSet {white: None, black: Some(_x)} => FilterWay::Black,
        SubSet {white: None, black: None}          => FilterWay::None,
        SubSet {white: Some(_x), black: Some(_y)}     => panic!("ERROR: black and white list define!"),
    };
    // if we have one or the other than we need to get the IDs 
    // next
    let filter_ids = match filtering {
        FilterWay::Black => get_ids(filter.black.unwrap()),
        FilterWay::White => get_ids(filter.white.unwrap()),
        FilterWay::None  => FxHashMap::default(),
    };
    for record in reader.records() {
        let record = record.expect("ERROR: could not read BAM entry!");
        if rq.is_some() && match record.aux(b"rq").unwrap() {
            rust_htslib::bam::record::Aux::Float(x) => x,
            _ => panic!("ERROR: no float in RQ field encountered!"),
        } < rq.unwrap(){
            continue
        }
        // if no filter, write directly
        // if blacklist and ID in list then skip
        // if whitelist and ID in list then keep
        match filtering {
            FilterWay::Black => {
                if !filter_ids.contains_key(str::from_utf8(record.qname()).unwrap()){
                    writer.write_record(&ubam2fq_record(record)).expect("ERROR: could not write out record!");
                }
            },
            FilterWay::White => {
                if filter_ids.contains_key(str::from_utf8(record.qname()).unwrap()){
                    writer.write_record(&ubam2fq_record(record)).expect("ERROR: could not write out record!");
                }
            },
            FilterWay::None => {
                let new_record = ubam2fq_record(record);
                writer.write_record(&new_record).expect("ERROR: could not write out record!");
            }
        }
        
    }
    writer.flush().expect("ERROR: could not close output file!");
}



/// converts uBAM to fasta which allows to keep id, sequences and will
/// remove all other meta-data (e.g. PacBio kinetics annotation).
/// Checks that output format is valid and provides capability to subset sequences.
pub fn ubam2fa (
    input: &str,
    output: &str ,
    filter: SubSet,
    threads: u32,
    rq: Option<f32>
){
    let mut reader = bam::Reader::from_path(input).expect("ERROR: could not open input file");

    let mut writer = match input.ends_with(".gz") {
        true  => io::fasta::Writer::new(niffler::to_path(&output, compression::Format::Gzip, niffler::Level::Nine).expect("ERROR: writer issue!")),
        false => io::fasta::Writer::new(niffler::to_path(&output, compression::Format::No, niffler::Level::Nine).expect("ERROR: writer issue!")),
    };    // generate a thread pool for the multithreaded reading and writing of BAM files
    let pool = rust_htslib::tpool::ThreadPool::new(threads).unwrap();
    reader.set_thread_pool(&pool).unwrap();
    // first verify if we have a black or a white list potentially
    let filtering  = match filter {
        SubSet {white: Some(_x), black: None} => FilterWay::White,
        SubSet {white: None, black: Some(_x)} => FilterWay::Black,
        SubSet {white: None, black: None}     => FilterWay::None,
        SubSet {white: Some(_x), black: Some(_y)}     => panic!("ERROR: black and white list define!"),
    };
    // if we have one or the other than we need to get the IDs 
    // next
    let filter_ids = match filtering {
        FilterWay::Black => get_ids(filter.black.unwrap()),
        FilterWay::White => get_ids(filter.white.unwrap()),
        FilterWay::None  => FxHashMap::default(),
    };
    for record in reader.records() {
        let record = record.expect("ERROR: could not read BAM entry!");
        if rq.is_some() && match record.aux(b"rq").unwrap() {
            rust_htslib::bam::record::Aux::Float(x) => x,
            _ => panic!("ERROR: no float in RQ field encountered!"),
        } < rq.unwrap(){
            continue
        }
        // if no filter, write directly
        // if blacklist and ID in list then skip
        // if whitelist and ID in list then keep
        match filtering {
            FilterWay::Black => {
                if !filter_ids.contains_key(str::from_utf8(record.qname()).unwrap()){
                    writer.write_record(&ubam2fa_record(record)).expect("ERROR: could not write out record!");
                }
            },
            FilterWay::White => {
                if filter_ids.contains_key(str::from_utf8(record.qname()).unwrap()){
                    writer.write_record(&ubam2fa_record(record)).expect("ERROR: could not write out record!");
                }
            },
            FilterWay::None => {
                writer.write_record(&ubam2fa_record(record)).expect("ERROR: could not write out record!");
            }
        }
        
    }
    writer.flush().expect("ERROR: could not close output file!");
}




/// converts fastq to uBAM which allows to keep id, sequences and quality.
/// All other aspects will be set to the default and will generate unmapped entries.
/// Checks that output format is valid and provides capability to subset sequences.
/// 
/// /// Example of PacBio uBAM entries:
/// m64030_200808_213327/1/0 4735   4  *  0  255  *  *  0  0  CT...ATGA !!!!..!!! 
pub fn fq2ubam (
    input: &str, 
    output: &str , 
    rg: &str, 
    filter: SubSet, 
    args: &str, 
    id:  &str, 
    vn:  &str
){
    let (unflate_stream, _compression_format) = niffler::from_path(input).expect("ERROR: could not open path for reading!");
    let reader     = io::fastq::Reader::new(unflate_stream);
    
    let mut header = bam::Header::default();
    header.push_record(
        bam::header::HeaderRecord::new(b"RG")
            .push_tag(b"ID", &"1")
            .push_tag(b"SM", &"unknown")
            .push_tag(b"PL", &"unknown")
    );
    header.push_record(
        bam::header::HeaderRecord::new(b"PG")
            .push_tag(b"ID", &id)
            .push_tag(b"PN", &id)
            .push_tag(b"VN", &vn)
            .push_tag(b"CL",&args)
    );
    let mut writer = bam::Writer::from_path(output, &header, bam::Format::Bam).expect("ERROR: could not open output file");
    // first verify if we have a black or a white list potentially
    let filtering  = match filter {
        SubSet {white: Some(_x), black: None} => FilterWay::White,
        SubSet {white: None, black: Some(_x)} => FilterWay::Black,
        SubSet {white: None, black: None}     => FilterWay::None,
        SubSet {white: Some(_x), black: Some(_y)}     => panic!("ERROR: black and white list define!"),
    };
    // if we have one or the other than we need to get the IDs 
    // next
    let filter_ids = match filtering {
        FilterWay::Black => get_ids(filter.black.unwrap()),
        FilterWay::White => get_ids(filter.white.unwrap()),
        FilterWay::None  => FxHashMap::default(),
    };
    for record in reader.records() {
        let record = record.expect("ERROR: could not read BAM entry!");
        
        // if no filter, write directly
        // if blacklist and ID in list then skip
        // if whitelist and ID in list then keep
        match filtering {
            FilterWay::Black => {
                if !filter_ids.contains_key(record.id()){
                    match check_fq(&record) {
                        Ok(_) => (),
                        Err(r) => println!("ERROR: {:?} --> QC NOT PASSED! {:?} ",record,r),
                    };
                    writer.write(&fq2ubam_record(record,rg)).expect("ERROR: could not write out record!");
                }
            },
            FilterWay::White => {
                if filter_ids.contains_key(record.id()){
                    match check_fq(&record) {
                        Ok(_) => (),
                        Err(r) => println!("ERROR: {:?} --> QC NOT PASSED! {:?} ",record,r),
                    };
                    writer.write(&fq2ubam_record(record,rg)).expect("ERROR: could not write out record!");
                }
            },
            FilterWay::None => {
                match check_fq(&record) {
                    Ok(_) => (),
                    Err(r) => println!("ERROR: {:?} --> QC NOT PASSED! {:?} ",record,r),
                };
                writer.write(&fq2ubam_record(record,rg)).expect("ERROR: could not write out record!");
            }
        }
        
    }
}




/// converts fastq to uBAM which allows to keep id, sequences and quality.
/// All other aspects will be set to the default and will generate unmapped entries.
/// A dummy quality value (symbol) must be provided
/// in order to fill the needed quality track.
/// Checks that output format is valid and provides capability to subset sequences.
pub fn fa2ubam (
    input: &str, 
    output: &str , 
    qual: &str, 
    rg: &str , 
    filter: SubSet, 
    args: &str, 
    id:  &str, 
    vn:  &str
) {
    let (unflate_stream, _compression_format) = niffler::from_path(input).expect("ERROR: could not open path for reading!");
    let reader     = io::fasta::Reader::new(unflate_stream);
    
    let mut header = bam::Header::default();
    header.push_record(
        bam::header::HeaderRecord::new(b"PG")
            .push_tag(b"ID", &id)
            .push_tag(b"PN", &id)
            .push_tag(b"VN", &vn)
            .push_tag(b"CL",&args)
    );
    let mut writer = bam::Writer::from_path(output, &header, bam::Format::Bam).expect("ERROR: could not open output file");
    // first verify if we have a black or a white list potentially
    let filtering  = match filter {
        SubSet {white: Some(_x), black: None} => FilterWay::White,
        SubSet {white: None, black: Some(_x)} => FilterWay::Black,
        SubSet {white: None, black: None}     => FilterWay::None,
        SubSet {white: Some(_x), black: Some(_y)}     => panic!("ERROR: black and white list define!"),
    };
    // if we have one or the other than we need to get the IDs 
    // next
    let filter_ids = match filtering {
        FilterWay::Black => get_ids(filter.black.unwrap()),
        FilterWay::White => get_ids(filter.white.unwrap()),
        FilterWay::None  => FxHashMap::default(),
    };
    for record in reader.records() {
        let record = record.expect("ERROR: could not read BAM entry!");
        
        // if no filter, write directly
        // if blacklist and ID in list then skip
        // if whitelist and ID in list then keep
        match filtering {
            FilterWay::Black => {
                if !filter_ids.contains_key(record.id()){
                    match check_fa(&record) {
                        Ok(_) => (),
                        Err(r) => println!("ERROR: {:?} --> QC NOT PASSED! {:?} ",record,r),
                    };
                    writer.write(&fa2ubam_record(record,qual,rg)).expect("ERROR: could not write out record!");
                }
            },
            FilterWay::White => {
                if filter_ids.contains_key(record.id()){
                    match check_fa(&record) {
                        Ok(_) => (),
                        Err(r) => println!("ERROR: {:?} --> QC NOT PASSED! {:?} ",record,r),
                    };
                    writer.write(&fa2ubam_record(record,qual,rg)).expect("ERROR: could not write out record!");
                }
            },
            FilterWay::None => {
                match check_fa(&record) {
                    Ok(_) => (),
                    Err(r) => println!("ERROR: {:?} --> QC NOT PASSED! {:?} ",record,r),
                };
                writer.write(&fa2ubam_record(record,qual,rg)).expect("ERROR: could not write out record!");
            }
        }
        
    }
}




/// gathers stats for fasta
/// sequnces
/// ```rust
/// use delta_biomed_seq::fa_stats;
/// use delta_biomed_seq::SeqInfo;
/// use std::convert::TryFrom;
/// use bio::io::fasta::Record;
/// use bio::io::fasta::Writer;
/// use tempfile::NamedTempFile;
/// 
///  let tmp = NamedTempFile::new().unwrap();
/// let path = tmp.path().to_str().unwrap();
/// let record = bio::io::fasta::Record::with_attrs("id_str", Some("desc"), b"TTTTTTTTTTCCCCCCCCCCAAAAAAAAAAGGGGGGGGGG");
/// let mut writer = bio::io::fasta::Writer::to_file(path).expect("ERROR: cant write to file!");
/// writer.write_record(&record).expect("ERROR: could not write out record!");
/// writer.flush().unwrap();
/// let result = fa_stats(path);
/// let mut verified : Vec<SeqInfo> = Vec::default();
/// let mut entry : SeqInfo = SeqInfo::default();
/// entry.seq_length = 40;
/// entry.seq_gc =0.5;
/// entry.id = String::from("id_str");
/// verified.push(entry);
/// assert_eq!(result[0],verified[0]);
/// ```
pub fn fa_stats (
    input: &str,
) -> Vec<SeqInfo> {
    let (unflate_stream, _compression_format) = niffler::from_path(input).expect("ERROR: could not open path for reading!");
    let reader     = io::fasta::Reader::new(unflate_stream);
    //let reader = io::fasta::Reader::from_file(input).expect("ERROR: could not open file for reading!");
    let mut collection : Vec<SeqInfo> = Vec::default();
    for element in reader.records() {
        
        let record = element.expect("ERROR: could not read fa entry!");
        debug!("Reading sequence: {}",record.id());
        // if no filter, write directly
        // if blacklist and ID in list then skip
        // if whitelist and ID in list then keep
        match check_fa(&record) {
            Ok(_) => (),
            Err(r) => println!("ERROR: {:?} --> QC NOT PASSED! {:?} ",record,r),
        };
        
        let mut infos : SeqInfo = SeqInfo::default();

        let id = record.id();
        infos.id = id.to_string();
        infos.seq_length = record.seq().len() as u64;
        let mut at = 0_u64;
        let mut gc = 0_u64;
        let mut _amg = 0_u64;
   
        for letter in record.seq() {
            match &[*letter] {
                b"A" => at += 1,
                b"T" => at += 1,
                b"a" => at += 1,
                b"t" => at += 1,
                b"G" => gc += 1,
                b"C" => gc += 1,
                b"g" => gc += 1,
                b"c" => gc += 1,
                b"N" => infos.seq_bp_gaps +=1,
                b"n" => infos.seq_bp_gaps +=1,
                _ => _amg += 1,
            };
        };
        if at!=0 && gc!=0 {
            // we dont actually care that much about exact precision
            // so instead of carrying around then 20 precisions, lets just round
            infos.seq_gc = format!("{:.4}",gc as f32 /(gc+at) as f32).parse().unwrap();
        }else if at==0 && gc!=0 {
            infos.seq_gc = 1_f32;
        }else if at!=0 && gc==0 {
            infos.seq_gc = 0_f32;
        }
        collection.push(infos);
        debug!("Finished sequence: {}",record.id());
    }
    debug!("{:?}",collection);
    collection
}



/// gathers stats for fastq
/// sequnces
/// ```rust
/// use delta_biomed_seq::fq_stats;
/// use delta_biomed_seq::SeqInfo;
/// use std::convert::TryFrom;
/// use bio::io::fasta::Record;
/// use bio::io::fasta::Writer;
/// use tempfile::NamedTempFile;
/// 
///  let tmp = NamedTempFile::new().unwrap();
/// let path = tmp.path().to_str().unwrap();
/// let record = bio::io::fastq::Record::with_attrs("id_str", Some("desc"), b"TTTTTTTTTTCCCCCCCCCCAAAAAAAAAAGGGGGGGGGG",b"????????????????????????????????????????");
/// let mut writer = bio::io::fastq::Writer::to_file(path).expect("ERROR: cant write to file!");
/// writer.write_record(&record).expect("ERROR: could not write out record!");
/// writer.flush().unwrap();
/// let result = fq_stats(path);
/// let mut verified : Vec<SeqInfo> = Vec::default();
/// let mut entry : SeqInfo = SeqInfo::default();
/// entry.seq_length = 40;
/// entry.seq_gc =0.5;
/// entry.seq_qual = 30;
/// entry.id = String::from("id_str");
/// verified.push(entry);
/// assert_eq!(result[0],verified[0]);
/// ```
pub fn fq_stats (
    input: &str,
) -> Vec<SeqInfo> {
    let (unflate_stream, _compression_format) = niffler::from_path(input).expect("ERROR: could not open path for reading!");
    let reader     = io::fastq::Reader::new(unflate_stream);
    //let reader = io::fastq::Reader::from_file(input).expect("ERROR: could not open file for reading!");
    let mut collection : Vec<SeqInfo> = Vec::default();
    for element in reader.records() {
        
        let record = element.expect("ERROR: could not read fa entry!");
        debug!("Reading sequence: {}",record.id());
        // if no filter, write directly
        // if blacklist and ID in list then skip
        // if whitelist and ID in list then keep
        match check_fq(&record) {
            Ok(_) => (),
            Err(r) => println!("ERROR: {:?} --> QC NOT PASSED! {:?} ",record,r),
        };

        let mut infos : SeqInfo = SeqInfo::default();

        let id = record.id();
        infos.id = id.to_string();
        infos.seq_length = record.seq().len() as u64;
        let mut at = 0_u64;
        let mut gc = 0_u64;
        let mut _amg = 0_u64;
   
        for letter in record.seq() {
            match &[*letter] {
                b"A" => at += 1,
                b"T" => at += 1,
                b"a" => at += 1,
                b"t" => at += 1,
                b"G" => gc += 1,
                b"C" => gc += 1,
                b"g" => gc += 1,
                b"c" => gc += 1,
                b"N" => infos.seq_bp_gaps +=1,
                b"n" => infos.seq_bp_gaps +=1,
                _ => _amg += 1,
            };
        };
        if at!=0 && gc!=0 {
            // we dont actually care that much about exact precision
            // so instead of carrying around then 20 precisions, lets just round
            infos.seq_gc = format!("{:.4}",gc as f32 /(gc+at) as f32).parse().unwrap();
        }else if at==0 && gc!=0 {
            infos.seq_gc = 1_f32;
        }else if at!=0 && gc==0 {
            infos.seq_gc = 0_f32;
        }
        let mut quals : u32 = 0;
        for entry in record.qual(){
            let tmp: u32  = entry.to_string().parse().unwrap();
            quals += tmp;
        }
        // here we need to correct for 33 offset when we want phred
        infos.seq_qual = quals/(record.qual().len() as u32)-33;
        collection.push(infos);
        debug!("Finished sequence: {}",record.id());
    }
    debug!("{:?}",collection);
    collection
}


/// gathers stats for uBAM
/// sequences
/// ```rust
/// use delta_biomed_seq::ubam_stats;
/// use delta_biomed_seq::SeqInfo;
/// use delta_biomed_seq::fq2ubam_record;
/// use rust_htslib::{bam, bam::Read, bam::Writer};
/// use std::fs::File;
/// use std::io::{self, Write};
/// use tempfile::tempdir;
/// use rust_htslib::bam::record::CigarString;
/// let mut record = bam::Record::default();
/// use std::convert::TryFrom;
/// use std::fs;
/// 
/// let tmp_dir = tempdir().unwrap();
/// let tmp_path = tmp_dir.into_path();
/// let file_path = &tmp_path.join("tmp.bam");
/// let mut tmp_file = File::create(file_path).unwrap();
/// /// here, half of them are phred score 20 = 5 and half of them are phred score 40 = I
/// let record = bio::io::fastq::Record::with_attrs("id_str", Some("desc"), b"TTTTTTTTTTCCCCCCCCCCAAAAAAAAAAGGGGGGGGGG",b"55555555555555555555IIIIIIIIIIIIIIIIIIII");
/// eprintln!("{:?}",record.qual());
/// let mut header = bam::Header::default();
///    header.push_record(
///     bam::header::HeaderRecord::new(b"PG")
///     .push_tag(b"ID", &"1")
///     .push_tag(b"SM", &"unknown")
///     .push_tag(b"PL", &"unknown")
///     );
/// let mut writer = bam::Writer::from_path(file_path, &header, bam::Format::Bam ).expect("ERROR: could not open output file");
/// let mut test = fq2ubam_record(record,"test");
/// let rg = "dummy";
/// test.push_aux(b"rq", rust_htslib::bam::record::Aux::Float(0.999)).expect("ERROR: could not push read quality into BAM record");
/// eprintln!("{:?}",test.qual());
/// writer.write(&test).unwrap();
/// drop(writer);
/// // Here we calculate now the quality first from rq field
/// let result  = ubam_stats(file_path.to_str().unwrap(),false);
/// // Next we take it from the PacBio rq field
/// let result2 = ubam_stats(file_path.to_str().unwrap(),true);
/// let mut verified : Vec<SeqInfo> = Vec::default();
/// let mut entry : SeqInfo = SeqInfo::default();
/// entry.seq_length = 40;
/// entry.seq_gc =0.5;
/// entry.seq_qual = 30;
/// entry.id = String::from("id_str");
/// verified.push(entry);
/// eprintln!("result: {:?} verified: {:?}", result,verified);
/// eprintln!("result: {:?} verified: {:?}", result2,verified);
/// assert_eq!(result[0],verified[0]);
/// assert_eq!(result2[0],verified[0]);
/// assert_eq!(result2[0],result[0]);
/// ```
pub fn ubam_stats (
    input: &str,
    force_qual: bool,
) -> Vec<SeqInfo> {
    
    let threads = 10;
    debug!("reading from {}",input);
    match metadata(input) {
        Ok(_) => println!("File exists!"),
        Err(_) => panic!("File does not exist!"),
    };
    let mut reader = bam::Reader::from_path(input).unwrap();//.expect("ERROR: could not open input file");
    // generate a thread pool for the multithreaded reading and writing of BAM files
    let pool = rust_htslib::tpool::ThreadPool::new(threads).unwrap();
    reader.set_thread_pool(&pool).unwrap();

    let mut collection : Vec<SeqInfo> = Vec::default();
    //eprintln!("reading now ");
    for element in reader.records() {
        let record = element.expect("ERROR: could not read BAM entry!");
        let rq : Option<f32> = match record.aux(b"rq").is_ok(){
            true => {
                if !force_qual {
                    let result = match record.aux(b"rq").unwrap() {
                        rust_htslib::bam::record::Aux::Float(x) => x,
                        _ => panic!("ERROR: no float in RQ field encountered!"),
                    };
                    Some(result)
                }else{
                    None
                }
            },
            false => None,
            
        };

        // we just convert it to fq and use the same routine
        // I tried as well alternatively to just read the values instead but there was 
        // no difference in speed
        let trans_record = ubam2fq_record(record);
        //eprintln!("{:?}",&trans_record);
        //let regex_gc = Regex::new(r"(?-u)[GC]").expect("ERROR: could not generate AT alphabet");
        //let regex_at = Regex::new(r"(?-u)[AT]").expect("ERROR: could not generate AT alphabet");
        //let regex_gap = Regex::new(r"(?-u)[N]").expect("ERROR: could not generate AT alphabet");
        let mut infos : SeqInfo = SeqInfo::default();
        let id = trans_record.id();
        infos.id = id.to_string();
        infos.seq_length = trans_record.seq().len() as u64;
        let mut at = 0_u64;
        let mut gc = 0_u64;
        let mut _amg = 0_u64;
   
        for letter in trans_record.seq() {
            match &[*letter] {
                b"A" => at += 1,
                b"T" => at += 1,
                b"a" => at += 1,
                b"t" => at += 1,
                b"G" => gc += 1,
                b"C" => gc += 1,
                b"g" => gc += 1,
                b"c" => gc += 1,
                b"N" => infos.seq_bp_gaps +=1,
                b"n" => infos.seq_bp_gaps +=1,
                _ => _amg += 1,
            };
        };
        if at!=0 && gc!=0 {
            // we dont actually care that much about exact precision
            // so instead of carrying around then 20 precisions, lets just round
            infos.seq_gc = format!("{:.4}",gc as f32 /(gc+at) as f32).parse().unwrap();
        }else if at==0 && gc!=0 {
            infos.seq_gc = 1_f32;
        }else if at!=0 && gc==0 {
            infos.seq_gc = 0_f32;
        }
        
        // pacbio is actually cheating and caps the rq 
        // at 1 which is not correct. As this returns otherwise just
        // the max of u32 I will instead cap to phred 60 for the moment,
        // which is  99.9999% quality
        if rq.is_none() {
            let mut quals : u32 = 0;
            //let mut quals = trans_record.qual().iter().sum::<u8>() as f64;
            //eprintln!("qual: {:?}",trans_record.qual().to_vec());
            for entry in trans_record.qual(){
                let mut tmp: u32  = entry.to_string().parse().unwrap();
                tmp = tmp-33;
                quals += tmp;
                debug!("entry: {}, sum: {}",&tmp,&quals);
            }
            
            let tmp_result =  quals as u64/infos.seq_length ;
             infos.seq_qual = tmp_result.try_into().unwrap();
        }else{
            let er = rq.unwrap();
            
            let phred : f32;
            if er < 0.0 {
                phred = 10.0 ;
            }else if er > 0.999999 {
                phred = 60.0 ;
            }else{
                phred = -10_f32 * f32::log10(1_f32-er);
                //eprintln!("{}",&phred);
            }
            infos.seq_qual = phred.round() as u32;
        }
        collection.push(infos);
        debug!("Finished sequence: {}",trans_record.id());
    }
    debug!("{:?}",collection);
    collection
}



#[derive(Default,Debug,Clone,PartialEq)]
pub struct SummaryData {
    // min seq length in bp
    pub min : u64,
    // max seq length in bp
    pub max : u64,
    // mean seq length in bp
    pub mean: u64,
    // median seq length in bp
    pub median: u64,
    // gc content 
    pub gc: f32,
    // sum of all sequences in bp
    pub total_bp: u64,
    // N of all sequences
    pub total_seq: usize,
    // sum of all gaps in bp
    pub total_nbp: u64,
    // longest gap in bp
    pub longest_n: u64,
    // NX values from 0.0-1.0 in 0.1 steps
    pub nx: Vec<u64>,
    // LX values from 0.0-1.0 in 0.1 steps
    pub lx: Vec<u64>,
    // NGX values from 0.0-1.0 in 0.1 steps
    pub ngx: Option<Vec<u64>>,
    // LGX values from 0.0-1.0 in 0.1 steps
    pub lgx: Option<Vec<u64>>,
}


/// since we dont want to write this ugly part every time again,
/// lets implement a standard display instead
impl std::fmt::Display for SummaryData {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match (self.ngx.as_ref(),self.lgx.as_ref()) {
            (Some(x),Some(y)) => {
                write!(f,"gc:\t{:.3}\nmin:\t{}\nmax:\t{}\nmean:\t{}\nmedian:\t{}\ntotal_bp:\t{}\ntotal_seq:\t{}\ntotal_nbp:\t{}\nnx:\t{}\nlx:\t{}\nngx:\t{}\nlgs:\t{:?}\n",
                self.gc,self.min,self.max,self.mean,self.median,self.total_bp,self.total_seq,self.total_nbp,self.nx.iter().format(","),self.lx.iter().format(","),x.iter().format(","),y.iter().format(","))
                },
            (Some(x),None) => {
                    write!(f,"gc:\t{:.3}\nmin:\t{}\nmax:\t{}\nmean:\t{}\nmedian:\t{}\ntotal_bp:\t{}\ntotal_seq:\t{}\ntotal_nbp:\t{}\nnnx:\t{}\nlx:\t{}\nngx:\t{}\nlgs:\tNA\n",
                    self.gc,self.min,self.max,self.mean,self.median,self.total_bp,self.total_seq,self.total_nbp,self.nx.iter().format(","),self.lx.iter().format(","),x.iter().format(","))
                },
            (None,Some(y)) => {
                    write!(f,"gc:\t{:.3}\nmin:\t{}\nmax:\t{}\nmean:\t{}\nmedian:\t{}\ntotal_bp:\t{}\ntotal_seq:\t{}\ntotal_nbp:\t{}\nnx:\t{}\nlx:\t{}\nngx:\tNA\nlgs:\t{}\n",
                    self.gc,self.min,self.max,self.mean,self.median,self.total_bp,self.total_seq,self.total_nbp,self.nx.iter().format(","),self.lx.iter().format(","),y.iter().format(","))
                    },
            (None,None) => {
                    write!(f,"gc:\t{:.3}\nmin:\t{}\nmax:\t{}\nmean:\t{}\nmedian:\t{}\ntotal_bp:\t{}\ntotal_seq:\t{}\ntotal_nbp:\t{}\nnx:\t{}\nlx:\t{}\nngx:\tNA\nlgs:\tNA\n",
                    self.gc,self.min,self.max,self.mean,self.median,self.total_bp,self.total_seq,self.total_nbp,self.nx.iter().format(","),self.lx.iter().format(","))
                    },
        }
    }   
}

/// this function expects a vector with our stats/record
/// It then generates summary statistics from these records
pub fn aggregated_stats (
    input: Vec<SeqInfo>,
    output: &str,
    g_size: &Option<u64>,
    exhaustive: bool,
    infos: &VersionInfo,
) {
    // define the output file names 
    let tmp_out1 = output.to_owned();
    let tmp_out2 = output.to_owned();
    let out_1 = tmp_out1 + "_report.txt";
    let out_2: String = tmp_out2 + "_exhaustive.txt";
    let mut ouput_1 = File::create(out_1).expect("ERROR: could not generate output file!");
    // get all obtained lengths
    let mut coll_length : Vec<u64> = input.clone().into_iter().map(|x| x.seq_length).collect();
    // get all obtained gap length
    let gap_length : Vec<u64> = input.clone().into_iter().map(|x| x.seq_bp_gaps).collect();
    // get all gc contents
    let mut gcs : Vec<f32> = input.clone().into_iter().map(|x| x.seq_gc).collect();
    // sort already our sequence length -> faster for all following steps
    coll_length.sort();
    // the NGX steps, currently in 10% intervals
    let steps : Vec<f64> = vec![0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0];
    // aggregate now all results
    let mut out = SummaryData::default();
    out.total_seq = coll_length.len();
    out.total_bp  = cumsum_u64(&coll_length);
    out.total_nbp = cumsum_u64(&gap_length);
    out.median    = median_u64(&mut coll_length);
    out.mean      = mean_u64(&coll_length);
    out.min       = coll_length[0];
    out.max       = coll_length[coll_length.len()-1];
    out.gc        = median_f32(&mut gcs);
    out.nx        = steps.clone().into_iter().map(|x| nx_u64(&mut coll_length,x,None)).collect();
    out.lx        = steps.clone().into_iter().map(|x| lx_u64(&mut coll_length,x,None)).collect();
    out.ngx       = match g_size {
        Some(genome_size) => Some(
            steps.clone()
                .into_iter()
                .map(|x| nx_u64(&mut coll_length,x,Some(*genome_size)))
                .collect()
            ),
        None => None,
    };
    out.lgx     = match g_size {
        Some(genome_size) => Some(
            steps
            .into_iter()
            .map(|x| lx_u64(&mut coll_length,x,Some(*genome_size)))
            .collect()
        ),
        None => None,
    };
    // write the obtained stats thanks to our display format definition
    writeln!(ouput_1,"# program:\t{}",infos.program).unwrap();
    writeln!(ouput_1,"# version:\t{}",infos.version).unwrap();
    writeln!(ouput_1,"# author:\t{}",infos.author).unwrap();
    writeln!(ouput_1,"# command:\t{}",infos.command).unwrap();
    writeln!(ouput_1,"{}",out).unwrap();
    
    // we can as well get all stats for all sequences 
    // e.g. for plotting quality vs length or stuff like that
    if exhaustive {
        let mut ouput_2 = File::create(out_2).expect("ERROR: could not generate output file!");
        writeln!(ouput_2,"# program:\t{}",infos.program).unwrap();
        writeln!(ouput_2,"# version:\t{}",infos.version).unwrap();
        writeln!(ouput_2,"# author:\t{}",infos.author).unwrap();
        writeln!(ouput_2,"# command:\t{}",infos.command).unwrap();
        writeln!(ouput_2,"# fields:\tID\tLengths\tGC\tQual\tGap").unwrap();
        for element in input {
            writeln!(ouput_2,"{}\t{}\t{:3}\t{}\t{}",element.id,element.seq_length,element.seq_gc,element.seq_qual,element.seq_bp_gaps).unwrap();
        }
    }

    
    
}

/// converts uBAM record to fastq record.
/// Currently only keeping the minimum necessary entries.
/// Potentially to be extended in the future to keep additional
/// values in tags such as ipd values
/// Example:
/// ```
/// use delta_biomed_seq::ubam2fq_record;
/// use rust_htslib::{bam, bam::Read};
/// use rust_htslib::bam::record::CigarString;
/// use std::convert::TryFrom;
/// use std::str;
/// let mut record = bam::Record::default();
/// let cigar_str = CigarString::try_from("".as_bytes()).expect("ERROR: could not create empty CIGAR string!");
/// // cigar can only be initiated with the set command
/// // which takes as well the qname, sequence and qual
/// record.set(
///     b"m00001_010101_000001/4194579/ccs", //qname
///     Some(&cigar_str), //Cigar
///     b"ACGT", //sequence 
///     &[30,30,30,30] //qual without +33 offset
/// );
/// record.set_tid(-1);
/// // flag 4 = unmapped
/// // note:         out_record.set_unmapped() does the same
/// record.set_flags(4);
/// // normally is 0 but BAM is 0 and not 1 based as SAM
/// record.set_pos(-1);
/// let mapq = 255;
/// record.set_mapq(mapq);
/// // no set_cigar available in library
/// // no out_record.set_mqname in library
/// // normally is 0 but BAM is 0 and not 1 based as SAM
/// record.set_mpos(-1);
/// // no out_record.set_tlen in library
/// record.set_insert_size(0);
/// // the next option is unexpectedly mandatory! 
/// // see https://github.com/rust-bio/rust-htslib/issues/339
/// record.set_mtid(-1);
/// let fq_record = ubam2fq_record(record);
/// assert_eq!("????",String::from_utf8_lossy(fq_record.qual()));
/// assert_eq!("@m00001_010101_000001/4194579/ccs\nACGT\n+\n????\n", fq_record.to_string());
///```
pub fn ubam2fq_record (bam_record: bam::Record) -> bio::io::fastq::Record {
    // this took me a lot of time
    // essentially the function qual() returns not the proper
    // value but the already corrected value by 33
    // So if we want the original values we need to correct first again
    let new_qual: Vec<u8> =  bam_record.qual().iter().map(|x| x + 33).collect();
    bio::io::fastq::Record::with_attrs(
        &str::from_utf8(bam_record.qname()).unwrap(), 
        None, 
        &bam_record.seq().as_bytes(),
        &new_qual
    )
}


/// converts uBAM record to fasta record.
/// looses evidently information
/// Example:
/// ```
/// use delta_biomed_seq::ubam2fa_record;
/// use rust_htslib::{bam, bam::Read};
/// use rust_htslib::bam::record::CigarString;
/// use std::convert::TryFrom;
/// let mut record = bam::Record::default();
/// let cigar_str = CigarString::try_from("".as_bytes()).expect("ERROR: could not create empty CIGAR string!");
/// // cigar can only be initiated with the set command
/// // which takes as well the qname, sequence and qual
/// record.set(
///     b"read1", //qname
///     Some(&cigar_str), //Cigar
///     b"ACGT", //sequence 
///     b"????" //qual
/// );
/// record.set_tid(-1);
/// // flag 4 = unmapped
/// // note:         out_record.set_unmapped() does the same
/// record.set_flags(4);
/// // normally is 0 but BAM is 0 and not 1 based as SAM
/// record.set_pos(-1);
/// let mapq = 255;
/// record.set_mapq(mapq);
/// // no set_cigar available in library
/// // no out_record.set_mqname in library
/// // normally is 0 but BAM is 0 and not 1 based as SAM
/// record.set_mpos(-1);
/// // no out_record.set_tlen in library
/// record.set_insert_size(0);
/// // the next option is unexpectedly mandatory! 
/// // see https://github.com/rust-bio/rust-htslib/issues/339
/// record.set_mtid(-1);
/// //let read_id = "read1";
/// //let description = None;
/// //let sequence = b"ACGT";
/// //let record = bio::io::fasta::Record::with_attrs(read_id, description, sequence);
/// let fa_record = ubam2fa_record(record);
/// assert_eq!(">read1\nACGT\n", fa_record.to_string());
/// ```
pub fn ubam2fa_record (bam_record: bam::Record) -> bio::io::fasta::Record {
    bio::io::fasta::Record::with_attrs(
        &str::from_utf8(bam_record.qname()).unwrap(), 
        None, 
        &bam_record.seq().as_bytes()
    )
}


/// converts fasta to fastq. A dummy quality value (symbol) must be provided
/// in order to fill the needed quality track.
/// Example:
/// ```
/// use delta_biomed_seq::fq2fa_record;
/// let read_id = "read1";
/// let description = Some("desc");
/// let sequence = b"ACGT";
/// let record = bio::io::fastq::Record::with_attrs(read_id, description, sequence, b"????" );
/// let fa_record = fq2fa_record(record);
/// assert_eq!(">read1 desc\nACGT\n", fa_record.to_string());
/// ```
pub fn fq2fa_record (fq_record: bio::io::fastq::Record) -> bio::io::fasta::Record {
    bio::io::fasta::Record::with_attrs(
        fq_record.id(),
        fq_record.desc(),
        fq_record.seq()
    )
}


/// converts fasta to fastq. A dummy quality value (symbol) must be provided
/// in order to fill the needed quality track.
/// Example:
/// 
/// ```
/// use delta_biomed_seq::fa2fq_record;
/// let read_id = "read1";
/// let description = None;
/// let sequence = b"ACGT";
/// let record = bio::io::fasta::Record::with_attrs(read_id, description, sequence);
/// let fq_record = fa2fq_record(record, &"?");
/// assert_eq!("@read1\nACGT\n+\n????\n", fq_record.to_string());
/// ```
pub fn fa2fq_record (fa_record: bio::io::fasta::Record, qual:&str) -> bio::io::fastq::Record {
    bio::io::fastq::Record::with_attrs(
        fa_record.id(),
        fa_record.desc(),
        fa_record.seq(),
        qual.repeat(fa_record.seq().len()).as_bytes() 
    )
}


/// converts fastq to uBAM record.
/// Currently only keeping the minimum necessary entries.
/// This is though unaware of the header which needs to
/// be accordingly specified correctly as well in the
/// "mother" function.
/// 
/// - Field2: Flag is 4 = unmapped
/// - Field3: targetID is * = unmapped 
/// - Field4: POS is 0 in SAM as SAM is 1-based
/// - Field5: MAPQ is 255 which is specific for unmapped
/// - Field6: CIGAR is * 
/// - Field7: RNEXT is * as not mate
/// - Field8: PNEXT is 0 as not mate 
/// - Field9: TLEN is 0 as no template length
/// 
/// Example:
/// ```
/// use delta_biomed_seq::fq2ubam_record;
/// use rust_htslib::bam::record::CigarString;
/// use rust_htslib::{bam, bam::Read};
/// use std::convert::TryFrom;
/// let mut bam_record = bam::Record::default();
/// let cigar_str = CigarString::try_from("".as_bytes()).expect("ERROR: could not create empty CIGAR string!");
/// // cigar can only be initiated with the set command
/// // which takes as well the qname, sequence and qual
/// bam_record.set(
///     b"read1", //qname
///     Some(&cigar_str), //Cigar
///     b"ACGT", //sequence 
///     b"????" //qual
/// );
/// bam_record.set_tid(-1);
/// // flag 4 = unmapped
/// // note:         out_record.set_unmapped() does the same
/// bam_record.set_flags(4);
/// // normally is 0 but BAM is 0 and not 1 based as SAM
/// bam_record.set_pos(-1);
/// let mapq = 255;
/// bam_record.set_mapq(mapq);
/// // no set_cigar available in library
/// // no out_record.set_mqname in library
/// // normally is 0 but BAM is 0 and not 1 based as SAM
/// bam_record.set_mpos(-1);
/// // no out_record.set_tlen in library
/// bam_record.set_insert_size(0);
/// // the next option is unexpectedly mandatory! 
/// // see https://github.com/rust-bio/rust-htslib/issues/339
/// bam_record.set_mtid(-1);
/// let read_id = "read1";
/// let description = Some("desc");
/// let sequence = b"ACGT";
/// let fq_record = bio::io::fastq::Record::with_attrs(read_id, description, sequence, b"????" );
/// let rg = "dummy";
/// bam_record.push_aux(b"RG", rust_htslib::bam::record::Aux::String(rg)).expect("ERROR: could not push RG into BAM record");
/// let new_bam_record = fq2ubam_record(fq_record,&rg);
/// assert_eq!(&vec![30,30,30,30], new_bam_record.qual());
/// ```
pub fn fq2ubam_record (fq_record: bio::io::fastq::Record, rg: &str ) -> bam::Record {
    let mut out_record = bam::Record::default();
    //let head_view  = std::rc::Rc::new(bam::HeaderView::from_header(&header));
    //out_record.set_header(head_view);
    // now we need to populate a bit and assure everything is according
    // to uBAM specifications (using PacBio as template)
    let cigar_str = CigarString::try_from("".as_bytes()).expect("ERROR: could not create empty CIGAR string!");
    // cigar can only be initiated with the set command
    // which takes as well the qname, sequence and qual
    // unfortunately we need to convert for htslib the score to one without the 33 offset....
    let new_qual : Vec<u8> = fq_record.qual().iter().map(|x| x - 33 ).collect();
    out_record.set(
        fq_record.id().as_bytes(), //qname
        Some(&cigar_str), //Cigar
        fq_record.seq(), //sequence 
        &new_qual //qual
    );
    out_record.set_tid(-1);
    // flag 4 = unmapped
    // note:         out_record.set_unmapped() does the same
    out_record.set_flags(4);
    // normally is 0 but BAM is 0 and not 1 based as SAM
    out_record.set_pos(-1);
    let mapq = 255;
    out_record.set_mapq(mapq);
    // no set_cigar available in library
    // no out_record.set_mqname in library
    // normally is 0 but BAM is 0 and not 1 based as SAM
    out_record.set_mpos(-1);
    // no out_record.set_tlen in library
    out_record.set_insert_size(0);
    // the next option is unexpectedly mandatory! 
    // see https://github.com/rust-bio/rust-htslib/issues/339
    out_record.set_mtid(-1);
    out_record.push_aux(b"RG", rust_htslib::bam::record::Aux::String(rg)).expect("ERROR: could not push RG into BAM record");
    out_record
}


/// converts fasta to uBAM record.
/// Adding a pseudo quality for the sequence.
/// Currently only keeping the minimum necessary entries.
/// This is though unaware of the header which needs to
/// be accordingly specified correctly as well in the
/// "mother" function.
/// 
/// - Field2: Flag is 4 = unmapped
/// - Field3: targetID is * = unmapped 
/// - Field4: POS is 0 in SAM as SAM is 1-based
/// - Field5: MAPQ is 255 which is specific for unmapped
/// - Field6: CIGAR is * 
/// - Field7: RNEXT is * as not mate
/// - Field8: PNEXT is 0 as not mate 
/// - Field9: TLEN is 0 as no template length
/// 
/// Example:
/// ```
/// use delta_biomed_seq::fa2ubam_record;
/// use rust_htslib::bam::record::CigarString;
/// use rust_htslib::{bam, bam::Read};
/// use std::convert::TryFrom;
/// let mut bam_record = bam::Record::default();
/// let cigar_str = CigarString::try_from("".as_bytes()).expect("ERROR: could not create empty CIGAR string!");
/// // cigar can only be initiated with the set command
/// // which takes as well the qname, sequence and qual
/// bam_record.set(
///     b"read1", //qname
///     Some(&cigar_str), //Cigar
///     b"ACGT", //sequence 
///     b"????" //qual
/// );
/// bam_record.set_tid(-1);
/// // flag 4 = unmapped
/// // note:         out_record.set_unmapped() does the same
/// bam_record.set_flags(4);
/// // normally is 0 but BAM is 0 and not 1 based as SAM
/// bam_record.set_pos(-1);
/// let mapq = 255;
/// bam_record.set_mapq(mapq);
/// // no set_cigar available in library
/// // no out_record.set_mqname in library
/// // normally is 0 but BAM is 0 and not 1 based as SAM
/// bam_record.set_mpos(-1);
/// // no out_record.set_tlen in library
/// bam_record.set_insert_size(0);
/// // the next option is unexpectedly mandatory! 
/// // see https://github.com/rust-bio/rust-htslib/issues/339
/// bam_record.set_mtid(-1);
/// let read_id = "read1";
/// let description = Some("desc");
/// let sequence = b"ACGT";
/// let qual = "?";
/// let fa_record = bio::io::fasta::Record::with_attrs(read_id, description, sequence);
/// let rg = "dummy";
/// bam_record.push_aux(b"RG", rust_htslib::bam::record::Aux::String(rg)).expect("ERROR: could not push RG into BAM record");
/// let new_bam_record = fa2ubam_record(fa_record,qual,rg);
/// assert_eq!(new_bam_record, bam_record);
/// ```
pub fn fa2ubam_record (fa_record: bio::io::fasta::Record, qual: &str, rg: &str  ) -> bam::Record {
    let mut out_record = bam::Record::default();
    //let head_view  = std::rc::Rc::new(bam::HeaderView::from_header(&header));
    //out_record.set_header(head_view);
    // now we need to populate a bit and assure everything is according
    // to uBAM specifications (using PacBio as template)
    let cigar_str = CigarString::try_from("".as_bytes()).expect("ERROR: could not create empty CIGAR string!");
    // cigar can only be initiated with the set command
    // which takes as well the qname, sequence and qual
    out_record.set(
        fa_record.id().as_bytes(), //id
        Some(&cigar_str), //Cigar
        fa_record.seq(),  //seq
        qual.repeat(fa_record.seq().len()).as_bytes() //pseudo qual
    );
    out_record.set_tid(-1);
    // flag 4 = unmapped
    // note:         out_record.set_unmapped() does the same
    out_record.set_flags(4);
    // normally is 0 but BAM is 0 and not 1 based as SAM
    out_record.set_pos(-1);
    let mapq = 255;
    out_record.set_mapq(mapq);
    // no set_cigar available in library
    // no out_record.set_mqname in library
    // normally is 0 but BAM is 0 and not 1 based as SAM
    out_record.set_mpos(-1);
    // no out_record.set_tlen in library
    out_record.set_insert_size(0);
    // the next option is unexpectedly mandatory! 
    // see https://github.com/rust-bio/rust-htslib/issues/339
    out_record.set_mtid(-1);
    out_record.push_aux(b"RG", rust_htslib::bam::record::Aux::String(rg)).expect("ERROR: could not push RG into BAM record");
    out_record
}

/// Check the validity of a FastQ record.
///
/// # Errors
/// This function will return an `Err` if one of the following conditions is met:
/// -   The record identifier is empty.
/// -   There is a non-ASCII character found in either the sequence or quality strings.
/// -   The sequence and quality strings do not have the same length.
/// 
/// This is directly derived from the original library and modified
/// to accommodate the check of proper alphabet
///
/// # Example
///
/// ```rust
/// use bio::io::fastq::Record;
/// use delta_biomed_seq::check_fq;
/// let mut record = Record::with_attrs("id", None,  b"ATGCGGG", b"QQQQQQQQQ");
/// let actual = check_fq(&record).unwrap_err();
/// let expected = "Unequal length of sequence an qualities.";
/// assert_eq!(actual, expected);
///
/// record = Record::with_attrs("id_str", Some("desc"), b"ATGCGGG", b"QQQQQQQ");
/// assert!(check_fq(&record).is_ok());
/// record = Record::with_attrs("id_str", Some("desc"), b"ATGCWGG", b"QQQQQQQ");
/// assert!(check_fq(&record).is_ok());
/// record = Record::with_attrs("id_str", Some("desc"), b"ATGC@GGG", b"QQQQ0QQQ");
/// let actual = check_fq(&record).unwrap_err();
/// let expected = "Non valid character in seq found!";
/// assert_eq!(actual, expected);
/// ```
pub fn check_fq(fq_record: &bio::io::fastq::Record) -> Result<(), &'static str> {
    if fq_record.id().is_empty() {
        return Err("Expecting id for FastQ record.");
    }
    if !fq_record.qual().is_ascii() {
        return Err("Non-ascii character found in qualities.");
    }
    if fq_record.seq().len() != fq_record.qual().len() {
        return Err("Unequal length of sequence an qualities.");
    }
    let alphabet = alphabets::dna::alphabet();
    let iuapac   = alphabets::dna::iupac_alphabet();
    for letter in fq_record.seq(){
        if !iuapac.is_word(&[*letter]) {
            return Err("Non valid character in seq found!");
        }else if !alphabet.is_word(&[*letter]){
            eprintln!("WARNING: your fasta sequence contained IUPAC conform but non ATCGN nucleotides");
        }
    }


    Ok(())
}


/// Check the validity of a FastA record.
///
/// # Errors
/// This function will return an `Err` if one of the following conditions is met:
/// -   The record identifier is empty.
/// -   There is a non-ASCII character found in the sequence 
/// 
/// This is directly derived from the original library and modified
/// to accommodate the check of proper alphabet
///
/// # Example
///
/// ```rust
/// use bio::io::fasta::Record;
/// use delta_biomed_seq::check_fa;
/// let mut record = Record::with_attrs("id", None, "Prüfung".as_ref());
/// let actual = check_fa(&record).unwrap_err();
/// let expected = "Non valid character in seq found!";
/// assert_eq!(actual, expected);
///
/// record = Record::with_attrs("id_str", Some("desc"), b"ATGCGGG");
/// assert!(check_fa(&record).is_ok());
/// record = Record::with_attrs("id_str", Some("desc"), b"ATGSGGG");
/// assert!(check_fa(&record).is_ok());
/// record = Record::with_attrs("id_str", Some("desc"), b"ATGC@GGG");
/// let actual = check_fa(&record).unwrap_err();
/// let expected = "Non valid character in seq found!";
/// assert_eq!(actual, expected);
/// ```
pub fn check_fa(fa_record: &bio::io::fasta::Record) -> Result<(), &'static str> {
    //let alphabet = Regex::new(r"(?-u)[A-z\*\-]").expect("ERROR: could not generate personal alphabet");
    let alphabet = alphabets::dna::alphabet();
    let iuapac   = alphabets::dna::iupac_alphabet();
    if fa_record.id().is_empty() {
        return Err("Expecting id for FastQ record.");
    }
    for letter in fa_record.seq(){
        if !iuapac.is_word(&[*letter]) {
            return Err("Non valid character in seq found!");
        }else if !alphabet.is_word(&[*letter]) && letter != &b'N' {
            eprintln!("WARNING: your fasta sequence contained IUPAC conform but non ATCGN nucleotides");
            debug!("{}",std::str::from_utf8(&[*letter]).unwrap().to_owned());
        }
    }

    Ok(())
}


/// Force modify extended IUPAC nbucs
///
/// We dont like for certain applications
/// if we have non ATCGN characters.
/// Therefore, even though extended IUPAC such as S or W
/// are potential fine, we might want to convert them for other purposes.
/// We do replace them by N since one cant otherwise decide between 
/// 2 suggested alternatives
///
/// # Example
///
/// ```rust
/// use bio::io::fasta::Record;
/// use delta_biomed_seq::simple_iupac;
/// let record = Record::with_attrs("id_str", Some("desc"), b"ATGCGGG");
/// let trans_record = simple_iupac(&record);
/// assert_eq!(record,trans_record);
/// let record  = Record::with_attrs("id_str", Some("desc"), b"ATGWSWS");
/// let record2 = Record::with_attrs("id_str", Some("desc"), b"ATGNNNN");
/// let trans_record = simple_iupac(&record);
/// assert_eq!(record2,trans_record);
/// ```
pub fn simple_iupac(
    record: &bio::io::fasta::Record
) -> bio::io::fasta::Record {
    let alphabet = alphabets::dna::alphabet();
    let tmp_seq = record.seq()
        .iter()
        .map(
            |letter| {
                if alphabet.is_word(&[*letter]){
                    std::str::from_utf8(&[*letter]).unwrap().to_owned()
                }else{
                    String::from("N")
                }
            }
        )
        .collect::<String>();
    bio::io::fasta::Record::with_attrs(record.id(), record.desc(), tmp_seq.as_bytes())
}