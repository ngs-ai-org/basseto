use bio::bio_types::sequence::SequenceRead;
use clap::{app_from_crate,crate_name,crate_description,crate_authors,crate_version,Arg, SubCommand};
extern crate bio;
use bio::io;
//use rust_htslib::bam::record::Seq;
use core::panic;
use std::convert::TryInto;
use std::{env};
//extern crate bambam;
extern crate rust_htslib;
use rust_htslib::{bam, bam::Read};
use std::str;
use rustc_hash::FxHashMap;
use std::path::Path;
use std::io::{BufRead, BufReader};
use std::fs::File;
extern crate delta_biomed_seq;
use delta_biomed_seq::{*};
use bio::alphabets;

use log::debug;




/// this function executes the format conversion function
/// it expects the arguments necessary for its functionality and
/// returns it success or error to the main function
fn run_conversion (matches: &clap::ArgMatches , arg_string: &str ) {
    let input_file     = matches.value_of("INPUT").expect("ERROR: no input file correctly provided!");
    let input_format   = matches.value_of("INPUTF").expect("ERROR: no input file format correctly provided!");
    let read_group     = matches.value_of("RG").expect("ERROR: could not get RG information!");
    let pseudo_qual    = matches.value_of("QUAL").expect("ERROR: could not parse pseudo-qual value!");
    let output_file    = matches.value_of("OUTPUT").expect("ERROR: no output file correctly provided!");
    let output_format  = matches.value_of("OUTPUTF").expect("ERROR: no output file format correctly provided!");
    let threads        = matches.value_of("THREADS").unwrap().parse::<u32>().expect("ERROR: could not parse threads option correctly!");
    let trash_file     = matches.value_of("TRASH").or(None);
    let rq_tmp         = matches. value_of("RQ").or(None);
    let force_upper = matches.is_present("UPPER");
    let force_simple = matches.is_present("SIMPLE");

    let rq : Option<f32> = rq_tmp.map(|x| x.parse::<f32>().expect("ERROR: could not parse read quality option correctly!"));
    let whitelist  : Option<&str>    = matches.value_of("WHITE").or(None);
    let blacklist  : Option<&str>    = matches.value_of("BLACK").or(None);
    let filtering = SubSet{
        black : blacklist,
        white : whitelist
    };
    let i_format = match input_format{
        "fasta" => InOutFormat::Fasta ,
        "fastq" => InOutFormat::Fastq ,
        "uBAM"  => InOutFormat::Ubam,
        _       => panic!("ERROR: input format not recognized!")
    };
   
    let o_format = match output_format{
        "fasta" => InOutFormat::Fasta ,
        "fastq" => InOutFormat::Fastq ,
        "uBAM"  => InOutFormat::Ubam,
        _       => panic!("ERROR: input format not recognized!")
    };

    debug!("INFO: Input format : {:?} , Output format: {:?}", i_format, o_format);
    
    match (i_format,o_format) {
        (InOutFormat::Fasta,InOutFormat::Fasta) => {
            if rq.is_some(){panic!("ERROR: RQ subsetting not applicable!")};
            fa2fa(
                input_file,
                output_file,
                filtering,
                force_upper,
                force_simple
            )
        },
        (InOutFormat::Fasta,InOutFormat::Fastq) => {
            if rq.is_some(){panic!("ERROR: RQ subsetting not applicable !")};
            fa2fq(
                input_file,
                output_file,
                pseudo_qual,
                filtering
            )
        }  ,
        (InOutFormat::Fasta,InOutFormat::Ubam)  => {
            if rq.is_some(){panic!("ERROR: RQ subsetting not applicable !")};
            fa2ubam(
            input_file,
            output_file,
            pseudo_qual,
            read_group,
            filtering,
            arg_string,
            crate_name!(),
            crate_version!()
            )
        },
        (InOutFormat::Fastq,InOutFormat::Fastq) => {
            fq2fq(
                input_file,
                output_file,
                filtering,
                rq
            )
        },
        (InOutFormat::Fastq,InOutFormat::Fasta) => {
            if rq.is_some(){panic!("ERROR: RQ subsetting not applicable !")};
            fq2fa(
            input_file,
            output_file,
            filtering
            )
        },
        (InOutFormat::Fastq,InOutFormat::Ubam)  => {
            if rq.is_some(){panic!("ERROR: RQ subsetting not applicable !")};
            fq2ubam(
            input_file,
            output_file,
            read_group,
            filtering,
            arg_string,
            crate_name!(),
            crate_version!()
            )
        },
        (InOutFormat::Ubam,InOutFormat::Ubam)  => ubam2ubam(
            input_file, 
            output_file,
            filtering,
            threads,
            arg_string,
            crate_name!(),
            crate_version!(),
            rq,
            trash_file,
        ),
        (InOutFormat::Ubam,InOutFormat::Fastq) => {
            if trash_file.is_some(){panic!("ERROR: keeping non RQ reads not applicable ")};
            ubam2fq(
            input_file,
            output_file,
            filtering,
            threads,
            rq
            )
        },
        (InOutFormat::Ubam,InOutFormat::Fasta) => {
            if trash_file.is_some(){panic!("ERROR: keeping non RQ reads not applicable ")};
            ubam2fa(
            input_file,
            output_file,
            filtering,
            threads,
            rq
            )
        },
    };
    
}

/// this function executes the statistical analysis function
/// it expects the arguments necessary for its functionality and
/// returns it success or error to the main function
fn run_stats (matches: &clap::ArgMatches , arg_string: &str )  {
    let input_file     = matches.value_of("INPUT").expect("ERROR: no input file correctly provided!");
    let output_suffix    = matches.value_of("OUTPUT").expect("ERROR: no output suffix provided!");
    let input_format   = matches.value_of("INPUTF").expect("ERROR: no input file format correctly provided!");
    let gsize_tmp   = matches.value_of("GSIZE").or(None);
    // here we convert as well immeditately as well the genome size in absolute base-pairs instead of Mbp
    let gsize : Option<u64>        = gsize_tmp.map(|x| x.parse::<u64>().expect("ERROR: could not genome size option correctly!")*1000000);
    let exhaustive = matches.is_present("EXT");
    let force_qual = matches.is_present("RQ");
    let infos = &VersionInfo {
        program: &String::from("basetto"),
        version: crate_version!(),
        author: crate_authors!(),
        command: arg_string,
    };

    let i_format = match input_format{
        "fasta" => InOutFormat::Fasta ,
        "fastq" => InOutFormat::Fastq ,
        "uBAM"  => InOutFormat::Ubam,
        _       => panic!("ERROR: input format not recognized!")
    };
    debug!("INFO: Input format : {:?} ", i_format);

    debug!("INFO: Input format : {:?} ", i_format);
    
    let seq_infos = match i_format {
        InOutFormat::Fasta => {
            fa_stats(
                input_file            )
        },
        InOutFormat::Fastq => {
            fq_stats(
                input_file            )
        },
        
        InOutFormat::Ubam => {
            ubam_stats(
                 input_file,
                 force_qual
            )
        }
    };
    aggregated_stats(seq_infos,output_suffix,&gsize, exhaustive,infos);
}


fn main() {
    env_logger::init();
    // now the next is not really for any argument 
    // parsing but simply to get the command which 
    // was used to execute as I cant get this from clap
    let args: Vec<String> = env::args().collect();
    let args_string = args.join(" ");
    let matches = app_from_crate!()
        .about(crate_description!())
        .arg(Arg::with_name("BAMBAM")
        .short("x")
        .long("bambam")
        .value_name("IDH")
        .takes_value(false)
        .required(false)
        .hidden(true))	
        .subcommand(SubCommand::with_name("convert")
            .about("convert between formats, subset and correct")
            .arg(Arg::with_name("INPUT")
                .short("i")
                .long("input")
                .value_name("FILE")
                .help("input sequences")
                .takes_value(true)
                .required(true))
            .arg(Arg::with_name("THREADS")
                .short("t")
                .long("threads")
                .value_name("INT")
                .help("number of threads for reading and writing\
                Note: beyond 20 IO will throttle likely")
                .takes_value(true)
                .default_value("1"))
            .arg(Arg::with_name("OUTPUT")
                .short("o")
                .long("output")
                .value_name("FILE")
                .help("the output sequences")
                .takes_value(true)
                .required(true))
            .arg(Arg::with_name("INPUTF")
                .short("n")
                .long("input-format")
                .value_name("STRING")
                .help("format of the input sequences")
                .possible_values(&["fasta","fastq","uBAM"])
                .takes_value(true)
                .required(true))
            .arg(Arg::with_name("OUTPUTF")
                .short("m")
                .long("output-format")
                .value_name("STRING")
                .help("format of the output sequences")
                .takes_value(true)
                .possible_values(&["fasta","fastq","uBAM"])
                .required(true))
            .arg(Arg::with_name("BAMBAM")
               .short("x")
               .long("bambam")
               .takes_value(false)
               .required(false)
               .hidden(true))
            .arg(Arg::with_name("UPPER")
                .short("s")
                .long("uppercase")
                .takes_value(false)
                .required(false)
                .help("force nucleotide sequence to upper case"))
            .arg(Arg::with_name("SIMPLE")
                .short("k")
                .long("simple-nuc")
                .takes_value(false)
                .required(false)
                .help("replaces extended IUPAC letters by N"))
            .arg(Arg::with_name("QUAL")
                .short("q")
                .long("pseudoqual")
                .value_name("SYMBOL")
                .takes_value(true)
                .default_value("?")
                .help("dummy phred score in case of fasta-in available symbols :\
                    https://support.illumina.com/help/BaseSpace_OLH_009008/Content/Source/Informatics/BS/QualityScoreEncoding_swBS.html")
                .required(false))
            .arg(Arg::with_name("RG")
                .short("r")
                .long("read-group")
                .value_name("STRING")
                .takes_value(true)
                .default_value("ID:1\tSM:unkown\tPL:unkown")
                .help("read group added in header and to each read for x2uBAM. Note: not for uBAM2uBAM")
                .required(false))
            .arg(Arg::with_name("RQ")
                .short("c")
                .long("quality")
                .value_name("FLOAT")
                .takes_value(true)
                .help("subsetting uBAM and fastq based on quality field (uBAM->uBAM; fastq->fastq")
                .conflicts_with_all(&["QUAL","RG"])
                .required(false))
            .arg(Arg::with_name("strict")
                .short("l")
                .long("strict")
                .help("strict concerning identifiers taking 1:1 string, otherwise stopping after space"))
            .arg(Arg::with_name("BLACK")
                .short("b")
                .long("blacklist")
                .value_name("FILE")
                .help("list of sequence identifiers to be ignored")
                .takes_value(true)  
                .required(false))
            .arg(Arg::with_name("WHITE")
                .short("w")
                .long("whitelist")
                .value_name("FILE")
                .help("list of sequence identifiers to be used, all others being ignored")
                .takes_value(true)
                .required(false))
            .arg(Arg::with_name("TRASH")
                .short("u")
                .long("unfiltered")
                .value_name("FILE")
                .help("will write unselected entries in 2nd file with same out format \
                currently only supported for uBAM-->uBAM")
                .takes_value(true)
                .required(false)))
        .subcommand(SubCommand::with_name("stats")
            .about("provides statistics for the input sequences ")
            .arg(Arg::with_name("INPUT")
                .short("i")
                .long("input")
                .value_name("FILE")
                .help("input sequences")
                .takes_value(true)
                .required(true))
            .arg(Arg::with_name("OUTPUT")
                .short("o")
                .long("output")
                .value_name("PREFIX")
                .help("output prefix")
                .takes_value(true)
                .required(true))
            .arg(Arg::with_name("INPUTF")
                .short("n")
                .long("input-format")
                .value_name("STRING")
                .help("format of the input sequences")
                .possible_values(&["fasta","fastq","uBAM"])
                .takes_value(true)
                .required(true))
            .arg(Arg::with_name("GSIZE")
                .short("g")
                .long("genomeSize")
                .takes_value(true)
                .required(false)
                .help("provide the genome size for NG50 in Mbp"))
            .arg(Arg::with_name("RQ")
                .short("rq")
                .long("read-quality")
                .takes_value(false)
                .required(false)
                .help("forces to calculate quality from read in presence of PacBio rq field (takes substantionally longer!)"))
            .arg(Arg::with_name("EXT")
                .short("e")
                .long("exhaustive")
                .takes_value(false)
                .required(false)
                .help("generate additional files including all values/sequence ")))
        .get_matches();
		
    
    if let Some(matches) = matches.subcommand_matches("convert") {
        debug!("INFO: Converting formats");
        run_conversion(matches,&args_string);

    }
    if let Some(matches) = matches.subcommand_matches("stats") {
        debug!("INFO: Obtaining stats");
        run_stats(matches,&args_string);
    }
}
