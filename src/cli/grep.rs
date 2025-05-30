use crate::{errors::FqkitError, utils::file_reader, utils::file_writer};
use bio::io::fastq;
use log::{error, info};
use std::io::BufRead;

pub fn grep_fastq(
    fq: Option<&String>,
    list: &String,
    full_name: bool,
    out: Option<&String>,
    compression_level: u32,
    stdout_type: char,
) -> Result<(), FqkitError> {
    let fq_reader = file_reader(fq).map(fastq::Reader::new)?;

    let mut num = 0usize;
    let mut ids = vec![];
    let fp_id = file_reader(Some(list))?;
    info!("reading reads id from file: {}", list);
    for id in fp_id.lines().map_while(Result::ok) {
        ids.push(id);
    }
    if ids.is_empty() {
        error!("{}", FqkitError::EmptyFile(list.to_string()));
        std::process::exit(1);
    }

    if let Some(out) = out {
        info!("reads write to file: {}", out);
    } else {
        info!("reads write to stdout");
    }
    let mut fo = file_writer(out, compression_level, stdout_type).map(fastq::Writer::new)?;
    for rec in fq_reader.records().map_while(Result::ok) {
        let name = if full_name {
            if let Some(desc) = rec.desc() {
                format!("{} {}", rec.id(), desc)
            } else {
                String::from(rec.id())
            }
        } else {
            String::from(rec.id())
        };
        if ids.contains(&name) {
            num += 1;
            fo.write_record(&rec)?;
        }
    }
    fo.flush()?;

    info!("total reads matched number: {}", num);
    Ok(())
}
