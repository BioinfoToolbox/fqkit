use crate::{errors::FqkitError, utils::file_reader, utils::file_writer};
use bio::io::fastq;
use log::info;

pub fn tail_n_records(
    input: Option<&String>,
    number: usize,
    rdc: bool,
    output: Option<&String>,
    compression_level: u32,
    stdout_type: char,
) -> Result<(), FqkitError> {
    let fp = fastq::Reader::new(file_reader(input)?);
    info!("get tail {} records", number);

    let mut fo = fastq::Writer::new(file_writer(output, compression_level, stdout_type)?);
    if rdc {
        let mut total = 0usize;
        for _ in fp.records() {
            total += 1;
        }
        info!("fastq file total reads number: {}", total);

        let skip_n = total - number;
        let fp2 = fastq::Reader::new(file_reader(input)?);
        for rec in fp2.records().skip(skip_n).map_while(Result::ok) {
            fo.write_record(&rec)?;
        }
    } else {
        let mut tail = vec![];
        for rec in fp.records().map_while(Result::ok) {
            tail.push(rec);
        }
        for rec in tail.iter().rev().take(number).rev() {
            fo.write_record(rec)?;
        }
    }
    fo.flush()?;

    Ok(())
}
