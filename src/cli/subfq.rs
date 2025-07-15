use super::misc::write_record;
use crate::{errors::FqkitError, utils::file_reader, utils::file_writer};
use log::{error, info};
use paraseq::{fastq, fastx::Record};
use rand::{Rng, prelude::*};
use rand_pcg::Pcg64;

pub fn sample_fastq(
    file: Option<&String>,
    n: usize,
    seed: u64,
    two_pass: bool,
    out: Option<&String>,
    compression_level: u32,
    stdout_type: char,
) -> Result<(), FqkitError> {
    if n == 0 {
        error!("n must be greater than 0");
        std::process::exit(1);
    }
    info!("subseq number: {}", n);

    let mut fq_reader = fastq::Reader::new(file_reader(file)?);
    let mut rset = fastq::RecordSet::default();
    info!("rand seed: {}", seed);
    let mut fq_writer = file_writer(out, compression_level, stdout_type)?;

    let mut rng: rand_pcg::Lcg128Xsl64 = Pcg64::seed_from_u64(seed);
    let mut order: usize = 0;
    if two_pass {
        info!("two pass mode enabled");
        let mut get: Vec<usize> = Vec::with_capacity(n);
        while rset.fill(&mut fq_reader)? {
            for _ in rset.iter().map_while(Result::ok) {
                if order < n {
                    get.push(order);
                } else {
                    let ret = rng.random_range(0..=order);
                    if ret < n {
                        get[ret] = order;
                    }
                }
                order += 1;
            }
        }

        let mut fq_reader2 = fastq::Reader::new(file_reader(file)?);
        let mut rset2 = fastq::RecordSet::default();
        order = 0;
        get.sort_unstable(); // keep the order
        let mut idx = 0usize;
        info!("all records has been readed into memory, start write to output ...");
        while rset2.fill(&mut fq_reader2)? {
            for rec in rset2.iter().map_while(Result::ok) {
                if idx < get.len() && order == get[idx] {
                    write_record(&mut fq_writer, rec.id(), rec.seq(), rec.qual())?;
                    idx += 1;
                }
                if idx >= get.len() {
                    break;
                }
                order += 1;
            }
        }
    } else {
        let mut get = Vec::with_capacity(n);
        while rset.fill(&mut fq_reader)? {
            for rec in rset.iter().map_while(Result::ok) {
                if order < n {
                    get.push((
                        order,
                        rec.id_str().to_owned(),
                        rec.seq_str().to_owned(),
                        rec.qual_str().to_owned(),
                    ));
                } else {
                    let ret = rng.random_range(0..=order);
                    if ret < n {
                        get[ret] = (
                            order,
                            rec.id_str().to_owned(),
                            rec.seq_str().to_owned(),
                            rec.qual_str().to_owned(),
                        );
                    }
                }
                order += 1;
            }
        }
        info!("all records has been readed into memory, start write to output ...");
        get.sort_unstable_by_key(|x| x.0); // sort by order to keep the raw order
        for (_, id, seq, qual) in get {
            write_record(
                &mut fq_writer,
                id.as_bytes(),
                seq.as_bytes(),
                qual.as_bytes(),
            )?;
        }
    }

    fq_writer.flush()?;

    Ok(())
}
