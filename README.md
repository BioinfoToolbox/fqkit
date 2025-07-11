![icon](https://github.com/sharkLoc/fqkit/blob/main/doc/fqkit_icon.PNG)
<!-- ![icon](doc/fqkit_icon.PNG) -->

# fqkit
![Static Badge](https://img.shields.io/badge/Author-sharkLoc-blue)
![Static Badge](https://img.shields.io/badge/Tool-fqkit-red)
![Crates.io (latest)](https://img.shields.io/crates/dv/fqkit?labelColor=rgb&color=hex&link=https%3A%2F%2Fcrates.io%2Fcrates%2Ffqkit)
![Crates.io](https://img.shields.io/crates/d/fqkit?label=Total%20download%20in%20crate.io)
![GitHub Gist last commit](https://img.shields.io/github/gist/last-commit/a4910923a230b8975218a188528463d7?logo=github)



🦀 a simple program for fastq file manipulation

## new
Starting from version 0.4.13, the paraseq library is used to parse FASTQ files, significantly improving data processing speed.


## install
##### setp1： install cargo first 
```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
```

##### step2:  on linux or windows
```bash
cargo install fqkit
# or

git clone https://github.com/sharkLoc/fqkit.git
cd fqkit
cargo b --release
# mv target/release/fqkit to anywhere you want 
```
##### install latest version

```bash
cargo install --git https://github.com/sharkLoc/fqkit.git
```

## usage

```bash
FqKit -- A simple and cross-platform program for fastq file manipulation

Version: 0.4.14

Authors: sharkLoc <mmtinfo@163.com> <mmtinfo@163.com>
Source code: https://github.com/sharkLoc/fqkit.git

Fqkit supports reading and writing gzip (.gz) format.
Bzip2 (.bz2) format is supported since v0.3.8.
Xz (.xz) format is supported since v0.3.9.
Under the same compression level, xz has the highest compression ratio but consumes more time. 

Compression level:
  format   range   default   crate
  gzip     1-9     6         https://crates.io/crates/flate2
  bzip2    1-9     6         https://crates.io/crates/bzip2
  xz       1-9     6         https://crates.io/crates/xz2
  zstd     1-4     2         roughly equals to zstd 1, 3, 7, 11, respectively


Usage: fqkit.exe [OPTIONS] <COMMAND>

Commands:
  topn     get first N records from fastq file [aliases: head]
  tail     get last N records from fastq file
  concat   concat fastq files from different lanes
  subfq    subsample sequences from big fastq file [aliases: sample]
  select   select pair-end reads by read id
  trim     trim fastq reads by position
  adapter  cut the adapter sequence on the reads
  filter   a simple filter for pair end fastq sqeuence
  join     join paired end reads that are overlapping into a single longer read
  range    print fastq records in a range
  search   search reads/motifs from fastq file
  grep     grep fastq sequence by read id or full name
  stats    summary for fastq format file [aliases: stat]
  kmer     a simple kmer counter
  shuffle  shuffle fastq sequences
  size     report the number sequences and bases
  slide    extract subsequences in sliding windows
  sort     sort fastq file by name/seq/gc/length
  plot     line plot for A T G C N percentage in read position
  fq2fa    translate fastq to fasta
  fq2sam   converts a fastq file to an unaligned SAM file
  fqscore  converts the fastq file quality scores
  flatten  flatten fastq sequences [aliases: flat]
  barcode  perform demultiplex for pair-end fastq reads [aliases: demux]
  remove   remove reads by read name [aliases: rm]
  rename   rename sequence id in fastq file [aliases: rn]
  reverse  get a reverse-complement of fastq file [aliases: rev]
  split    split interleaved fastq file
  merge    merge PE reads as interleaved fastq file
  mask     convert any low quality base to 'N' or other chars
  split2   split fastq file by records number
  gcplot   get GC content result and plot [aliases: gc]
  length   get reads length count [aliases: len]
  view     view fastq file page by page
  help     Print this message or the help of the given subcommand(s)

Global Arguments:
  -@, --threads <INT>            threads number [default: 4]
      --compress-level <INT>     set gzip/bzip2/xz/zstd compression level 1 (compress faster) - 9 (compress better) for output file, just work with option -o/--out [default: 6]
      --output-type <u|g|b|x|z>  output type for stdout: 'g' gzip; 'b' bzip2; 'x' xz; 'z' zstd; 'u' uncompressed txt format [default: u]
      --log <FILE>               if file name specified, write log message to this file, or write to stderr
  -v, --verbosity...             control verbosity of logging, [-v: Error, -vv: Warn, -vvv: Info, -vvvv: Debug, -vvvvv: Trace, defalut: Debug]

Global FLAGS:
  -q, --quiet    be quiet and do not show any extra information
  -h, --help     prints help information
  -V, --version  prints version information

Use "fqkit help [command]" for more information about a command
```

#### ** any bugs please report issues **💖
