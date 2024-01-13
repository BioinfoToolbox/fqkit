use anyhow::{Ok, Error};
use bio::io::fastq;
use std::collections::HashMap;
use log::*;
use plotters::prelude::*;
use lowcharts::plot;
use crate::utils::*;
use std::time::Instant;

pub fn gc_content(
    fqin: &Option<&str>,
    output: &Option<&str>,
    show: bool,
    prefix: String, 
    width: usize, 
    height: usize,
    ylim: usize,
    types: &str,
    compression_level: u32,
) -> Result<(),Error> {
    if let Some(inp) = fqin {
        info!("reading from file: {}", inp);
    } else {
        info!("reading from stdin");
    }
    let start = Instant::now();

    let fq_reader = file_reader(fqin).map(fastq::Reader::new)?;
    let mut fo = file_writer(output, compression_level)?;
    let mut df_hash = HashMap::new() ;

    for rec in fq_reader.records().flatten() {
        let gc_count = rec.seq().iter().filter(|x| *x ==&b'G'|| *x==&b'C').count();
        let gc_ratio = (gc_count as f64 / rec.seq().len() as f64 *100.0).round() as u64;
        *df_hash.entry(gc_ratio).or_insert(0) +=1usize;
       
    }

    fo.write("GC(%)\tReads\tRatio(%)\n".as_bytes())?;
    let mut df_ret = vec![]; // data for PNG / SVG
    let mut df_num = vec![]; // data for histogram in terminal
    let total = df_hash.values().sum::<usize>() as f32;

    for i in 0..=100 {
        let num = *df_hash.get(&i).unwrap_or(&0);
        let v = (num as f32 *10000.0 / total ).round() / 100.0;
        df_ret.push(v);
        fo.write(format!("{}\t{}\t{}\n", i, num, v).as_bytes())?;
        if show {
            for _ in 0..num{
                df_num.push(i as f64)
            }
        }
    }
    fo.flush()?;
    
    //plot_gc(df_ret, prefix, width, height, ylim, types, quiet)?;
    if show {
        info!("{}",plot::Histogram::new(&df_num, plot::HistogramOptions { intervals: 20, ..Default::default() }));
    }
    plot_gc(df_ret, prefix, width, height, ylim, types)?;

    info!("time elapsed is: {:?}",start.elapsed());
    Ok(())
}

// GC content line plot
fn plot_gc(
    data: Vec<f32>, 
    prefix: String, 
    width: usize, 
    height: usize,
    ylim: usize,
    types: &str,
) -> Result<(),Error> {
    if !["svg", "png"].contains(&types) {
        error!("invalid args types.");
        std::process::exit(1);
    }
    if ylim > 100 {
        error!("invalid args ylim.");
        std::process::exit(1);
    }
    let name = if types == "png" {format!("{}.png",prefix)} else {format!("{}.svg",prefix)};
    info!("output gc content plot: {}",name);

    if types == "png" {
        let png = BitMapBackend::new(&name, (width as u32, height as u32)).into_drawing_area();
        png.fill(&WHITE)?;

        let mut charts = ChartBuilder::on(&png)
            .margin(10)
            .caption("GC distrbution plot", ("sans-serif", 30).into_font())
            .x_label_area_size(40)
            .y_label_area_size(40)
            //.build_cartesian_2d(0..100, 0..ylim)?;
            .build_cartesian_2d(0.0..100f32 , 0.0..ylim as f32)?;

        charts
            .configure_mesh()
            .x_labels(20)
            .x_desc("GC CONTENT")
            .x_label_formatter(&|x| format!("{:.0}%",x))
            .y_labels(20)
            .y_label_formatter(&|x| format!("{:.0}%", x))
            .y_desc("percent")
            .draw()?;

        charts.draw_series(
            AreaSeries::new((0..).zip(data.iter()).map(|(x,y)| (x as f32, *y)), 0., &GREEN.mix(0.5)).border_style(&GREEN),
        )?
            .label("GC content")
            .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &GREEN));

        //charts.draw_series( LineSeries::new((0..).zip(data.iter()).map(|(x,y)| (x as f32, *y)), &BLACK) )?;

        charts
            .configure_series_labels()
            .background_style(&WHITE.mix(0.9))
            .border_style(&BLACK)
            .position(SeriesLabelPosition::UpperRight)
            .draw()?;

    } else {
        let svg = SVGBackend::new(&name, (width as u32, height as u32)).into_drawing_area();
        svg.fill(&WHITE)?;

        let mut charts = ChartBuilder::on(&svg)
            .margin(10)
            .caption("GC distrbution plot", ("sans-serif", 30).into_font())
            .x_label_area_size(40)
            .y_label_area_size(40)
            .build_cartesian_2d(0.0..100f32 , 0.0..ylim as f32)?;

        charts
            .configure_mesh()
            .x_labels(20)
            .x_desc("GC CONTENT")
            .x_label_formatter(&|x| format!("{:.0}%",x))
            .y_labels(20)
            .y_label_formatter(&|x| format!("{:.0}%", x))
            .y_desc("percent")
            .draw()?;

        charts.draw_series(
            AreaSeries::new((0..).zip(data.iter()).map(|(x,y)| (x as f32, *y)), 0., &GREEN.mix(0.5)).border_style(&GREEN),
        )?
            .label("GC content")
            .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &GREEN));

        charts
            .configure_series_labels()
            .background_style(&WHITE.mix(0.9))
            .border_style(&BLACK)
            .position(SeriesLabelPosition::UpperRight)
            .draw()?;
    }
    Ok(())
}