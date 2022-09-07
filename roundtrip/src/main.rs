use std::{
    fs::{self, File},
    io::BufWriter,
    path::PathBuf,
};

mod util;
use clap::Parser;
use util::*;

use hound::{SampleFormat, WavReader, WavSpec, WavWriter};
use pulsejet_rs::Std;

#[derive(Parser)]
enum Opts {
    Encode {
        #[clap(short, long)]
        bit_rate: f64,
        input: PathBuf,
        output: PathBuf,
    },
    Decode {
        input: PathBuf,
        output: PathBuf,
    },
}

fn main() {
    eprintln!("start");
    match Opts::parse() {
        Opts::Encode {
            bit_rate,
            input,
            output,
        } => {
            let reader = WavReader::open(input).unwrap();
            let spec = reader.spec();
            assert_eq!(spec.channels, 2);
            assert_eq!(spec.sample_rate, 44100);
            let max_value = 1i32
                .checked_shl(spec.bits_per_sample as u32)
                .and_then(|x| x.checked_sub(1))
                .unwrap_or(i32::MAX);
            let samples: Vec<f32> = match spec.sample_format {
                SampleFormat::Float => reader
                    .into_samples::<f32>()
                    .array_chunks::<2>()
                    .map(|[a, b]| Ok::<_, hound::Error>((a? + b?) / 2.0))
                    .collect::<Result<_, _>>(),
                SampleFormat::Int => reader
                    .into_samples::<i32>()
                    .map(|x| x.map(|x| x as f32 / max_value as f32))
                    .array_chunks::<2>()
                    .map(|[a, b]| Ok((a? + b?) / 2.0))
                    .collect::<Result<_, _>>(),
            }
            .unwrap();
            let result =
                pulsejet_rs::encode::<Std>(&samples, 44100.0, bit_rate);
            eprintln!("Total bits estimate: {}", result.total_bits_estimate);
            let encoded = result.stream;
            fs::write(output, &encoded).unwrap();
        }
        Opts::Decode { input, output } => {
            let encoded = fs::read(input).unwrap();
            let decoded = pulsejet_rs::decode::<Std>(&encoded);
            let output = BufWriter::new(File::create(output).unwrap());
            let mut writer = WavWriter::new(
                output,
                WavSpec {
                    channels: 1,
                    sample_rate: 44100,
                    bits_per_sample: 32,
                    sample_format: SampleFormat::Float,
                },
            )
            .unwrap();
            for sample in decoded {
                writer.write_sample(sample).unwrap();
            }
            writer.finalize().unwrap();
            eprintln!("finished");
        }
    }
}
