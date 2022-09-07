use std::{collections::BTreeMap, f32::consts::PI};

use crate::{
    cmath::CMath,
    common::*,
    encode::helpers::{order0_bits_estimate, BAND_BIN_QUANTIZE_SCALE_BASES},
};

mod helpers {
    use crate::common::NUM_BANDS;

    use super::Map;

    pub const BAND_BIN_QUANTIZE_SCALE_BASES: [u8; NUM_BANDS] = [
        200, 200, 200, 200, 200, 200, 200, 200, 198, 193, 188, 183, 178, 173,
        168, 163, 158, 153, 148, 129,
    ];

    pub fn order0_bits_estimate<Key>(freqs: &Map<Key, u32>) -> f64 {
        let num_symbols: u32 = freqs.values().sum();
        let mut bits_estimate: f64 = 0.0;
        for &freq in freqs.values() {
            let prob = freq as f64 / num_symbols as f64;
            bits_estimate += -f64::log2(prob) * freq as f64;
        }
        bits_estimate
    }
}

type Map<K, V> = BTreeMap<K, V>;

#[non_exhaustive]
pub struct EncodeResult {
    pub stream: Vec<u8>,
    pub total_bits_estimate: f64,
}

/// Encodes a raw sample stream into a newly-allocated vector.
///
/// Like `Decode`, this function expects `CosF` and `SinF` to be defined
/// by the user in the `Pulsejet::Shims` namespace before including the
/// relevant pulsejet header(s). See the documentation for `Decode` for
/// more information.
///
/// @param sampleStream Input sample stream.
/// @param sampleStreamSize Input sample stream size in samples.
/// @param sampleRate Input sample rate in samples per second (hz).
///        pulsejet is designed for 44100hz samples only, and its
///        psychoacoustics are tuned to that rate. However, other rates
///        may do something useful/interesting, so this rate is not
///        enforced, and the encoder will happily try to match a target
///        bit rate at another sample rate if desired.
/// @param targetBitRate Target bit rate in kilobits per second (kbps).
///        There's no enforced lower/upper bound, but due to codec format
///        details, the resulting bit rate will often plateau around
///        128kbps (or lower, depending on the material). ~64kbps is
///        typically transparent, ~32-64kbps is typically high quality.
///        For anything lower, it depends on the material, but it's not
///        uncommon for rates <=16kbps to actually be useful. <=0kbps
///        will usually end up around 2-3kbps.
/// @param[out] outTotalBitsEstimate Total bits estimate for the
///             encoded sample. This will typically differ slightly
///             from the actual size after compression, but on average
///             is accurate enough to be useful.
/// @return Encoded sample stream.
pub fn encode<M: CMath>(
    sample_stream: &[f32],
    sample_rate: f64,
    target_bit_rate: f64,
) -> EncodeResult {
    let mut v = vec![];
    let mut total_bits_estimate = 0.0;

    // Determine target bits/frame
    let target_bits_per_frame =
        target_bit_rate * 1000.0 * ((FRAME_SIZE as f64) / sample_rate);

    // Write out tag+version number
    v.extend_from_slice(&SAMPLE_TAG);
    v.extend_from_slice(&CODEC_VERSION_MAJOR.to_le_bytes());
    v.extend_from_slice(&CODEC_VERSION_MINOR.to_le_bytes());

    // Determine and output number of frames
    let mut num_frames =
        (sample_stream.len() as u32 + FRAME_SIZE - 1) / FRAME_SIZE;
    v.extend_from_slice(&(num_frames as u16).to_le_bytes());

    // We're going to decode one more frame than we output, so adjust the frame count
    num_frames += 1;

    // Allocate internal sample buffer including padding, fill it with silence, and copy input data into it
    let num_samples = num_frames * FRAME_SIZE;
    let num_padded_samples = num_samples + FRAME_SIZE * 2;
    let mut padded_samples = vec![0.0f32; num_padded_samples as usize];
    padded_samples[FRAME_SIZE as usize..][..sample_stream.len()].copy_from_slice(sample_stream);

    // Fill padding regions with mirrored frames from the original sample
    for i in 0..FRAME_SIZE {
        // Head padding
        padded_samples[(FRAME_SIZE - 1 - i) as usize] =
            padded_samples[(FRAME_SIZE + i) as usize];
        // Tail padding
        padded_samples[(num_padded_samples - FRAME_SIZE + i) as usize] =
            padded_samples[(num_padded_samples - FRAME_SIZE - 1 - i) as usize];
    }

    // Allocate separate streams to group correlated data
    let mut window_mode_stream: Vec<u8> = vec![];
    let mut band_energy_stream: Vec<u8> = vec![];
    let mut bin_qstream: Vec<i8> = vec![];

    // Clear quantized band energy predictions
    let mut quantized_band_energy_predictions = vec![0u8; NUM_BANDS];

    // Clear slack bits
    let mut slack_bits: f64 = 0.0;

    // Build transient frame map
    let mut is_transient_frame_map: Vec<bool> = vec![];
    let mut last_frame_energy: f32 = 0.0f32;
    for frame_index in 0..num_frames {
        // Conceptually, frames are centered around the center of each long window
        let frame_offset = FRAME_SIZE / 2 + frame_index * FRAME_SIZE;
        let mut frame_energy: f32 = 0.0f32;
        for i in 0..FRAME_SIZE {
            let sample = padded_samples[(frame_offset + i) as usize];
            frame_energy += sample * sample;
        }
        is_transient_frame_map.push(frame_energy >= last_frame_energy * 2.0f32);
        last_frame_energy = frame_energy;
    }

    // Encode frames
    for frame_index in 0..num_frames {
        // Determine and output window mode for this frame
        let is_transient_frame = is_transient_frame_map[frame_index as usize];
        let mut window_mode: WindowMode = WindowMode::Long;
        if target_bit_rate > 8.0 {
            let is_prev_frame_transient_frame = frame_index > 0
                && is_transient_frame_map[(frame_index - 1) as usize];
            let is_next_frame_transient_frame = frame_index < num_frames - 1
                && is_transient_frame_map[(frame_index + 1) as usize];
            if is_transient_frame
                || (is_prev_frame_transient_frame
                    && is_next_frame_transient_frame)
            {
                window_mode = WindowMode::Short;
            } else if is_next_frame_transient_frame {
                window_mode = WindowMode::Start;
            } else if is_prev_frame_transient_frame {
                window_mode = WindowMode::Stop;
            }
        }
        window_mode_stream.push(window_mode as u8);

        // Determine subframe configuration from window mode
        let mut num_subframes: u32 = 1;
        let mut subframe_window_offset: u32 = 0;
        let mut subframe_window_size: u32 = LONG_WINDOW_SIZE;
        if window_mode == WindowMode::Short {
            num_subframes = NUM_SHORT_WINDOWS_PER_FRAME;
            subframe_window_offset =
                LONG_WINDOW_SIZE / 4 - SHORT_WINDOW_SIZE / 4;
            subframe_window_size = SHORT_WINDOW_SIZE;
        }
        let subframe_size = subframe_window_size / 2;

        let target_bits_per_subframe =
            target_bits_per_frame / (num_subframes as f64);

        // Encode subframe(s)
        for subframe_index in 0..num_subframes {
            let mut window_bins: Vec<f32> = vec![];
            window_bins.reserve(subframe_size as usize);
            {
                // Apply window
                let frame_offset = frame_index * FRAME_SIZE;
                let window_offset =
                    subframe_window_offset + subframe_index * subframe_size;
                let mut windowed_samples: Vec<f32> = vec![];
                windowed_samples.reserve(subframe_window_size as usize);
                for n in 0..subframe_window_size {
                    let sample = padded_samples
                        [(frame_offset + window_offset + n) as usize];
                    let window =
                        mdct_window(n, subframe_window_size, window_mode);
                    windowed_samples.push(sample * window);
                }

                // Perform MDCT
                for k in 0..subframe_size {
                    let mut bin: f32 = 0.0;
                    for n in 0..subframe_window_size {
                        bin += windowed_samples[n as usize]
                            * M::cos(
                                PI / (subframe_size as f32)
                                    * (n as f32
                                        + 0.5
                                        + (subframe_size / 2) as f32)
                                    * (k as f32 + 0.5),
                            );
                    }
                    window_bins.push(bin);
                }
            }

            // Search (exhaustively) for an appropriate bin quantization scaling factor
            let mut best_quantized_band_energies: Vec<u8> = vec![];
            let mut best_band_energy_stream: Vec<u8> = vec![];
            let mut best_bin_qstream: Vec<i8> = vec![];
            let mut best_subframe_bits_estimate: f64 = 0.0;

            let min_scaling_factor: u32 = 1;
            let max_scaling_factor: u32 = 500;
            for scaling_factor in min_scaling_factor..=max_scaling_factor {
                let mut candidate_quantized_band_energies: Vec<u8> = vec![];
                let mut candidate_band_energy_stream: Vec<u8> = vec![];
                candidate_quantized_band_energies.reserve(NUM_BANDS);
                candidate_band_energy_stream.reserve(NUM_BANDS);
                let mut candidate_band_energy_freqs: Map<u8, u32> =
                    Default::default();
                let mut candidate_bin_qstream: Vec<i8> = vec![];
                candidate_bin_qstream.reserve(subframe_size as usize);
                let mut candidate_bin_qfreqs: Map<i8, u32> = Default::default();

                // Encode bands
                let mut band_bins = &window_bins[..];
                for band_index in 0..NUM_BANDS {
                    let num_bins =
                        BAND_TO_NUM_BINS[band_index] as u32 / num_subframes;

                    // Calculate band energy
                    let epsilon: f32 = 1e-27;
                    let mut band_energy: f32 = epsilon;
                    for bin_index in 0..num_bins {
                        let bin = band_bins[bin_index as usize];
                        band_energy += bin * bin;
                    }
                    band_energy = M::sqrt(band_energy);

                    // Quantize and encode band energy
                    let linear_band_energy = (f32::clamp(
                        f32::log2(band_energy / num_bins as f32),
                        -20.0f32,
                        20.0f32,
                    ) + 20.0f32)
                        / 40.0f32;
                    let quantized_band_energy =
                        (f32::round(linear_band_energy * 64.0f32)) as u8;
                    candidate_quantized_band_energies
                        .push(quantized_band_energy);
                    let quantized_band_energy_residual: u8 =
                        quantized_band_energy
                            - quantized_band_energy_predictions[band_index];
                    candidate_band_energy_stream
                        .push(quantized_band_energy_residual);
                    *candidate_band_energy_freqs
                        .entry(quantized_band_energy_residual)
                        .or_default() += 1;

                    // Determine band bin quantization scale
                    let band_bin_quantize_scale = f32::powf(
                        (BAND_BIN_QUANTIZE_SCALE_BASES[band_index]) as f32
                            / 200.0f32,
                        3.0f32,
                    ) * (scaling_factor as f32)
                        / (max_scaling_factor as f32)
                        * 127.0f32
                        * linear_band_energy
                        * linear_band_energy;

                    // Normalize, quantize, and encode band bins
                    for bin_index in 0..num_bins {
                        let bin = band_bins[bin_index as usize];
                        let bin_q = (f32::round(
                            bin / (band_energy + epsilon)
                                * band_bin_quantize_scale,
                        )) as i8;
                        candidate_bin_qstream.push(bin_q);
                        *candidate_bin_qfreqs.entry(bin_q).or_default() += 1;
                    }

                    band_bins = &band_bins[num_bins as usize..];
                }

                // Model the order 0 entropy of the quantized stream symbols in order to estimate the total bits used for encoding
                //  Also adjust estimate slightly, as squishy (and likely other compressors) tend to find additional correlations not captured by this simple model
                let band_energy_bits_estimate =
                    order0_bits_estimate(&candidate_band_energy_freqs);
                let bin_qbits_estimate =
                    order0_bits_estimate(&candidate_bin_qfreqs);
                let estimate_adjustment: f64 = 0.83;
                let subframe_bits_estimate = (band_energy_bits_estimate
                    + bin_qbits_estimate)
                    * estimate_adjustment;

                // Accept these candidate streams if this bit count estimate is closest to the target for the subframe
                let target_bits_per_subframe_with_slack_bits =
                    target_bits_per_subframe + slack_bits;
                if scaling_factor == min_scaling_factor
                    || f64::abs(
                        subframe_bits_estimate
                            - target_bits_per_subframe_with_slack_bits,
                    ) < f64::abs(
                        best_subframe_bits_estimate
                            - target_bits_per_subframe_with_slack_bits,
                    )
                {
                    best_quantized_band_energies =
                        candidate_quantized_band_energies;
                    best_band_energy_stream = candidate_band_energy_stream;
                    best_bin_qstream = candidate_bin_qstream;
                    best_subframe_bits_estimate = subframe_bits_estimate;
                }
            }

            // Update quantized band energy predictions for next subframe
            quantized_band_energy_predictions = best_quantized_band_energies;

            // Output the best-performing parameters/coefficients to their respective streams
            band_energy_stream.append(&mut best_band_energy_stream);
            bin_qstream.append(&mut best_bin_qstream);

            // Adjust slack bits depending on our estimated bits used for this subframe
            slack_bits +=
                target_bits_per_subframe - best_subframe_bits_estimate;

            // Update total bits estimate
            total_bits_estimate += best_subframe_bits_estimate;
        }
    }

    // Concatenate streams
    v.append(&mut window_mode_stream);
    v.extend(bin_qstream.into_iter().map(|x| x as u8));
    v.append(&mut band_energy_stream);

    EncodeResult {
        stream: v,
        total_bits_estimate,
    }
}
