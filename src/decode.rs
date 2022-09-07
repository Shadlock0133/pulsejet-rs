use std::{f32::consts::PI, mem::size_of};

use crate::{cmath::CMath, common::*};

trait SkipForward {
    fn skip_forward(&mut self, offset: usize);
}

impl<T> SkipForward for &[T] {
    fn skip_forward(&mut self, offset: usize) {
        *self = &self[offset..];
    }
}

/// Decodes an encoded pulsejet sample into a newly-allocated buffer.
///
/// This function is optimized for size and designed to be compiled in a
/// size-constrained environment. In such environments, it's common not
/// to have access to all of the required math functions, and instead
/// implement them by hand. For this reason, this decoder does not
/// depend on any such functions directly, and instead expects that
/// `cos`, `exp2`, `sin`, and `sqrt` functions are defined using a [`CMath`]
/// trait. Default implemention is provided through [`Std`](`crate::Std`) type.
/// pulsejet expects that these functions behave similarly
/// to the corresponding similarly-named cmath functions. This shim
/// mechanism can also be used to provide less accurate, speed-optimized
/// versions of these functions if desired.
///
/// Additionally, this function will not perform any error checking or
/// handling. The included metadata API can be used for high-level error
/// checking before decoding takes place if required (albeit not in a
/// non-size-constrained environment).
///
/// # Returns
///
/// Decoded samples in the [-1, 1] range (normalized).
pub fn decode<M: CMath>(mut input_stream: &[u8]) -> Vec<f32> {
    // Skip tag and codec version
    input_stream.skip_forward(8);

    // Read frame count, determine number of samples, and allocate
    // output sample buffer
    let mut num_frames = u16::from_le_bytes(
        input_stream[..size_of::<u16>()].try_into().unwrap(),
    ) as u32;
    input_stream.skip_forward(size_of::<u16>());
    let num_samples = num_frames * FRAME_SIZE;
    let mut samples = vec![0f32; num_samples as usize];

    // We're going to decode one more frame than we output, so adjust
    // the frame count
    num_frames += 1;

    // Set up and skip window mode stream
    let mut window_mode_stream = input_stream;
    input_stream.skip_forward(num_frames as usize);

    // Set up and skip quantized band bin stream
    let mut quantized_band_bin_stream = input_stream;
    input_stream.skip_forward((num_frames * NUM_TOTAL_BINS) as usize);

    // Allocate padded sample buffer, and fill with silence
    let num_padded_samples = num_samples + FRAME_SIZE * 2;
    let mut padded_samples = vec![0f32; num_padded_samples as usize];

    // Initialize LCG
    let mut lcg_state: u32 = 0;

    // Clear quantized band energy predictions
    let mut quantized_band_energy_predictions = vec![0u8; NUM_BANDS];

    // Decode frames
    for frame_index in 0..num_frames {
        // Read window mode for this frame
        let window_mode = match window_mode_stream[0] {
            0 => WindowMode::Long,
            1 => WindowMode::Short,
            2 => WindowMode::Start,
            3 => WindowMode::Stop,
            _ => unreachable!(),
        };
        window_mode_stream.skip_forward(1);

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

        // Decode subframe(s)
        for subframe_index in 0..num_subframes {
            // Decode bands
            let mut window_bins = vec![0f32; FRAME_SIZE as usize];
            let mut band_bins = &mut window_bins[..];
            for band_index in 0..NUM_BANDS {
                // Decode band bins
                let num_bins =
                    BAND_TO_NUM_BINS[band_index] as u32 / num_subframes;
                let mut num_nonzero_bins: u32 = 0;
                for bin_index in 0..num_bins {
                    let bin_q = quantized_band_bin_stream[0] as i8;
                    quantized_band_bin_stream.skip_forward(1);
                    if bin_q != 0 {
                        num_nonzero_bins += 1;
                    }
                    let bin = bin_q as f32;
                    band_bins[bin_index as usize] = bin;
                }

                // If this band is significantly sparse, fill in (nearly)
                // spectrally flat noise
                let bin_fill = (num_nonzero_bins as f32) / (num_bins as f32);
                let noise_fill_threshold = 0.1f32;
                if bin_fill < noise_fill_threshold {
                    let bin_sparsity = (noise_fill_threshold - bin_fill)
                        / noise_fill_threshold;
                    let noise_fill_gain = bin_sparsity * bin_sparsity;
                    for bin_index in 0..num_bins {
                        let noise_sample =
                            ((lcg_state >> 16) as i8) as f32 / 127.0f32;
                        band_bins[bin_index as usize] +=
                            noise_sample * noise_fill_gain;

                        // Transition LCG state using
                        // Numerical Recipes parameters
                        lcg_state = lcg_state
                            .wrapping_mul(1664525)
                            .wrapping_add(1013904223);
                    }
                }

                // Decode band energy
                let quantized_band_energy_residual = input_stream[0];
                input_stream.skip_forward(1);
                let quantized_band_energy: u8 =
                    quantized_band_energy_predictions[band_index]
                        .wrapping_add(quantized_band_energy_residual);
                quantized_band_energy_predictions[band_index] =
                    quantized_band_energy;
                let band_energy = M::exp2(
                    (quantized_band_energy) as f32 / 64.0f32 * 40.0f32
                        - 20.0f32,
                ) * (num_bins) as f32;

                // Normalize band bins and scale by band energy
                let epsilon: f32 = 1e-27f32;
                let mut band_bin_energy = epsilon;
                for bin_index in 0..num_bins {
                    let bin = band_bins[bin_index as usize];
                    band_bin_energy += bin * bin;
                }
                band_bin_energy = M::sqrt(band_bin_energy);
                let bin_scale = band_energy / band_bin_energy;
                for bin_index in 0..num_bins {
                    band_bins[bin_index as usize] *= bin_scale;
                }

                band_bins = &mut band_bins[num_bins as usize..];
            }

            // Apply the IMDCT to the subframe bins, then apply the appropriate
            // window to the resulting samples, and finally accumulate them into
            // the padded output buffer
            let frame_offset = frame_index * FRAME_SIZE;
            let window_offset = subframe_window_offset
                + subframe_index * subframe_window_size / 2;
            for n in 0..subframe_window_size {
                let n_plus_half = n as f32 + 0.5f32;

                let mut sample = 0.0f32;
                for k in 0..subframe_window_size / 2 {
                    sample += (2.0f32 / (subframe_window_size / 2) as f32)
                        * window_bins[k as usize]
                        * M::cos(
                            PI / (subframe_window_size / 2) as f32
                                * (n_plus_half
                                    + (subframe_window_size / 4) as f32)
                                * ((k) as f32 + 0.5f32),
                        );
                }

                let window = mdct_window(n, subframe_window_size, window_mode);
                padded_samples[(frame_offset + window_offset + n) as usize] +=
                    sample * window;
            }
        }
    }

    let size = num_samples as usize;
    samples[..size]
        .copy_from_slice(&padded_samples[FRAME_SIZE as usize..][..size]);

    samples
}
