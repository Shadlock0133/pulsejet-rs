#![allow(dead_code)]

use std::f32::consts::{FRAC_PI_2, PI};

pub const SAMPLE_TAG: [u8; 4] = *b"PLSJ";

pub const CODEC_VERSION_MAJOR: u16 = 0;
pub const CODEC_VERSION_MINOR: u16 = 1;

pub const FRAME_SIZE: u32 = 1024;
pub const NUM_SHORT_WINDOWS_PER_FRAME: u32 = 8;
pub const LONG_WINDOW_SIZE: u32 = FRAME_SIZE * 2;
pub const SHORT_WINDOW_SIZE: u32 =
    LONG_WINDOW_SIZE / NUM_SHORT_WINDOWS_PER_FRAME;

pub const NUM_BANDS: usize = 20;
pub const NUM_TOTAL_BINS: u32 = 856;

#[derive(PartialEq, Eq, Clone, Copy)]
pub enum WindowMode {
    Long = 0,
    Short = 1,
    Start = 2,
    Stop = 3,
}

pub const BAND_TO_NUM_BINS: [u8; NUM_BANDS] = [
    8, 8, 8, 8, 8, 8, 8, 8, 16, 16, 24, 32, 32, 40, 48, 64, 80, 120, 144, 176,
];

pub fn vorbis_window(n_plus_half: f32, size: u32) -> f32 {
    let sine_window = f32::sin(PI / size as f32 * n_plus_half);
    f32::sin(FRAC_PI_2 * sine_window * sine_window)
}

pub fn mdct_window(n: u32, size: u32, mode: WindowMode) -> f32 {
    let n_plus_half = (n as f32) + 0.5f32;
    if mode == WindowMode::Start {
        let short_window_offset =
            LONG_WINDOW_SIZE * 3 / 4 - SHORT_WINDOW_SIZE / 4;
        if n >= short_window_offset + SHORT_WINDOW_SIZE / 2 {
            return 0.0f32;
        } else if n >= short_window_offset {
            return 1.0f32
                - vorbis_window(
                    n_plus_half - (short_window_offset as f32),
                    SHORT_WINDOW_SIZE,
                );
        } else if n >= LONG_WINDOW_SIZE / 2 {
            return 1.0f32;
        }
    } else if mode == WindowMode::Stop {
        let short_window_offset = LONG_WINDOW_SIZE / 4 - SHORT_WINDOW_SIZE / 4;
        if n < short_window_offset {
            return 0.0f32;
        } else if n < short_window_offset + SHORT_WINDOW_SIZE / 2 {
            return vorbis_window(
                n_plus_half - (short_window_offset as f32),
                SHORT_WINDOW_SIZE,
            );
        } else if n < LONG_WINDOW_SIZE / 2 {
            return 1.0f32;
        }
    }
    vorbis_window(n_plus_half, size)
}
