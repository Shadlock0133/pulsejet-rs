pub trait CMath {
    fn cos(x: f32) -> f32;
    fn exp2(x: f32) -> f32;
    fn sin(x: f32) -> f32;
    fn sqrt(x: f32) -> f32;
}

pub struct Std;

impl CMath for Std {
    fn cos(x: f32) -> f32 {
        f32::cos(x)
    }

    fn exp2(x: f32) -> f32 {
        f32::exp2(x)
    }

    fn sin(x: f32) -> f32 {
        f32::sin(x)
    }

    fn sqrt(x: f32) -> f32 {
        f32::sqrt(x)
    }
}
