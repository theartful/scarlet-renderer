pub use crate::scalar::Float;

// sample points from the hypercube [0, 1)^dimension
pub trait Sampler {
    fn get_1d(&mut self) -> Float;
    fn get_2d(&mut self) -> (Float, Float);
    fn advance_sample(&mut self);
    fn advance_pixel(&mut self);
}
