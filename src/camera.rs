pub use crate::ray::*;
pub use crate::vector::*;

pub struct CameraSample {
    p_film: Point2f,
}

pub trait Camera {
    fn generate_ray(&self, sample: CameraSample) -> Rayf;
}
