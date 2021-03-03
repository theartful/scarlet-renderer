pub use super::sampler::Sampler;

use crate::scalar::Float;
use rand::distributions::{Distribution, Uniform};
use rand::rngs::StdRng;
use rand::SeedableRng;

pub struct UniformSampler {
    rng: StdRng,
    dist: Uniform<Float>,
}

impl UniformSampler {
    pub fn init(seed: u64) -> Self {
        let rng = StdRng::seed_from_u64(seed);
        let dist = Uniform::new(0., 1.);
        Self { rng, dist }
    }
}

impl Sampler for UniformSampler {
    fn get_1d(&mut self) -> Float {
        self.dist.sample(&mut self.rng)
    }
    fn get_2d(&mut self) -> (Float, Float) {
        (
            self.dist.sample(&mut self.rng),
            self.dist.sample(&mut self.rng),
        )
    }
    fn advance_sample(&mut self) {}
    fn advance_pixel(&mut self) {}
}
