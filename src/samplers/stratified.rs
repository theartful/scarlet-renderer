pub use super::sampler::Sampler;
pub use crate::scalar::Float;
use rand::distributions::{Distribution, Uniform};
use rand::rngs::StdRng;
use rand::{RngCore, SeedableRng};

// based on the paper "Correlated Multi-Jittered Sampling" from pixar
// https://graphics.pixar.com/library/MultiJitteredSampling/paper.pdf

fn permute_index(index: u32, len: u32, key: u32) -> u32 {
    let mut w = len - 1;
    w |= w >> 1;
    w |= w >> 2;
    w |= w >> 4;
    w |= w >> 8;
    w |= w >> 16;

    let mut i = index;
    loop {
        i ^= key;
        i *= 0xe170893d;
        i ^= key >> 16;
        i ^= (i & w) >> 4;
        i ^= key >> 8;
        i *= 0x0929eb3f;
        i ^= key >> 23;
        i ^= (i & w) >> 1;
        i *= 1 | key >> 27;
        i *= 0x6935fa69;
        i ^= (i & w) >> 11;
        i *= 0x74dcb303;
        i ^= (i & w) >> 2;
        i *= 0x9e501cc3;
        i ^= (i & w) >> 2;
        i *= 0xc860a3df;
        i &= w;
        i ^= i >> 5;

        if i < len {
            break;
        }
    }
    (i + key) % len
}

pub struct StratifiedSampler {
    // sample_count should be a square number
    sample_count: u32,
    sample_count_sqrt: u32,
    current_sample: u32,
    current_dimension: u32,
    permutation_key: u32,
    rng: StdRng,
    jitter_dist: Uniform<Float>,
}

impl StratifiedSampler {
    pub fn init(sample_count: u32, seed: u64) -> Self {
        let mut rng = StdRng::seed_from_u64(seed);
        let permutation_key = rng.next_u32();
        let jitter_dist = Uniform::new(0., 1.);
        let sample_count_sqrt = (sample_count as Float).sqrt().round() as u32;

        Self {
            sample_count,
            sample_count_sqrt,
            current_sample: 0,
            current_dimension: 0,
            permutation_key,
            rng,
            jitter_dist,
        }
    }

    /// Correlated multi-jittered sampling for sampling 2D points
    ///
    /// Works as follows:
    /// 1 - Sample points in canonical arrangements
    /// 2 - Permute "x" of each row using the same permutation
    /// 3 - Permute "y" of each column using the same permuttation
    /// 4 - Add jitter
    /// 5 - Return jittered & permuted x and y
    ///
    /// Since we use full (sample_count * sample_count) grid, the variance
    /// resulting from jitter is reduced compared to using a
    /// (sample_count_sqrt * sample_count_sqrt) grid
    fn cmj(&mut self) -> (Float, Float) {
        let cmj_ij = |index: u32| -> (u32, u32) {
            (
                index / self.sample_count_sqrt,
                index % self.sample_count_sqrt,
            )
        };
        let cmj_xy = |i: u32, j: u32| -> (u32, u32) {
            // returns (x, y) of the "canonical arrangement"
            // (x, y) in [0, self.sample_count]^2
            (
                (i * self.sample_count_sqrt + j),
                (j * self.sample_count_sqrt + i),
            )
        };

        let (i, j) = cmj_ij(self.current_sample);
        let i_permuted = permute_index(
            i,
            self.sample_count_sqrt,
            self.permutation_key.wrapping_add(self.current_dimension),
        );
        let j_permuted = permute_index(
            j,
            self.sample_count_sqrt,
            self.permutation_key
                .wrapping_add(self.current_dimension + 1),
        );
        let (x, y) = (cmj_xy(i_permuted, j).0, cmj_xy(i, j_permuted).1);

        let jitter_x = self.jitter_dist.sample(&mut self.rng);
        let jitter_y = self.jitter_dist.sample(&mut self.rng);
        self.current_dimension += 2;
        (
            (x as Float + jitter_x) / self.sample_count as Float,
            (y as Float + jitter_y) / self.sample_count as Float,
        )
    }

    /// Samples 2D points by dividing the space to (sample_count_sqrt * sample_count_sqrt) grid,
    /// choosing a stratum randomly, and adding jitter
    #[allow(dead_code)]
    fn random_stratum_sampling(&mut self) -> (Float, Float) {
        let stratum = permute_index(
            self.current_sample,
            self.sample_count,
            self.permutation_key.wrapping_add(self.current_dimension),
        );
        let x = stratum / self.sample_count_sqrt;
        let y = stratum % self.sample_count_sqrt;
        let jitter_x = self.jitter_dist.sample(&mut self.rng);
        let jitter_y = self.jitter_dist.sample(&mut self.rng);

        self.current_dimension += 2;
        (
            (x as Float + jitter_x) / self.sample_count_sqrt as Float,
            (y as Float + jitter_y) / self.sample_count_sqrt as Float,
        )
    }
}

impl Sampler for StratifiedSampler {
    fn get_1d(&mut self) -> Float {
        let stratum = permute_index(
            self.current_sample,
            self.sample_count,
            self.permutation_key.wrapping_add(self.current_dimension),
        );
        let jitter = self.jitter_dist.sample(&mut self.rng);

        self.current_dimension += 1;
        (stratum as Float + jitter) / self.sample_count as Float
    }

    fn get_2d(&mut self) -> (Float, Float) {
        self.cmj()
    }
    fn advance_sample(&mut self) {
        self.current_dimension = 0;
        self.current_sample += 1;
    }
    fn advance_pixel(&mut self) {
        self.current_dimension = 0;
        self.current_sample = 0;
        self.permutation_key = self.rng.next_u32();
    }
}
