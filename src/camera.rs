pub use crate::ray::*;
pub use crate::transform::*;
pub use crate::vector::*;

pub struct CameraSample {
    // p_film in [-1, 1] x [-1/aspect_ratio, 1/aspect_ratio]
    p_film: Point2f,
    // p_aperture in [-1, 1]^2
    p_aperture: Point2f,
}

impl CameraSample {
    pub fn init(p_film: Point2f, p_aperture: Point2f) -> Self {
        Self { p_film, p_aperture }
    }
}

pub trait Camera {
    fn generate_ray(&self, sample: CameraSample) -> Rayf;
}

// http://psgraphics.blogspot.com/2011/01/improved-code-for-concentric-map.html
// http://l2program.co.uk/900/concentric-disk-sampling
// maps unit square to unit circle
pub fn sample_concentric_disk(p: Point2f) -> Point2f {
    if p.x == 0. && p.y == 0. {
        return Point2f::init(0., 0.);
    }
    let (r, theta) = if p.x.abs() > p.y.abs() {
        (p.x, Float::PI * 0.25 * p.y / p.x)
    } else {
        (p.y, Float::PI * 0.5 - (p.x / p.y) * Float::PI * 0.25)
    };

    Point2f::init(r * theta.cos(), r * theta.sin())
}

pub struct PerspectiveCamera {
    view: Transformf,
    focal_length: Float,
    fov_tan: Float,
    aperture_radius: Float,
}

impl PerspectiveCamera {
    pub fn init(view: Transformf, focal_length: Float, fov: Float, aperture_radius: Float) -> Self {
        Self {
            view,
            focal_length,
            fov_tan: (fov / 2.).tan(),
            aperture_radius,
        }
    }

    pub fn camera_center(&self) -> Point3f {
        Point3f::from(Vector3f::from(self.view.matrix().col_vec(3)))
    }
}

impl Camera for PerspectiveCamera {
    fn generate_ray(&self, sample: CameraSample) -> Rayf {
        // convert p_film to camera coordinates
        let dir = self
            .view
            .map_inverse(Vector3f::init(
                self.focal_length * self.fov_tan * sample.p_film.x,
                self.focal_length * self.fov_tan * sample.p_film.y,
                self.focal_length,
            ))
            .normalize();
        let origin_2 = sample_concentric_disk(sample.p_aperture) * self.aperture_radius;
        let origin = Point3f::init(origin_2.x, origin_2.y, 0.);
        Rayf::init(self.view.map_inverse(origin), dir, Float::highest())
    }
}

// TODO
pub struct OrthographicCamera {
    view: Transformf,
    focal_length: Float,
    aperture_radius: Float,
}
