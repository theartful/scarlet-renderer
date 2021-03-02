pub use crate::ray::*;
pub use crate::transform::*;
pub use crate::vector::*;

pub struct CameraSample {
    // p_film in [-1, 1] x [-1/aspect_ratio, 1/aspect_ratio]
    p_film: Point2f,
    // p_lens in [-1, 1]^2
    p_lens: Point2f,
}

impl CameraSample {
    pub fn init(p_film: Point2f, p_lens: Point2f) -> Self {
        Self { p_film, p_lens }
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

pub struct ThinLensCamera {
    view: Transformf,
    projection: Transformf,
    focus_distance: Float,
    lens_radius: Float,
}

impl ThinLensCamera {
    pub fn init(
        view: Transformf,
        focus_distance: Float,
        fov: Float,
        lens_radius: Float,
        near: Float,
        far: Float,
    ) -> Self {
        let projection = Transformf::perspective(fov, near, far);
        Self {
            view,
            projection,
            focus_distance,
            lens_radius,
        }
    }
}

impl Camera for ThinLensCamera {
    fn generate_ray(&self, sample: CameraSample) -> Rayf {
        // all rays passing from a single point on the film converge to the same
        // point on the focal plane
        // we compute this point by constructing a ray from the film point
        // to the center of the lens and evaluate the ray at z=focus_distance
        // and construct a ray from p_lens to that point

        // map p_film (which is in NDC) to the near plane
        let p_film_camera =
            self.projection
                .map_inverse(Point3f::init(sample.p_film.x, sample.p_film.y, 0.));
        // construct ray from center of the lens to the point on the film
        let ray = Rayf::init(Point3f::init(0., 0., 0.), Vector3f::from(p_film_camera), 0.);

        // find t at z=focus_distance
        let t = self.focus_distance / ray.direction.z;

        // evaluate ray intersection at the focal plane
        let p_focus_plane = ray.eval(t);

        let lens_xy = sample_concentric_disk(sample.p_lens) * self.lens_radius;
        let p_lens_camera = Point3f::init(lens_xy.x, lens_xy.y, 0.);

        self.view.map_inverse(Rayf::init(
            p_lens_camera,
            (p_focus_plane - p_lens_camera).normalize(),
            Float::highest(),
        ))
    }
}

// TODO
pub struct OrthographicCamera {
    view: Transformf,
    focus_distance: Float,
    lens_radius: Float,
}
