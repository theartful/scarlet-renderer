pub use crate::primitive::*;

// A sphere is parameterized by the equations
// x = r * sin theta cos phi
// y = r * sin theta sin phi
// z = r * cos theta
// where phi in [0, 2 pi] and theta in [0, pi]
//
// this means that
// phi = atan2(y, x)
// theta = acos(z)

pub struct Sphere {
    radius: Float,
    z_min: Float,
    z_max: Float,
    theta_min: Float,
    theta_max: Float,
    phi_max: Float,
}

struct IntersectionParams {
    t: Intervalf,
    // intersection point in object space
    p: Point3fi,
    phi: Float,
}

impl Bounded3<Float> for Sphere {
    fn bbox(&self) -> Bbox3f {
        Bbox3f::init(
            Point3f::init(-self.radius, -self.radius, self.z_min),
            Point3f::init(self.radius, self.radius, self.z_max),
        )
    }
}

impl Sphere {
    pub fn init(radius: Float, z_min: Float, z_max: Float, phi_max: Float) -> Self {
        assert!(z_min < z_max);

        let theta_min = clamp(z_min / radius, Float::zero(), Float::one()).acos();
        let theta_max = clamp(z_max / radius, Float::zero(), Float::one()).acos();
        Self {
            radius,
            z_min,
            z_max,
            theta_min,
            theta_max,
            phi_max,
        }
    }
}

// intersection related
impl Sphere {
    fn get_intersection_params(&self, ray: Rayfi) -> Option<IntersectionParams> {
        // solve quadratic equation a t^2 + b t + c = 0 which is the result of
        // substituting parameteric line equation in the sphere equation
        // with a = d_x^2 + d_y^2 + d_z^2 ~= 1
        //      b = 2(o_x d_x + o_y d_y + o_z d_z)
        //      c = o_x^2 + o_y^2 + o_z^2 - r^2

        let a = ray.direction.square_norm();
        let b = Vector3fi::from(ray.origin).dot(ray.direction) * 2.;
        let c = Vector3fi::from(ray.origin).square_norm() - self.radius * self.radius;

        let (t0, t1) = solve_quadratic(a, b, c)?;

        let get_intersection = |t: Intervalf| -> Option<IntersectionParams> {
            if !t.in_range(Float::zero(), ray.tmax.inf) {
                return None;
            }

            // calculate hit point
            let hit_point = ray.origin + ray.direction * t;
            // should I refine hitpoint by projecting it on the sphere instead of using interval
            // arithmetic?

            if !hit_point.z.in_range(self.z_min, self.z_max) {
                return None;
            }

            let phi = {
                // atan2 returns in the range [-pi, pi]
                let tmp = hit_point.y.approx().atan2(hit_point.x.approx());
                if tmp.is_sign_negative() {
                    tmp + 2. * Float::PI
                } else {
                    tmp
                }
            };

            if phi > self.phi_max {
                return None;
            }

            Some(IntersectionParams {
                t: t,
                p: hit_point,
                phi: phi,
            })
        };

        match get_intersection(t0) {
            None => get_intersection(t1),
            Some(intersection) => Some(intersection),
        }
    }
}

impl Primitive for Sphere {
    fn intersecti(&self, ray: &mut Rayfi) -> Option<SurfaceInteraction> {
        let intersection_params = self.get_intersection_params(*ray)?;

        let t = intersection_params.t;
        let p = Point3f::from(intersection_params.p);
        let phi = intersection_params.phi;
        let theta = clamp(p.z / self.radius, Float::zero(), Float::one()).acos();

        let sin_theta = theta.sin();
        let (sin_phi, cos_phi) = phi.sin_cos();

        // (u, v) in [0, 1]^2 while (phi, theta) in [0, 2 pi] x [0, pi]
        // so we transform them
        let u = phi / self.phi_max;
        let v = (theta - self.theta_min) / (self.theta_max - self.theta_min);

        // now we calculate dp/du and dp/dv
        // phi   = u * phi_max
        // theta = v * theta_max + theta_min
        //
        // taking the derivatives we get

        // dxdu = -phi_max * radius * sin(theta) * sin(phi);
        // dydu = phi_max * radius * sin(theta) * cos(phi);
        // dzdu = 0

        // dxdv = (theta_max - theta_min) * radius * cos(theta) * cos(phi);
        // dydv = (theta_max - theta_min) * radius * cos(theta) * sin(phi);
        // dzdv = -(theta_max - theta_min) * radius * sin(theta);

        // which is
        let dpdu = Vector3f::init(-p.y, p.x, Float::zero()) * self.phi_max;
        let dpdv = Vector3f::init(p.z * cos_phi, p.z * sin_phi, -self.radius * sin_theta)
            * (self.theta_max - self.theta_min);

        ray.tmax = t;
        Some(SurfaceInteraction::init(
            intersection_params.p,
            Point2f::init(u, v),
            dpdu,
            dpdv,
        ))
    }
    fn area(&self) -> Float {
        // TODO: derive this formula
        self.phi_max * self.radius * (self.z_max - self.z_min)
    }
}
