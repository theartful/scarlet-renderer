pub use crate::scalar::*;
pub use crate::vector::*;

#[derive(Copy, Clone, Debug)]
pub struct Ray<T: Scalar> {
    pub origin: Point3<T>,
    pub direction: Vector3<T>,
    pub tmax: T,
}

#[derive(Clone, Copy)]
pub struct RayInfo<T: GFloat> {
    pub reciprocal: Vector3<T>,
    pub sign: Vector3<u32>,
}

#[derive(Clone, Copy)]
pub struct RayWithInfo<T: GFloat> {
    pub ray: Ray<T>,
    pub info: RayInfo<T>,
}

impl<T: Scalar> Ray<T> {
    pub fn new() -> Self {
        Ray::init(Point3::<T>::new(), Vector3::<T>::new(), T::highest())
    }

    pub fn init(origin: Point3<T>, direction: Vector3<T>, tmax: T) -> Self {
        Ray {
            origin,
            direction,
            tmax,
        }
    }

    pub fn eval(&self, t: T) -> Point3<T> {
        self.origin + self.direction * t
    }
}

impl<T: GFloat> RayInfo<T> {
    pub fn init(ray: Ray<T>) -> Self {
        RayInfo {
            reciprocal: Vector3::init(
                ray.direction.x.recip(),
                ray.direction.y.recip(),
                ray.direction.z.recip(),
            ),
            sign: Vector3::init(
                u32::sign(ray.direction.x),
                u32::sign(ray.direction.y),
                u32::sign(ray.direction.z),
            ),
        }
    }
}

impl<T: GFloat> RayWithInfo<T> {
    pub fn init(ray: Ray<T>) -> Self {
        RayWithInfo {
            ray,
            info: RayInfo::init(ray),
        }
    }
}

impl<T> From<Ray<T>> for Ray<Interval<T>>
where
    T: GFloat,
{
    fn from(ray: Ray<T>) -> Self {
        Self::init(
            Point3::from(ray.origin),
            Vector3::from(ray.direction),
            Interval::from(ray.tmax),
        )
    }
}
impl<T> From<Ray<Interval<T>>> for Ray<T>
where
    T: GFloat,
{
    fn from(ray: Ray<Interval<T>>) -> Self {
        Self::init(
            Point3::from(ray.origin),
            Vector3::from(ray.direction),
            ray.tmax.approx(),
        )
    }
}

pub type Rayf = Ray<Float>;
pub type Rayfi = Ray<Intervalf>;
pub type RayWithInfof = RayWithInfo<Float>;
