use crate::scalar::*;
use crate::vector_traits::*;

pub fn lerp<U, T>(a: U, b: U, t: T) -> U
where
    T: Scalar,
    U: std::ops::Mul<T, Output = U> + std::ops::Add<U, Output = U>,
{
    a * (T::one() - t) + b * t
}

pub fn slerp<T, U>(a: U, b: U, t: T) -> U
where
    T: GFloat,
    U: Norm + VectorTraits<ScalarType = T>,
{
    let cos_theta = a.dot(b);
    if cos_theta > T::from(0.98).unwrap() {
        lerp(a, b, t).normalize()
    } else {
        let theta = clamp(cos_theta, -T::one(), T::one()).acos();
        let t_theta = t * theta;
        let perp = (a - b * cos_theta).normalize();
        a * t_theta.cos() + perp * t_theta.sin()
    }
}
