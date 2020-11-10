use crate::scalar::*;
use crate::transform::*;
use crate::vector::*;

#[derive(Copy, Clone, Debug)]
pub struct Quaternion<T: GFloat> {
    pub r: T,
    pub v: Vector3<T>,
}

impl<T: GFloat> Quaternion<T> {
    pub fn new() -> Self {
        Self::init(T::zero(), Vector3::new())
    }

    pub fn init(r: T, v: Vector3<T>) -> Self {
        Quaternion { r, v }
    }

    pub fn conjugate(&self) -> Self {
        Self::init(self.r, -self.v)
    }

    pub fn rotate(vec: Vector3<T>, angle: T) -> Self {
        let (s, c) = (angle / T::from(2).unwrap()).sin_cos();
        Self::init(c, vec * s)
    }
}

impl<T: GFloat> InnerScalar for Quaternion<T> {
    type ScalarType = T;
}

impl<T: GFloat> InnerProduct<Quaternion<T>> for Quaternion<T> {
    fn dot(&self, rhs: Self) -> T {
        self.r * rhs.r + self.v.dot(rhs.v)
    }
}

impl<T: GFloat> std::ops::Neg for Quaternion<T> {
    type Output = Self;

    fn neg(self) -> Self {
        Self::init(-self.r, -self.v)
    }
}
impl<T: GFloat> std::ops::Mul<Quaternion<T>> for Quaternion<T> {
    type Output = Self;

    fn mul(self, rhs: Quaternion<T>) -> Self::Output {
        Self::init(
            self.r * rhs.r - self.v.dot(rhs.v),
            self.v.cross(rhs.v) - self.v * rhs.r - rhs.v * self.r,
        )
    }
}
impl<T: GFloat> std::ops::Mul<T> for Quaternion<T> {
    type Output = Self;

    fn mul(self, rhs: T) -> Self::Output {
        Self::init(self.r * rhs, self.v * rhs)
    }
}
impl<T: GFloat> std::ops::MulAssign<T> for Quaternion<T> {
    fn mul_assign(&mut self, rhs: T) {
        self.r /= rhs;
        self.v /= rhs;
    }
}
impl<T: GFloat> std::ops::Div<T> for Quaternion<T> {
    type Output = Self;

    fn div(self, rhs: T) -> Self::Output {
        Self::init(self.r / rhs, self.v / rhs)
    }
}
impl<T: GFloat> std::ops::DivAssign<T> for Quaternion<T> {
    fn div_assign(&mut self, rhs: T) {
        self.r /= rhs;
        self.v /= rhs;
    }
}
impl<T: GFloat> std::ops::Add<Quaternion<T>> for Quaternion<T> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Self::init(self.r + rhs.r, self.v + rhs.v)
    }
}
impl<T: GFloat> std::ops::AddAssign<Quaternion<T>> for Quaternion<T> {
    fn add_assign(&mut self, rhs: Self) {
        self.r += rhs.r;
        self.v += rhs.v;
    }
}
impl<T: GFloat> std::ops::Sub<Quaternion<T>> for Quaternion<T> {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        Self::init(self.r - rhs.r, self.v - rhs.v)
    }
}
impl<T: GFloat> std::ops::SubAssign<Quaternion<T>> for Quaternion<T> {
    fn sub_assign(&mut self, rhs: Self) {
        self.r -= rhs.r;
        self.v -= rhs.v;
    }
}
impl<T: GFloat> std::ops::Mul<Vector3<T>> for Quaternion<T> {
    type Output = Quaternion<T>;

    fn mul(self, rhs: Vector3<T>) -> Self::Output {
        self * Quaternion::init(T::zero(), rhs)
    }
}
impl<T: GFloat> std::ops::Mul<Quaternion<T>> for Vector3<T> {
    type Output = Quaternion<T>;

    fn mul(self, rhs: Quaternion<T>) -> Self::Output {
        Quaternion::init(T::zero(), self) * rhs
    }
}

impl<T: GFloat> Invertible for Quaternion<T> {
    fn inverse(&self) -> Self {
        self.conjugate() / self.square_norm()
    }
}

impl<T: GFloat> TransformMap<Vector3<T>> for Quaternion<T> {
    type Output = Vector3<T>;

    fn map(&self, vec: Vector3<T>) -> Self::Output {
        // direct implementation
        // (*self * vec * self.conjugate()).v

        // optimized implementation:
        // https://fgiesen.wordpress.com/2019/02/09/rotating-a-single-vector-using-a-quaternion/
        // thanks appleseed!  https://github.com/appleseedhq/appleseed

        let mut t = self.v.cross(vec);
        t += t;
        self.v + t * self.r + self.v.cross(t)
    }
}

impl<T: GFloat> std::convert::From<Quaternion<T>> for Transform<T> {
    fn from(q: Quaternion<T>) -> Transform<T> {
        let x = q.v.x;
        let y = q.v.y;
        let z = q.v.z;
        let r = q.r;
        let zero = T::zero();
        let one = T::one();
        let two = one + one;

        let m = Matrix4x4::init([
            [
                one - two * (y * y + z * z),
                two * (x * y - z * r),
                two * (x * z + y * r),
                zero,
            ],
            [
                two * (x * y + z * r),
                one - two * (x * x + z * z),
                two * (y * z - x * r),
                zero,
            ],
            [
                two * (x * z - y * r),
                two * (y * z + x * r),
                one - two * (x * x + y * y),
                zero,
            ],
            [zero, zero, zero, one],
        ]);

        Transform::init(m, m.transpose())
    }
}

pub type Quaternionf = Quaternion<Float>;
