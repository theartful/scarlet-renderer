pub use crate::scalar::*;
pub use crate::vector_traits::*;

#[derive(Debug)]
pub struct VectorMarker {}
#[derive(Debug)]
pub struct PointMarker {}
#[derive(Debug)]
pub struct NormalMarker {}

pub trait HasSub {}
impl HasSub for VectorMarker {}

pub type Vector4<T> = GenericVector4<T, VectorMarker>;
pub type Vector3<T> = GenericVector3<T, VectorMarker>;
pub type Vector2<T> = GenericVector2<T, VectorMarker>;

pub type Point3<T> = GenericVector3<T, PointMarker>;
pub type Point2<T> = GenericVector2<T, PointMarker>;

pub type Normal3<T> = GenericVector3<T, NormalMarker>;

#[derive(Debug, PartialEq)]
pub struct GenericVector4<T: Scalar, U> {
    pub x: T,
    pub y: T,
    pub z: T,
    pub w: T,
    p: std::marker::PhantomData<U>,
}

#[derive(Debug, PartialEq)]
pub struct GenericVector3<T: Scalar, U> {
    pub x: T,
    pub y: T,
    pub z: T,
    p: std::marker::PhantomData<U>,
}

#[derive(Debug, PartialEq)]
pub struct GenericVector2<T: Scalar, U> {
    pub x: T,
    pub y: T,
    p: std::marker::PhantomData<U>,
}

impl<T: Scalar, U> InnerScalar for GenericVector4<T, U> {
    type ScalarType = T;
}
impl<T: Scalar, U> InnerScalar for GenericVector3<T, U> {
    type ScalarType = T;
}
impl<T: Scalar, U> InnerScalar for GenericVector2<T, U> {
    type ScalarType = T;
}

// Operator overloading
// Addition and subtraction
// GenericVector4
impl<T: SignedScalar, U> std::ops::Neg for GenericVector4<T, U> {
    type Output = GenericVector4<T, U>;

    fn neg(self) -> Self::Output {
        GenericVector4::init(-self.x, -self.y, -self.z, -self.w)
    }
}

impl<T: Scalar, U> std::ops::Add<GenericVector4<T, U>> for GenericVector4<T, U> {
    type Output = GenericVector4<T, U>;

    fn add(self, rhs: GenericVector4<T, U>) -> Self::Output {
        GenericVector4::init(
            self.x + rhs.x,
            self.y + rhs.y,
            self.z + rhs.z,
            self.w + rhs.w,
        )
    }
}

impl<T: Scalar, U> std::ops::AddAssign<GenericVector4<T, U>> for GenericVector4<T, U> {
    fn add_assign(&mut self, rhs: GenericVector4<T, U>) {
        self.x += rhs.x;
        self.y += rhs.y;
        self.z += rhs.z;
        self.w += rhs.w;
    }
}

impl<T: Scalar, U: HasSub> std::ops::Sub<GenericVector4<T, U>> for GenericVector4<T, U> {
    type Output = GenericVector4<T, U>;

    fn sub(self, rhs: GenericVector4<T, U>) -> Self::Output {
        GenericVector4::init(
            self.x - rhs.x,
            self.y - rhs.y,
            self.z - rhs.z,
            self.w - rhs.w,
        )
    }
}

impl<T: Scalar, U: HasSub> std::ops::SubAssign<GenericVector4<T, U>> for GenericVector4<T, U> {
    fn sub_assign(&mut self, rhs: GenericVector4<T, U>) {
        self.x -= rhs.x;
        self.y -= rhs.y;
        self.z -= rhs.z;
        self.w -= rhs.w;
    }
}

// GenericVector3
impl<T: SignedScalar, U> std::ops::Neg for GenericVector3<T, U> {
    type Output = GenericVector3<T, U>;

    fn neg(self) -> Self::Output {
        GenericVector3::init(-self.x, -self.y, -self.z)
    }
}

impl<T: Scalar, U> std::ops::Add<GenericVector3<T, U>> for GenericVector3<T, U> {
    type Output = GenericVector3<T, U>;

    fn add(self, rhs: GenericVector3<T, U>) -> Self::Output {
        GenericVector3::init(self.x + rhs.x, self.y + rhs.y, self.z + rhs.z)
    }
}

impl<T: Scalar, U> std::ops::AddAssign<GenericVector3<T, U>> for GenericVector3<T, U> {
    fn add_assign(&mut self, rhs: GenericVector3<T, U>) {
        self.x += rhs.x;
        self.y += rhs.y;
        self.z += rhs.z;
    }
}

impl<T: Scalar, U: HasSub> std::ops::Sub<GenericVector3<T, U>> for GenericVector3<T, U> {
    type Output = GenericVector3<T, U>;

    fn sub(self, rhs: GenericVector3<T, U>) -> Self::Output {
        GenericVector3::init(self.x - rhs.x, self.y - rhs.y, self.z - rhs.z)
    }
}

impl<T: Scalar, U: HasSub> std::ops::SubAssign<GenericVector3<T, U>> for GenericVector3<T, U> {
    fn sub_assign(&mut self, rhs: GenericVector3<T, U>) {
        self.x -= rhs.x;
        self.y -= rhs.y;
        self.z -= rhs.z;
    }
}

// GenericVector2
impl<T: SignedScalar, U> std::ops::Neg for GenericVector2<T, U> {
    type Output = GenericVector2<T, U>;

    fn neg(self) -> Self::Output {
        GenericVector2::init(-self.x, -self.y)
    }
}

impl<T: Scalar, U> std::ops::Add<GenericVector2<T, U>> for GenericVector2<T, U> {
    type Output = GenericVector2<T, U>;

    fn add(self, rhs: GenericVector2<T, U>) -> Self::Output {
        GenericVector2::init(self.x + rhs.x, self.y + rhs.y)
    }
}

impl<T: Scalar, U> std::ops::AddAssign<GenericVector2<T, U>> for GenericVector2<T, U> {
    fn add_assign(&mut self, rhs: GenericVector2<T, U>) {
        self.x += rhs.x;
        self.y += rhs.y;
    }
}

impl<T: Scalar, U: HasSub> std::ops::Sub<GenericVector2<T, U>> for GenericVector2<T, U> {
    type Output = GenericVector2<T, U>;

    fn sub(self, rhs: GenericVector2<T, U>) -> Self::Output {
        GenericVector2::init(self.x - rhs.x, self.y - rhs.y)
    }
}

impl<T: Scalar, U: HasSub> std::ops::SubAssign<GenericVector2<T, U>> for GenericVector2<T, U> {
    fn sub_assign(&mut self, rhs: GenericVector2<T, U>) {
        self.x -= rhs.x;
        self.y -= rhs.y;
    }
}

// Scaling
// GenericVector4
impl<T: Scalar, U> std::ops::Mul<T> for GenericVector4<T, U> {
    type Output = GenericVector4<T, U>;

    fn mul(self, rhs: T) -> Self::Output {
        GenericVector4::init(self.x * rhs, self.y * rhs, self.z * rhs, self.w * rhs)
    }
}

impl<T: Scalar, U> std::ops::MulAssign<T> for GenericVector4<T, U> {
    fn mul_assign(&mut self, rhs: T) {
        self.x *= rhs;
        self.y *= rhs;
        self.z *= rhs;
        self.w *= rhs;
    }
}

impl<T: Scalar, U> std::ops::Div<T> for GenericVector4<T, U> {
    type Output = GenericVector4<T, U>;

    fn div(self, rhs: T) -> Self::Output {
        GenericVector4::init(self.x / rhs, self.y / rhs, self.z / rhs, self.w / rhs)
    }
}

impl<T: Scalar, U> std::ops::DivAssign<T> for GenericVector4<T, U> {
    fn div_assign(&mut self, rhs: T) {
        self.x /= rhs;
        self.y /= rhs;
        self.z /= rhs;
        self.w /= rhs;
    }
}

// GenericVector3
impl<T: Scalar, U> std::ops::Mul<T> for GenericVector3<T, U> {
    type Output = GenericVector3<T, U>;

    fn mul(self, rhs: T) -> Self::Output {
        GenericVector3::init(self.x * rhs, self.y * rhs, self.z * rhs)
    }
}

impl<T: Scalar, U> std::ops::MulAssign<T> for GenericVector3<T, U> {
    fn mul_assign(&mut self, rhs: T) {
        self.x *= rhs;
        self.y *= rhs;
        self.z *= rhs;
    }
}

impl<T: Scalar, U> std::ops::Div<T> for GenericVector3<T, U> {
    type Output = GenericVector3<T, U>;

    fn div(self, rhs: T) -> Self::Output {
        GenericVector3::init(self.x / rhs, self.y / rhs, self.z / rhs)
    }
}

impl<T: Scalar, U> std::ops::DivAssign<T> for GenericVector3<T, U> {
    fn div_assign(&mut self, rhs: T) {
        self.x /= rhs;
        self.y /= rhs;
        self.z /= rhs;
    }
}

// GenericVector2
impl<T: Scalar, U> std::ops::Mul<T> for GenericVector2<T, U> {
    type Output = GenericVector2<T, U>;

    fn mul(self, rhs: T) -> Self::Output {
        GenericVector2::init(self.x * rhs, self.y * rhs)
    }
}

impl<T: Scalar, U> std::ops::MulAssign<T> for GenericVector2<T, U> {
    fn mul_assign(&mut self, rhs: T) {
        self.x *= rhs;
        self.y *= rhs;
    }
}

impl<T: Scalar, U> std::ops::Div<T> for GenericVector2<T, U> {
    type Output = GenericVector2<T, U>;

    fn div(self, rhs: T) -> Self::Output {
        GenericVector2::init(self.x / rhs, self.y / rhs)
    }
}

impl<T: Scalar, U> std::ops::DivAssign<T> for GenericVector2<T, U> {
    fn div_assign(&mut self, rhs: T) {
        self.x /= rhs;
        self.y /= rhs;
    }
}

impl<T: SignedScalar, U: HasSub, I> InnerProduct<GenericVector4<T, I>> for GenericVector4<T, U> {
    fn dot(&self, rhs: GenericVector4<T, I>) -> Self::ScalarType {
        self.x * rhs.x + self.y * rhs.y + self.z * rhs.z + self.w * rhs.w
    }
}
impl<T: SignedScalar, U: HasSub, I> InnerProduct<GenericVector3<T, I>> for GenericVector3<T, U> {
    fn dot(&self, rhs: GenericVector3<T, I>) -> Self::ScalarType {
        self.x * rhs.x + self.y * rhs.y + self.z * rhs.z
    }
}
impl<T: SignedScalar, U: HasSub, I> InnerProduct<GenericVector2<T, I>> for GenericVector2<T, U> {
    fn dot(&self, rhs: GenericVector2<T, I>) -> Self::ScalarType {
        self.x * rhs.x + self.y * rhs.y
    }
}

// main traits
impl<T: Scalar, U> GenericVector4<T, U> {
    pub fn new() -> Self {
        GenericVector4::init(T::zero(), T::zero(), T::zero(), T::zero())
    }

    pub fn init(x: T, y: T, z: T, w: T) -> Self {
        let p = std::marker::PhantomData;
        GenericVector4 { x, y, z, w, p }
    }

    pub fn get(&self, axis: u8) -> T {
        if axis == 0 {
            self.x
        } else if axis == 1 {
            self.y
        } else if axis == 2 {
            self.z
        } else if axis == 3 {
            self.w
        } else {
            panic!();
        }
    }
}

impl<T: Scalar, U> GenericVector3<T, U> {
    pub fn new() -> Self {
        GenericVector3::init(T::zero(), T::zero(), T::zero())
    }

    pub fn init(x: T, y: T, z: T) -> Self {
        let p = std::marker::PhantomData;
        GenericVector3 { x, y, z, p }
    }

    pub fn cross(&self, vec: GenericVector3<T, U>) -> Self {
        GenericVector3::init(
            self.y * vec.z - self.z * vec.y,
            self.z * vec.x - self.x * vec.z,
            self.x * vec.y - self.y * vec.x,
        )
    }

    pub fn get(&self, axis: u8) -> T {
        if axis == 0 {
            self.x
        } else if axis == 1 {
            self.y
        } else if axis == 2 {
            self.z
        } else {
            panic!();
        }
    }
}

impl<T: Scalar, U> GenericVector2<T, U> {
    pub fn new() -> Self {
        GenericVector2::init(T::zero(), T::zero())
    }

    pub fn init(x: T, y: T) -> Self {
        let p = std::marker::PhantomData;
        GenericVector2 { x, y, p }
    }

    pub fn get(&self, axis: u8) -> T {
        if axis == 0 {
            self.x
        } else if axis == 1 {
            self.y
        } else {
            panic!();
        }
    }
}

impl<T: GFloat, U> From<GenericVector4<T, U>> for GenericVector3<T, U> {
    fn from(vec: GenericVector4<T, U>) -> GenericVector3<T, U> {
        // TODO: remove comment after fequals is implemented
        // std::assert!(fequals(vec.w, T::zero()));
        GenericVector3::init(vec.x, vec.y, vec.z)
    }
}

// interactions between vector/point/normal
impl<T: Scalar> std::ops::Sub<Point3<T>> for Point3<T> {
    type Output = Vector3<T>;

    fn sub(self, rhs: Point3<T>) -> Self::Output {
        Vector3::init(self.x - rhs.x, self.y - rhs.y, self.z - rhs.z)
    }
}
impl<T: Scalar> std::ops::Sub<Point2<T>> for Point2<T> {
    type Output = Vector2<T>;

    fn sub(self, rhs: Point2<T>) -> Self::Output {
        Vector2::init(self.x - rhs.x, self.y - rhs.y)
    }
}
impl<T: Scalar> std::ops::Add<Point3<T>> for Vector3<T> {
    type Output = Point3<T>;

    fn add(self, rhs: Point3<T>) -> Self::Output {
        Point3::init(self.x + rhs.x, self.y + rhs.y, self.z + rhs.z)
    }
}
impl<T: Scalar> std::ops::Add<Vector3<T>> for Point3<T> {
    type Output = Point3<T>;

    fn add(self, rhs: Vector3<T>) -> Self::Output {
        Point3::init(self.x + rhs.x, self.y + rhs.y, self.z + rhs.z)
    }
}
impl<T: Scalar> std::ops::Sub<Vector3<T>> for Point3<T> {
    type Output = Point3<T>;

    fn sub(self, rhs: Vector3<T>) -> Self::Output {
        Point3::init(self.x - rhs.x, self.y - rhs.y, self.z - rhs.z)
    }
}
impl<T: Scalar> From<Vector3<T>> for Normal3<T> {
    fn from(p: Vector3<T>) -> Normal3<T> {
        Normal3::init(p.x, p.y, p.z)
    }
}
impl<T: Scalar> From<Point3<T>> for Vector3<T> {
    fn from(p: Point3<T>) -> Vector3<T> {
        Vector3::init(p.x, p.y, p.z)
    }
}
impl<T: Scalar> From<Vector3<T>> for Point3<T> {
    fn from(p: Vector3<T>) -> Point3<T> {
        Point3::init(p.x, p.y, p.z)
    }
}
impl<T: Scalar> From<Normal3<T>> for Vector3<T> {
    fn from(p: Normal3<T>) -> Vector3<T> {
        Vector3::init(p.x, p.y, p.z)
    }
}

// useless stuff
impl<T: Scalar, U> Copy for GenericVector4<T, U> {}
impl<T: Scalar, U> Clone for GenericVector4<T, U> {
    fn clone(&self) -> Self {
        *self
    }
}
impl<T: Scalar, U> Copy for GenericVector3<T, U> {}
impl<T: Scalar, U> Clone for GenericVector3<T, U> {
    fn clone(&self) -> Self {
        *self
    }
}
impl<T: Scalar, U> Copy for GenericVector2<T, U> {}
impl<T: Scalar, U> Clone for GenericVector2<T, U> {
    fn clone(&self) -> Self {
        *self
    }
}

pub type Vector4f = Vector4<Float>;
pub type Vector4i = Vector4<Int>;
pub type Vector4u = Vector4<UInt>;

pub type Vector3f = Vector3<Float>;
pub type Vector3i = Vector3<Int>;
pub type Vector3u = Vector3<UInt>;

pub type Vector2f = Vector2<Float>;
pub type Vector2i = Vector2<Int>;
pub type Vector2u = Vector2<UInt>;

pub type Point3f = Point3<Float>;
pub type Point3i = Point3<Int>;
pub type Point3u = Point3<UInt>;

pub type Point2f = Point2<Float>;
pub type Point2i = Point2<Int>;
pub type Point2u = Point2<UInt>;

pub type Normal3f = Normal3<Float>;
pub type Normal3i = Normal3<Int>;

// is this the place for these to be here?
pub use crate::interval::*;
pub type Vector3fi = Vector3<Interval<Float>>;
pub type Vector3ii = Vector3<Interval<Int>>;
pub type Point3fi = Point3<Interval<Float>>;
pub type Point3ii = Point3<Interval<Int>>;

impl<T: GFloat, U> std::ops::Mul<T> for GenericVector3<Interval<T>, U> {
    type Output = GenericVector3<Interval<T>, U>;

    fn mul(self, rhs: T) -> Self::Output {
        GenericVector3::init(self.x * rhs, self.y * rhs, self.z * rhs)
    }
}
impl<T: GFloat, U> std::ops::Div<T> for GenericVector3<Interval<T>, U> {
    type Output = GenericVector3<Interval<T>, U>;

    fn div(self, rhs: T) -> Self::Output {
        GenericVector3::init(self.x / rhs, self.y / rhs, self.z / rhs)
    }
}

impl<T: GFloat, U> From<GenericVector3<T, U>> for GenericVector3<Interval<T>, U> {
    fn from(vec: GenericVector3<T, U>) -> GenericVector3<Interval<T>, U> {
        Self::init(
            Interval::from(vec.x),
            Interval::from(vec.y),
            Interval::from(vec.z),
        )
    }
}
impl<T: GFloat, U> From<GenericVector3<Interval<T>, U>> for GenericVector3<T, U> {
    fn from(vec: GenericVector3<Interval<T>, U>) -> GenericVector3<T, U> {
        Self::init(vec.x.approx(), vec.y.approx(), vec.z.approx())
    }
}

// this is more accurate than self.dot(self) because of interval arithmetic
impl Norm for Vector3fi {
    fn square_norm(self) -> Self::ScalarType {
        self.x.square() + self.y.square() + self.z.square()
    }
}
