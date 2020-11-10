use crate::bbox::*;
use crate::ray::*;
use crate::transform_impl::*;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Matrix4x4<T: Scalar> {
    pub data: [[T; 4]; 4],
}

impl<T: Scalar> Matrix4x4<T> {
    pub fn new() -> Self {
        Self::identity()
    }

    pub fn init(data: [[T; 4]; 4]) -> Self {
        Matrix4x4 { data }
    }

    pub fn identity() -> Self {
        Matrix4x4 {
            data: [
                [T::one(), T::zero(), T::zero(), T::zero()],
                [T::zero(), T::one(), T::zero(), T::zero()],
                [T::zero(), T::zero(), T::one(), T::zero()],
                [T::zero(), T::zero(), T::zero(), T::one()],
            ],
        }
    }

    pub fn zero() -> Self {
        Matrix4x4 {
            data: [
                [T::zero(), T::zero(), T::zero(), T::zero()],
                [T::zero(), T::zero(), T::zero(), T::zero()],
                [T::zero(), T::zero(), T::zero(), T::zero()],
                [T::zero(), T::zero(), T::zero(), T::zero()],
            ],
        }
    }

    pub fn transpose(&self) -> Matrix4x4<T> {
        let d = &self.data;
        Matrix4x4 {
            data: [
                [d[0][0], d[1][0], d[2][0], d[3][0]],
                [d[0][1], d[1][1], d[2][1], d[3][1]],
                [d[0][2], d[1][2], d[2][2], d[3][2]],
                [d[0][3], d[1][3], d[2][3], d[3][3]],
            ],
        }
    }

    pub fn mul_scalar(&self, rhs: T) -> Self {
        let d = &self.data;
        Matrix4x4 {
            data: [
                [rhs * d[0][0], rhs * d[0][1], rhs * d[0][2], rhs * d[0][3]],
                [rhs * d[1][0], rhs * d[1][1], rhs * d[1][2], rhs * d[1][3]],
                [rhs * d[2][0], rhs * d[2][1], rhs * d[2][2], rhs * d[2][3]],
                [rhs * d[3][0], rhs * d[3][1], rhs * d[3][2], rhs * d[3][3]],
            ],
        }
    }

    pub fn mul(&self, rhs: &Self) -> Self {
        let mut res = Self::zero();
        for i in 0..4 {
            for j in 0..4 {
                for k in 0..4 {
                    res.data[i][j] += self.data[i][k] * rhs.data[k][j];
                }
            }
        }
        res
    }

    pub fn row(&self, r: usize) -> [T; 4] {
        self.data[r]
    }

    pub fn col(&self, c: usize) -> [T; 4] {
        [
            self.data[0][c],
            self.data[1][c],
            self.data[2][c],
            self.data[3][c],
        ]
    }
    pub fn row_vec(&self, r: usize) -> Vector4<T> {
        Vector4::init(
            self.data[r][0],
            self.data[r][1],
            self.data[r][2],
            self.data[r][3],
        )
    }

    pub fn col_vec(&self, c: usize) -> Vector4<T> {
        Vector4::init(
            self.data[0][c],
            self.data[1][c],
            self.data[2][c],
            self.data[3][c],
        )
    }
}

impl<T: Scalar> std::ops::Index<(usize, usize)> for Matrix4x4<T> {
    type Output = T;

    fn index(&self, idx: (usize, usize)) -> &Self::Output {
        &self.data[idx.0][idx.1]
    }
}

impl<T: Scalar> std::ops::Index<usize> for Matrix4x4<T> {
    type Output = [T; 4];

    fn index(&self, idx: usize) -> &Self::Output {
        &self.data[idx]
    }
}

impl<T: Scalar> std::ops::IndexMut<(usize, usize)> for Matrix4x4<T> {
    fn index_mut(&mut self, idx: (usize, usize)) -> &mut Self::Output {
        &mut self.data[idx.0][idx.1]
    }
}

impl<T: Scalar> std::ops::IndexMut<usize> for Matrix4x4<T> {
    fn index_mut(&mut self, idx: usize) -> &mut Self::Output {
        &mut self.data[idx]
    }
}

impl<T: Scalar> std::ops::Mul<T> for &Matrix4x4<T> {
    type Output = Matrix4x4<T>;
    fn mul(self, rhs: T) -> Self::Output {
        self.mul_scalar(rhs)
    }
}

impl<'a, T: Scalar> std::ops::Mul<&'a Matrix4x4<T>> for &Matrix4x4<T> {
    type Output = Matrix4x4<T>;
    fn mul(self, rhs: &'a Matrix4x4<T>) -> Self::Output {
        self.mul(rhs)
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Transform<T: GFloat> {
    m: Matrix4x4<T>,
    minv: Matrix4x4<T>,
}

impl<T: GFloat> Transform<T> {
    pub fn new() -> Self {
        Transform {
            m: Matrix4x4::<T>::identity(),
            minv: Matrix4x4::<T>::identity(),
        }
    }

    pub fn init(m: Matrix4x4<T>, minv: Matrix4x4<T>) -> Self {
        Transform { m, minv }
    }

    fn init_from_orthogonal_basis(
        camera_pos: Point3<T>,
        x: Vector3<T>,
        y: Vector3<T>,
        z: Vector3<T>,
    ) -> Self {
        Transform {
            m: Matrix4x4::<T>::init([
                [x.x, x.y, x.z, -x.dot(camera_pos)],
                [y.x, y.y, y.z, -y.dot(camera_pos)],
                [z.x, z.y, z.z, -z.dot(camera_pos)],
                [T::zero(), T::zero(), T::zero(), T::one()],
            ]),
            minv: Matrix4x4::<T>::init([
                [x.x, y.x, z.x, camera_pos.x],
                [x.y, y.y, z.y, camera_pos.y],
                [x.z, y.z, z.z, camera_pos.z],
                [T::zero(), T::zero(), T::zero(), T::one()],
            ]),
        }
    }

    pub fn translate(vec: Vector3<T>) -> Self {
        Self::init(
            Matrix4x4::init([
                [T::one(), T::zero(), T::zero(), vec.x],
                [T::zero(), T::one(), T::zero(), vec.y],
                [T::zero(), T::zero(), T::one(), vec.z],
                [T::zero(), T::zero(), T::zero(), T::one()],
            ]),
            Matrix4x4::init([
                [T::one(), T::zero(), T::zero(), -vec.x],
                [T::zero(), T::one(), T::zero(), -vec.y],
                [T::zero(), T::zero(), T::one(), -vec.z],
                [T::zero(), T::zero(), T::zero(), T::one()],
            ]),
        )
    }

    pub fn scale(vec: Vector3<T>) -> Self {
        Self::init(
            Matrix4x4::init([
                [vec.x, T::zero(), T::zero(), T::zero()],
                [T::zero(), vec.y, T::zero(), T::zero()],
                [T::zero(), T::zero(), vec.z, T::zero()],
                [T::zero(), T::zero(), T::zero(), T::one()],
            ]),
            Matrix4x4::init([
                [T::one() / vec.x, T::zero(), T::zero(), T::zero()],
                [T::zero(), T::one() / vec.y, T::zero(), T::zero()],
                [T::zero(), T::zero(), T::one() / vec.z, T::zero()],
                [T::zero(), T::zero(), T::zero(), T::one()],
            ]),
        )
    }

    pub fn rotate_x(angle: T) -> Self {
        let c = angle.cos();
        let s = angle.sin();
        Self::init(
            Matrix4x4::init([
                [T::one(), T::zero(), T::zero(), T::zero()],
                [T::zero(), c, -s, T::zero()],
                [T::zero(), s, c, T::zero()],
                [T::zero(), T::zero(), T::zero(), T::one()],
            ]),
            Matrix4x4::init([
                [T::one(), T::zero(), T::zero(), T::zero()],
                [T::zero(), c, s, T::zero()],
                [T::zero(), -s, c, T::zero()],
                [T::zero(), T::zero(), T::zero(), T::one()],
            ]),
        )
    }

    pub fn rotate_y(angle: T) -> Self {
        let c = angle.cos();
        let s = angle.sin();
        Self::init(
            Matrix4x4::init([
                [c, T::zero(), s, T::zero()],
                [T::zero(), T::one(), T::zero(), T::zero()],
                [-s, T::zero(), c, T::zero()],
                [T::zero(), T::zero(), T::zero(), T::one()],
            ]),
            Matrix4x4::init([
                [c, T::zero(), -s, T::zero()],
                [T::zero(), T::one(), T::zero(), T::zero()],
                [s, T::zero(), c, T::zero()],
                [T::zero(), T::zero(), T::zero(), T::one()],
            ]),
        )
    }

    pub fn rotate_z(angle: T) -> Self {
        let c = angle.cos();
        let s = angle.sin();
        Self::init(
            Matrix4x4::init([
                [c, -s, T::zero(), T::zero()],
                [s, c, T::zero(), T::zero()],
                [T::zero(), T::zero(), T::one(), T::zero()],
                [T::zero(), T::zero(), T::zero(), T::one()],
            ]),
            Matrix4x4::init([
                [c, s, T::zero(), T::zero()],
                [-s, c, T::zero(), T::zero()],
                [T::zero(), T::zero(), T::one(), T::zero()],
                [T::zero(), T::zero(), T::zero(), T::one()],
            ]),
        )
    }

    pub fn rotate(vec: Vector3<T>, angle: T) -> Self {
        // Rodrigues' rotation formula
        let u = vec.normalize();

        let x = Vector3::init(T::one(), T::zero(), T::zero());
        let y = Vector3::init(T::zero(), T::one(), T::zero());
        let z = Vector3::init(T::zero(), T::zero(), T::one());

        let c = angle.cos();
        let s = angle.sin();

        let x_rotated = x * c + u * (T::one() - c) * (u.dot(x)) + u.cross(x) * s;
        let y_rotated = y * c + u * (T::one() - c) * (u.dot(y)) + u.cross(y) * s;
        let z_rotated = z * c + u * (T::one() - c) * (u.dot(z)) + u.cross(z) * s;

        let mat = Matrix4x4::init([
            [x_rotated.x, y_rotated.x, z_rotated.x, T::zero()],
            [x_rotated.y, y_rotated.y, z_rotated.y, T::zero()],
            [x_rotated.z, y_rotated.z, z_rotated.z, T::zero()],
            [T::zero(), T::zero(), T::zero(), T::one()],
        ]);

        Self::init(mat, mat.transpose())
    }

    pub fn mul(&self, rhs: &Self) -> Self {
        Self::init(self.m.mul(&rhs.m), rhs.minv.mul(&self.minv))
    }

    pub fn look_at(pos: Point3<T>, target: Point3<T>, up: Vector3<T>) -> Self {
        let up_u = up.normalize();
        let z_u = (target - pos).normalize();
        let x_u = up_u.cross(z_u).normalize();
        let y_u = z_u.cross(x_u);
        Self::init_from_orthogonal_basis(target, x_u, y_u, z_u)
    }

    pub fn matrix(&self) -> Matrix4x4<T> {
        self.m
    }

    pub fn inverse_matrix(&self) -> Matrix4x4<T> {
        self.minv
    }
}

impl<T: GFloat> std::ops::Mul<&Transform<T>> for &Transform<T> {
    type Output = Transform<T>;

    fn mul(self, rhs: &Transform<T>) -> Self::Output {
        self.mul(rhs)
    }
}

pub trait Invertible {
    fn inverse(&self) -> Self;
}

impl<T: GFloat> Invertible for Matrix4x4<T> {
    fn inverse(&self) -> Self {
        Self::init(inverse_impl(&self.data))
    }
}

impl<T: GFloat> Invertible for Transform<T> {
    fn inverse(&self) -> Self {
        Self::init(self.minv, self.m)
    }
}

pub trait TransformMap<T> {
    type Output;

    fn map(&self, obj: T) -> Self::Output;

    fn map_inverse(&self, obj: T) -> Self::Output
    where
        Self: Invertible + Sized,
    {
        self.inverse().map(obj)
    }
}

impl<T, U> TransformMap<Vector3<U>> for Transform<T>
where
    T: GFloat,
    U: Scalar + std::ops::Mul<T, Output = U>,
{
    type Output = Vector3<U>;

    fn map(&self, vec: Vector3<U>) -> Self::Output {
        map_vector(&self.m, &vec)
    }

    fn map_inverse(&self, vec: Vector3<U>) -> Self::Output {
        map_vector(&self.minv, &vec)
    }
}

fn map_vector<T, U>(m: &Matrix4x4<T>, vec: &Vector3<U>) -> Vector3<U>
where
    T: GFloat,
    U: Scalar + std::ops::Mul<T, Output = U>,
{
    // TODO: assert that w is nearly 0
    // let w = m[(3, 0)] * p.x + m[(3, 1)] * p.y + m[(3, 2)] * p.z;
    Vector3::init(
        vec.x * m[(0, 0)] + vec.y * m[(0, 1)] + vec.z * m[(0, 2)],
        vec.x * m[(1, 0)] + vec.y * m[(1, 1)] + vec.z * m[(1, 2)],
        vec.x * m[(2, 0)] + vec.y * m[(2, 1)] + vec.z * m[(2, 2)],
    )
}

impl<T, U> TransformMap<Point3<U>> for Transform<T>
where
    T: GFloat,
    U: Scalar
        + std::ops::Mul<T, Output = U>
        + std::ops::Div<T, Output = U>
        + std::ops::Add<T, Output = U>,
{
    type Output = Point3<U>;

    fn map(&self, p: Point3<U>) -> Self::Output {
        map_point(&self.m, &p)
    }

    fn map_inverse(&self, p: Point3<U>) -> Self::Output {
        map_point(&self.minv, &p)
    }
}

fn map_point<T, U>(m: &Matrix4x4<T>, p: &Point3<U>) -> Point3<U>
where
    T: GFloat,
    U: Scalar
        + std::ops::Mul<T, Output = U>
        + std::ops::Div<T, Output = U>
        + std::ops::Add<T, Output = U>,
{
    let w = p.x * m[(3, 0)] + p.y * m[(3, 1)] + p.z * m[(3, 2)] + m[(3, 3)];
    assert!(w != U::zero());
    Point3::init(
        (p.x * m[(0, 0)] + p.y * m[(0, 1)] + p.z * m[(0, 2)] + m[(0, 3)]) / w,
        (p.x * m[(1, 0)] + p.y * m[(1, 1)] + p.z * m[(1, 2)] + m[(1, 3)]) / w,
        (p.x * m[(2, 0)] + p.y * m[(2, 1)] + p.z * m[(2, 2)] + m[(2, 3)]) / w,
    )
}

impl<T: GFloat> TransformMap<Bbox3<T>> for Transform<T> {
    type Output = Bbox3<T>;

    fn map(&self, bbox: Bbox3<T>) -> Self::Output {
        map_bbox(&self.m, &bbox)
    }

    fn map_inverse(&self, bbox: Bbox3<T>) -> Self::Output {
        map_bbox(&self.minv, &bbox)
    }
}

fn map_bbox<T: GFloat>(m: &Matrix4x4<T>, bbox: &Bbox3<T>) -> Bbox3<T> {
    // Graphics Gems: Transforming Axis-aligned Bounding Boxes
    // by James Arvo 1990
    let mut pmin = m.col(3);
    let mut pmax = pmin.clone();

    let bbox_pmin: [T; 3] = [bbox.pmin.x, bbox.pmin.y, bbox.pmin.z];
    let bbox_pmax: [T; 3] = [bbox.pmax.x, bbox.pmax.y, bbox.pmax.z];

    for i in 0..3 {
        for j in 0..3 {
            let a = m[i][j] * bbox_pmin[j];
            let b = m[i][j] * bbox_pmax[j];

            if a < b {
                pmin[j] += a;
                pmax[j] += b;
            } else {
                pmin[j] += b;
                pmax[j] += a;
            }
        }
    }

    Bbox3::init(
        Point3::init(pmin[0], pmin[1], pmin[2]),
        Point3::init(pmax[0], pmax[1], pmax[2]),
    )
}
impl<T, U> TransformMap<Normal3<U>> for Transform<T>
where
    T: GFloat,
    U: Scalar + std::ops::Mul<T, Output = U>,
{
    type Output = Normal3<U>;

    fn map(&self, n: Normal3<U>) -> Self::Output {
        map_normal(&self.minv, &n)
    }

    fn map_inverse(&self, n: Normal3<U>) -> Self::Output {
        map_normal(&self.m, &n)
    }
}

fn map_normal<T, U>(minv: &Matrix4x4<T>, n: &Normal3<U>) -> Normal3<U>
where
    T: GFloat,
    U: Scalar + std::ops::Mul<T, Output = U>,
{
    // normals are mapped by the inverse transpose of the transformation matrix
    Normal3::init(
        n.x * minv[(0, 0)] + n.y * minv[(1, 0)] + n.z * minv[(2, 0)] * n.z,
        n.x * minv[(0, 1)] + n.y * minv[(1, 1)] + n.z * minv[(2, 1)] * n.z,
        n.x * minv[(0, 2)] + n.y * minv[(1, 2)] + n.z * minv[(2, 2)] * n.z,
    )
}

impl<T, U> TransformMap<Ray<U>> for Transform<T>
where
    T: GFloat,
    U: Scalar
        + std::ops::Mul<T, Output = U>
        + std::ops::Div<T, Output = U>
        + std::ops::Add<T, Output = U>,
{
    type Output = Ray<U>;

    fn map(&self, ray: Ray<U>) -> Self::Output {
        Ray::init(self.map(ray.origin), self.map(ray.direction), ray.tmax)
    }

    fn map_inverse(&self, ray: Ray<U>) -> Self::Output {
        Ray::init(
            self.map_inverse(ray.origin),
            self.map_inverse(ray.direction),
            ray.tmax,
        )
    }
}

pub type Matrix4x4f = Matrix4x4<Float>;
pub type Matrix4x4i = Matrix4x4<Int>;
pub type Transformf = Transform<Float>;
