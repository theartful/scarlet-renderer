use crate::ray::*;
pub use crate::scalar::*;
use crate::vector::*;

#[derive(Debug, Clone, Copy)]
pub struct Bbox3<T: Scalar> {
    pub pmin: Point3<T>,
    pub pmax: Point3<T>,
}

pub trait Bounded3<T: Scalar> {
    fn bbox(&self) -> Bbox3<T>;

    fn xmin(&self) -> T {
        self.bbox().pmin.x
    }
    fn ymin(&self) -> T {
        self.bbox().pmin.y
    }
    fn zmin(&self) -> T {
        self.bbox().pmin.z
    }
    fn xmax(&self) -> T {
        self.bbox().pmax.x
    }
    fn ymax(&self) -> T {
        self.bbox().pmax.y
    }
    fn zmax(&self) -> T {
        self.bbox().pmax.z
    }
}

impl<T: Scalar> Bounded3<T> for Bbox3<T> {
    fn bbox(&self) -> Bbox3<T> {
        *self
    }
    fn xmin(&self) -> T {
        self.pmin.x
    }
    fn ymin(&self) -> T {
        self.pmin.y
    }
    fn zmin(&self) -> T {
        self.pmin.z
    }
    fn xmax(&self) -> T {
        self.pmax.x
    }
    fn ymax(&self) -> T {
        self.pmax.y
    }
    fn zmax(&self) -> T {
        self.pmax.z
    }
}

impl<T: Scalar> Bounded3<T> for Point3<T> {
    fn bbox(&self) -> Bbox3<T> {
        Bbox3::<T>::init(
            Point3::<T>::init(self.x, self.y, self.z),
            Point3::<T>::init(self.x, self.y, self.z),
        )
    }

    fn xmin(&self) -> T {
        self.x
    }
    fn ymin(&self) -> T {
        self.y
    }
    fn zmin(&self) -> T {
        self.z
    }
    fn xmax(&self) -> T {
        self.x
    }
    fn ymax(&self) -> T {
        self.y
    }
    fn zmax(&self) -> T {
        self.z
    }
}

impl<T: Scalar> Bbox3<T> {
    pub fn new() -> Self {
        let highest = T::highest();
        let lowest = T::lowest();
        Bbox3 {
            pmin: Point3::<T>::init(highest, highest, highest),
            pmax: Point3::<T>::init(lowest, lowest, lowest),
        }
    }

    pub fn init(p1: Point3<T>, p2: Point3<T>) -> Self {
        Bbox3 {
            pmin: Point3::<T>::init(min(p1.x, p2.x), min(p1.y, p2.y), min(p1.z, p2.z)),
            pmax: Point3::<T>::init(max(p1.x, p2.x), max(p1.y, p2.y), max(p1.z, p2.z)),
        }
    }

    pub fn union(&self, rhs: Bbox3<T>) -> Self {
        Bbox3 {
            pmin: Point3::<T>::init(
                min(self.xmin(), rhs.xmin()),
                min(self.ymin(), rhs.ymin()),
                min(self.zmin(), rhs.zmin()),
            ),
            pmax: Point3::<T>::init(
                max(self.xmax(), rhs.xmax()),
                max(self.ymax(), rhs.ymax()),
                max(self.zmax(), rhs.zmax()),
            ),
        }
    }

    pub fn intersect(&self, rhs: Bbox3<T>) -> Self {
        Bbox3 {
            pmin: Point3::<T>::init(
                max(self.xmin(), rhs.xmin()),
                max(self.ymin(), rhs.ymin()),
                max(self.zmin(), rhs.zmin()),
            ),
            pmax: Point3::<T>::init(
                min(self.xmax(), rhs.xmax()),
                min(self.ymax(), rhs.ymax()),
                min(self.zmax(), rhs.zmax()),
            ),
        }
    }

    pub fn overlaps(&self, rhs: &Bbox3<T>) -> bool {
        self.xmax() > rhs.xmin()
            && rhs.xmax() > self.xmin()
            && self.ymax() > rhs.ymin()
            && rhs.ymax() > self.ymin()
            && self.zmax() > rhs.zmin()
            && rhs.zmax() > self.zmin()
    }

    pub fn inside(&self, rhs: Point3<T>) -> bool {
        rhs.x >= self.xmin()
            && rhs.x <= self.xmax()
            && rhs.y >= self.ymin()
            && rhs.y <= self.ymax()
            && rhs.z >= self.zmin()
            && rhs.z <= self.zmax()
    }

    pub fn diagonal(&self) -> Vector3<T> {
        self.pmax - self.pmin
    }

    pub fn center(&self) -> Point3<T>
    where
        T: GFloat,
    {
        (self.pmax + self.pmin) * T::from(0.5).unwrap()
    }

    pub fn max_extent(&self) -> ((T, T), u8) {
        let x_extent = self.pmax.x - self.pmin.x;
        let y_extent = self.pmax.y - self.pmin.y;
        let z_extent = self.pmax.z - self.pmin.z;

        if x_extent > y_extent {
            if x_extent > z_extent {
                ((self.pmin.x, self.pmax.x), 0)
            } else {
                ((self.pmin.z, self.pmax.z), 2)
            }
        } else {
            if y_extent > z_extent {
                ((self.pmin.y, self.pmax.y), 1)
            } else {
                ((self.pmin.z, self.pmax.z), 2)
            }
        }
    }
}

impl<T: Scalar> std::ops::Add<Bbox3<T>> for Bbox3<T> {
    type Output = Bbox3<T>;

    fn add(self, rhs: Bbox3<T>) -> Self::Output {
        self.union(rhs)
    }
}

impl<T: Scalar> std::ops::AddAssign<Bbox3<T>> for Bbox3<T> {
    fn add_assign(&mut self, rhs: Bbox3<T>) {
        *self = self.union(rhs);
    }
}

pub type Bbox3f = Bbox3<Float>;
pub type Bbox3i = Bbox3<Int>;

impl Bbox3f {
    pub fn intersect_ray(&self, ray: Rayf) -> Option<(Float, Float)> {
        match Self::intersect_slab(
            (Float::zero(), ray.tmax),
            (self.pmin.x, self.pmax.x),
            ray.origin.x,
            ray.direction.x,
        ) {
            Some(t_intx) => {
                match Self::intersect_slab(
                    t_intx,
                    (self.pmin.y, self.pmax.y),
                    ray.origin.y,
                    ray.direction.y,
                ) {
                    Some(t_intxy) => Self::intersect_slab(
                        t_intxy,
                        (self.pmin.z, self.pmax.z),
                        ray.origin.z,
                        ray.direction.z,
                    ),
                    None => None,
                }
            }
            None => None,
        }
    }

    pub fn intersect_raywithinfo(&self, raywithinfo: &RayWithInfof) -> Option<(Float, Float)> {
        let ray = &raywithinfo.ray;
        let ray_info = &raywithinfo.info;
        let pminmax = [self.pmin, self.pmax];
        let ray_t_int = (Float::zero(), ray.tmax);

        let t_intx_tmp = unsafe {
            (
                (pminmax.get_unchecked(ray_info.sign.x as usize).x - ray.origin.x)
                    * ray_info.reciprocal.x,
                (pminmax.get_unchecked((1 - ray_info.sign.x) as usize).x - ray.origin.x)
                    * ray_info.reciprocal.x,
            )
        };

        if ray_t_int.0 > t_intx_tmp.1 || ray_t_int.1 < t_intx_tmp.0 {
            return None;
        }

        let t_intx = (
            max(t_intx_tmp.0, ray_t_int.0),
            min(t_intx_tmp.1, ray_t_int.1),
        );

        let t_inty = unsafe {
            (
                (pminmax.get_unchecked(ray_info.sign.y as usize).y - ray.origin.y)
                    * ray_info.reciprocal.y,
                (pminmax.get_unchecked((1 - ray_info.sign.y) as usize).y - ray.origin.y)
                    * ray_info.reciprocal.y,
            )
        };

        if t_intx.0 > t_inty.1 || t_intx.1 < t_inty.0 {
            return None;
        }
        let t_intxy = (max(t_intx.0, t_inty.0), min(t_intx.1, t_inty.1));

        let t_intz = unsafe {
            (
                (pminmax.get_unchecked(ray_info.sign.z as usize).z - ray.origin.z)
                    * ray_info.reciprocal.z,
                (pminmax.get_unchecked((1 - ray_info.sign.z) as usize).z - ray.origin.z)
                    * ray_info.reciprocal.z,
            )
        };

        if t_intxy.0 > t_intz.1 || t_intxy.1 < t_intz.0 {
            return None;
        }

        Some((max(t_intxy.0, t_intz.0), min(t_intxy.1, t_intz.1)))
    }
    fn intersect_slab(
        t_int: (Float, Float),
        p_int: (Float, Float),
        o: Float,
        d: Float,
    ) -> Option<(Float, Float)> {
        let t_int2 = min_max((p_int.0 - o) / d, (p_int.1 - o) / d);
        let intersection = (max(t_int2.0, t_int.0), min(t_int2.1, t_int.1));
        if intersection.0 > intersection.1 {
            None
        } else {
            Some(intersection)
        }
    }
}
