use crate::scalar::*;

#[derive(Debug, Clone, Copy)]
pub struct Interval<T: GFloat>(T);

pub type Intervalf = Interval<Float>;

impl<T: GFloat> Interval<T> {
    pub fn init(inf: T, sup: T) -> Self {
        assert!(sup >= inf);
        if inf == sup {
            Self { 0: inf }
        } else {
            let val = (inf + sup) / T::from(2).unwrap();
            Self { 0: val }
        }
    }
    pub fn is_positive(&self) -> bool {
        self.0 > T::zero()
    }
    pub fn is_nonnegative(&self) -> bool {
        self.0 >= T::zero()
    }
    pub fn is_negative(&self) -> bool {
        self.0 < T::zero()
    }
    pub fn square(&self) -> Self {
        let res = self.0 * self.0;
        Self { 0: res }
    }
    pub fn sqrt(&self) -> Option<Self> {
        if self.is_nonnegative() {
            let res = self.0.sqrt();
            Some(Self { 0: res })
        } else {
            None
        }
    }
    pub fn approx(&self) -> T {
        self.0
    }
    pub fn in_range(&self, t0: T, t1: T) -> bool {
        self.0 >= t0 && self.0 < t1
    }
}

impl<T: GFloat> From<T> for Interval<T> {
    fn from(t: T) -> Interval<T> {
        Interval::init(t, t)
    }
}
impl<T: GFloat> std::ops::Add<Self> for Interval<T> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        let res = self.0 + rhs.0;
        Self { 0: res }
    }
}
impl<T: GFloat> std::ops::Sub<Self> for Interval<T> {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        let res = self.0 - rhs.0;
        Self { 0: res }
    }
}
impl<T: GFloat> std::ops::Add<T> for Interval<T> {
    type Output = Self;

    fn add(self, rhs: T) -> Self::Output {
        let res = self.0 + rhs;
        Self { 0: res }
    }
}
impl<T: GFloat> std::ops::Sub<T> for Interval<T> {
    type Output = Self;

    fn sub(self, rhs: T) -> Self::Output {
        let res = self.0 - rhs;
        Self { 0: res }
    }
}
impl<T: GFloat> std::ops::Mul<Self> for Interval<T> {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        let res = self.0 * rhs.0;
        Self { 0: res }
    }
}
impl<T: GFloat> std::ops::Mul<T> for Interval<T> {
    type Output = Self;

    fn mul(self, rhs: T) -> Self::Output {
        let res = self.0 * rhs;
        Self { 0: res }
    }
}
impl<T: GFloat> std::ops::Div<Self> for Interval<T> {
    type Output = Self;

    fn div(self, rhs: Self) -> Self::Output {
        let res = self.0 / rhs.0;
        Self { 0: res }
    }
}
impl<T: GFloat> std::ops::Div<T> for Interval<T> {
    type Output = Self;

    fn div(self, rhs: T) -> Self::Output {
        let res = self.0 / rhs;
        Self { 0: res }
    }
}
impl<T: GFloat> std::ops::AddAssign<Self> for Interval<T> {
    fn add_assign(&mut self, rhs: Self) {
        *self = *self + rhs;
    }
}
impl<T: GFloat> std::ops::SubAssign<Self> for Interval<T> {
    fn sub_assign(&mut self, rhs: Self) {
        *self = *self - rhs;
    }
}
impl<T: GFloat> std::ops::MulAssign<Self> for Interval<T> {
    fn mul_assign(&mut self, rhs: Self) {
        *self = *self * rhs;
    }
}
impl<T: GFloat> std::ops::DivAssign<Self> for Interval<T> {
    fn div_assign(&mut self, rhs: Self) {
        *self = *self / rhs;
    }
}
impl<T: GFloat> std::ops::Neg for Interval<T> {
    type Output = Self;

    fn neg(self) -> Self::Output {
        let res = -self.0;
        Self { 0: res }
    }
}
impl<T: GFloat> std::cmp::PartialEq<Self> for Interval<T> {
    fn eq(&self, other: &Self) -> bool {
        self.0 == other.0
    }
}
impl<T: GFloat> std::cmp::PartialOrd<Self> for Interval<T> {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        if self == other {
            Some(std::cmp::Ordering::Equal)
        } else if self.0 < other.0 {
            Some(std::cmp::Ordering::Less)
        } else if other.0 < self.0 {
            Some(std::cmp::Ordering::Greater)
        } else {
            None
        }
    }
}
impl<T: GFloat> num_traits::Zero for Interval<T> {
    fn zero() -> Self {
        Self::init(T::zero(), T::zero())
    }
    fn is_zero(&self) -> bool {
        self.0 == T::zero()
    }
}
impl<T: GFloat> num_traits::One for Interval<T> {
    fn one() -> Self {
        Self::init(T::one(), T::one())
    }
}
impl<T: GFloat> LowestHighest for Interval<T> {
    fn lowest() -> Self {
        Self::from(T::lowest())
    }
    fn highest() -> Self {
        Self::from(T::highest())
    }
}
impl<T: GFloat> num_traits::Num for Interval<T> {
    type FromStrRadixErr = T::FromStrRadixErr;
    fn from_str_radix(s: &str, r: u32) -> std::result::Result<Self, Self::FromStrRadixErr> {
        match T::from_str_radix(s, r) {
            Err(e) => Err(e),
            Ok(s) => Ok(Self::from(s)),
        }
    }
}
impl<T: GFloat> std::ops::Rem<Self> for Interval<T> {
    type Output = Self;
    fn rem(self, _rhs: Self) -> Self::Output {
        std::unimplemented!();
    }
}
impl<T: GFloat> std::ops::RemAssign<Self> for Interval<T> {
    fn rem_assign(&mut self, _rhs: Self) -> () {
        std::unimplemented!();
    }
}

// Should this be here?
pub fn solve_quadratic<T: GFloat>(
    a: Interval<T>,
    b: Interval<T>,
    c: Interval<T>,
) -> Option<(Interval<T>, Interval<T>)> {
    // discriminant
    let d = b.square() - a * c * T::from(4).unwrap();
    match d.sqrt() {
        None => None,
        Some(sqrt_d) => {
            // https://people.csail.mit.edu/bkph/articles/Quadratics.pdf
            if b.is_nonnegative() {
                Some(min_max(
                    -(b + sqrt_d) * T::from(0.5).unwrap() / a,
                    -(c + c) / (b + sqrt_d),
                ))
            } else {
                Some(min_max(
                    (-b + sqrt_d) * T::from(0.5).unwrap() / a,
                    (c + c) / (-b + sqrt_d),
                ))
            }
        }
    }
}
