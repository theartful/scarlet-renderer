use crate::scalar::*;

pub trait InnerScalar {
    type ScalarType;
}

pub trait VectorTraits:
    InnerScalar
    + std::ops::Mul<<Self as InnerScalar>::ScalarType, Output = Self>
    + std::ops::Div<<Self as InnerScalar>::ScalarType, Output = Self>
    + std::ops::Add<Self, Output = Self>
    + std::ops::Sub<Self, Output = Self>
    + std::ops::Neg
    + Copy
    + Sized
{
}

impl<T> VectorTraits for T where
    T: InnerScalar
        + std::ops::Mul<<Self as InnerScalar>::ScalarType, Output = Self>
        + std::ops::Div<<Self as InnerScalar>::ScalarType, Output = Self>
        + std::ops::Add<Self, Output = Self>
        + std::ops::Sub<Self, Output = Self>
        + std::ops::Neg
        + Copy
        + Sized
{
}

pub trait InnerProduct<Rhs>: VectorTraits {
    fn dot(&self, rhs: Rhs) -> Self::ScalarType;
}

pub trait Norm: InnerProduct<Self> {
    fn square_norm(self) -> Self::ScalarType;

    fn norm(self) -> Self::ScalarType
    where
        Self::ScalarType: GFloat,
    {
        self.square_norm().sqrt()
    }

    fn square_distance(self, vec: Self) -> Self::ScalarType {
        (self - vec).square_norm()
    }

    fn distance(self, vec: Self) -> Self::ScalarType
    where
        Self::ScalarType: GFloat,
        Self: std::ops::Sub<Self, Output = Self>,
    {
        self.square_distance(vec).sqrt()
    }

    fn normalize(self) -> Self
    where
        Self::ScalarType: GFloat,
    {
        self / self.norm()
    }
}

impl<T> Norm for T
where
    T: InnerProduct<T>,
{
    fn square_norm(self) -> Self::ScalarType {
        self.dot(self)
    }
}
