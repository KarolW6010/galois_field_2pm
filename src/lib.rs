use std::fmt::{Debug, Display};
use std::ops::{Add, Sub, Mul, Div, AddAssign, SubAssign, MulAssign, DivAssign};

pub mod gf2_lut;

pub trait GaloisField: Clone + Copy + Eq + PartialEq + Debug + Display +
                       Add + Sub + Mul + Div + AddAssign + SubAssign + MulAssign + DivAssign {}

