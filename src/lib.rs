use std::fmt::{Debug, Display};
use std::ops::{Add, Sub, Mul, Div, AddAssign, SubAssign, MulAssign, DivAssign};

pub mod gf2_lut;

pub trait GaloisField: Clone + Copy + Eq + PartialEq + Debug + Display +
                       Add + Sub + Mul + Div + AddAssign + SubAssign + MulAssign + DivAssign {
    fn validate_value(&self) -> bool;
}

#[allow(dead_code)]
const fn calc_degree(x: u128) -> i8 {
    if x > 0 {
        return calc_degree(x >> 1) + 1;
    } else {
        return -1;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn calc_degree_test() {
        assert_eq!(calc_degree(0), -1);
        for i in 0..=127 {
            assert_eq!(calc_degree(1u128 << i), i);
        }
    }
}
