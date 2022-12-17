#[cfg(test)]
use paste::paste;

use std::fmt::{Debug, Display};
use std::ops::{Add, Sub, Mul, Div, AddAssign, SubAssign, MulAssign, DivAssign};

pub mod gf2_lut;

pub trait GaloisField: Clone + Copy + Eq + PartialEq + Debug + Display +
                       Add + Sub + Mul + Div + AddAssign + SubAssign + MulAssign + DivAssign {
    type StorageType;

    const M: u128;
    const NUM_ELEM: u128;

    const ZERO: Self;
    const ONE: Self;

    fn inverse(&self) -> Self;

    fn new(value: Self::StorageType) -> Self;

    fn validate(&self) -> bool;

    fn validate_poly();
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

    macro_rules! associative_test {
        ($type:ty, $op:tt) => {
            for i in 0..GF::NUM_ELEM {
                for j in 0..GF::NUM_ELEM {
                    for k in 0..GF::NUM_ELEM {
                        let a = GF::new(i as $type);
                        let b = GF::new(j as $type);
                        let c = GF::new(k as $type);
                        assert_eq!(a $op (b $op c), (a $op b) $op c);
                    }
                }
            }
        }
    }

    macro_rules! commutative_test {
        ($type:ty, $op:tt) => {
            for i in 0..GF::NUM_ELEM {
                for j in 0..GF::NUM_ELEM {
                    let a = GF::new(i as $type);
                    let b = GF::new(j as $type);
                    assert_eq!(a $op b, b $op a);
                }
            }
        }
    }

    macro_rules! identity_test {
        ($type:ty, $op:tt, $identity:tt) => {
            for i in 0..GF::NUM_ELEM {
                let a = GF::new(i as $type);
                assert_eq!(a $op GF::$identity, a);
            }
        }
    }

    macro_rules! inverse_addition_test {
        ($type:ty) => {
            for i in 0..GF::NUM_ELEM {
                let a = GF::new(i as $type);
                assert_eq!(a + a, GF::ZERO);
            }
        }
    }

    macro_rules! inverse_multiplication_test {
        ($type:ty) => {
            for i in 1..GF::NUM_ELEM {
                let a = GF::new(i as $type);
                assert_eq!(a * a.inverse(), GF::ONE);
                assert_eq!(a / a, GF::ONE);
            }
        }
    }

    macro_rules! addition_test {
        ($type:ty) => {
            associative_test!($type, +);
            commutative_test!($type, +);
            identity_test!($type, +, ZERO);
            inverse_addition_test!($type);
        }
    }

    macro_rules! multiplication_test {
        ($type:ty) => {
            associative_test!($type, *);
            commutative_test!($type, *);
            identity_test!($type, *, ONE);
            inverse_multiplication_test!($type);
        }
    }

    macro_rules! distributive_test {
        ($type:ty) => {
            for i in 0..GF::NUM_ELEM {
                for j in 0..GF::NUM_ELEM {
                    for k in 0..GF::NUM_ELEM {
                        let a = GF::new(i as $type);
                        let b = GF::new(j as $type);
                        let c = GF::new(k as $type);
                        assert_eq!(a * (b + c), a * b + a * c);
                    }
                }
            }
        }
    }

    macro_rules! field_test {
        ($($type:ty: $poly:expr,)*) => {
        $(
            paste! {
                #[test]
                fn [<gf_ $poly>]() {
                    type GF = gf2_lut::[<GF $type>]<$poly>;
                    GF::validate_poly();
                    addition_test!($type);
                    multiplication_test!($type);
                    distributive_test!($type);
                }

                #[test]
                #[should_panic]
                fn [<div_by_0_ $poly>]() {
                    type GF = gf2_lut::[<GF $type>]<$poly>;
                    let _ = GF::ONE / GF::ZERO;
                }

                #[test]
                #[should_panic]
                fn [<inverse_of_0_ $poly>]() {
                    type GF = gf2_lut::[<GF $type>]<$poly>;
                    let _ = GF::ZERO.inverse();
                }
            }
        )*
        }
    }

    field_test! {
        u16: 0x3,
        u16: 0x7,
        u16: 0xb,
        u16: 0x13,
        u8: 0x25,
        u8: 0x43,
        u8: 0x83,
        u8: 0x163,
    }

    macro_rules! check_not_primitive {
        ($($type:ty: $poly:expr,)*) => {
        $(
            paste! {
                #[test]
                #[should_panic]
                fn [<not_primitive_ $poly>]() {
                    gf2_lut::[<GF $type>]::<$poly>::validate_poly();
                }
            }
        )*
        }
    }

    check_not_primitive! {
        u8: 0x4,    // Zero is a root
        u8: 0x5,    // One is a root
        u8: 0x15,   // Irreducible but not primitive    
        u16: 0x203, // Irreducible but not primitive
    }
}
