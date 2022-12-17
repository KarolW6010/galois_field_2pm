//! Galois Field (2<sup>M</sup>) arithmetic support.
//! 
//! galois_field_2pm provides wrappers for the u8, u16, u32, u64, and u128 types that implement Galois Field arithmetic. There are two implementations:
//! - A look up table based implementation for multiplication and division for Galois Fields defined by a primitive polynomial
//! - A computation based implementation for multiplication and division for Galois Fields defined by an irreducible polynomial
//! 
//! # Getting Started
//! 
//! To start using a Galois Field defined by the irreducible polynomial p(x) ∈ GF(2)\[x\] start by representing p(x) as a u128. 
//! An irrecuible polynomial of degree M defines the field GF(2<sup>M</sup>) = GF(2)\[x\] / p(x)
//! Next we need to choose how the field is implemented:
//!     If p(x) is primitive and M <= 16 then the look up table implementation can be used (module gf2_lut)
//!     Else use the computation based implementation (module gf2)
//! Lastly we must use one of the structs to represent the elements in the field. The struct GFuX can be used for M <= X.
//! At this point we can define the GF type and use it.
//! 
//! Consider the example below
//! ```
//! use galois_field_2pm::{GaloisField, gf2_lut};
//! 
//! // First lets define out Galois Field
//! // For irreducible polynomial p(x) = x<sup>3</sup> + x + 1 the u128 representation is given by 0xB.
//! // Note that this p(x) is primitive so we can use either gf2_lut of gf2
//! // Note that GFu8 can represent the field generate by 0xB as the degree of 0xB is 3 and 3 <= 8. For p(x) = 0x211 (degree 9), GFu8 can not be used. GFu16 must be used instead
//! type GF = gf2_lut::GFu8::<0xB>;
//! let a = GF::ONE;
//! let b = GF::new(2);
//! let c = GF {value: 3};
//! 
//! let d = (a + b * c).inverse();
//! ```
//! 
//! For example p(x) = x<sup>3</sup> + x + 1 is represented as 0xB.

#[cfg(test)]
use paste::paste;

use std::fmt::{Debug, Display};
use std::ops::{Add, Sub, Mul, Div, AddAssign, SubAssign, MulAssign, DivAssign};

pub mod gf2_lut;

/// A trait used to indicate that a type can be used to represent the elements of a Galois Field.
pub trait GaloisField: Clone + Copy + Eq + PartialEq + Debug + Display +
                       Add + Sub + Mul + Div + AddAssign + SubAssign + MulAssign + DivAssign {
    /// The underlying type used to store the representation of an element in the field
    type StorageType;

    /// The degree of the polynomial used to define the field. Specifies GF(2<sup>M</sup>)
    const M: u128;

    /// The number of elements in the field. This value is 2<sup>M</sup>
    const NUM_ELEM: u128;

    /// The additive identity of the field
    const ZERO: Self;

    /// The multiplicative identity of the field
    const ONE: Self;

    /// Takes the inverse of an element in the field. Panics when inverting the additive identity
    fn inverse(&self) -> Self;

    /// Constructs a GF element using the underlying storage type
    fn new(value: Self::StorageType) -> Self;

    /// Used to check if the value stored is a valid element in the current field
    fn validate(&self) -> bool;
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
        u8: 0x11d,
    }
}
