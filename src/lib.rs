//! Galois Field (2<sup>M</sup>) arithmetic support.
//!
//! galois_field_2pm provides wrappers for the u8, u16, u32, u64, and u128 types that implement Galois Field arithmetic. There are two implementations:
//! - A look up table based implementation for multiplication and division for Galois Fields defined by a primitive polynomial
//! - A computation based implementation for multiplication and division for Galois Fields defined by an irreducible polynomial
//!
//! # Getting Started
//!
//! To start using a Galois Field defined by the irreducible polynomial p(x) ∈ GF(2)\[x\] start by representing p(x) as a u128.
//!
//! An irrecuible polynomial of degree M defines the field GF(2<sup>M</sup>) = GF(2)\[x\] / p(x)
//!
//! Next we need to choose how the field is implemented:
//!
//!   - If p(x) is primitive and M <= 16 then the look up table implementation can be used (module gf2_lut)
//!
//!   - Else use the computation based implementation (module gf2)
//!
//! Lastly we must use one of the structs to represent the elements in the field. The struct GFuX can be used for M ≤ X.
//!
//! At this point we can define the GF type and use it.
//!
//! Consider the example below
//! ```
//! use galois_field_2pm::{GaloisField, gf2_lut};
//!
//! // First lets define out Galois Field
//! // For the irreducible polynomial p(x) = x^3 + x + 1 the u128 representation is given by 0xB.
//! // Note that this p(x) is primitive so we can use either gf2_lut of gf2
//! // Note that GFu8 can represent the field generate by 0xB as the degree of 0xB is 3 and 3 ≤ 8.
//! //     For p(x) = 0x211 (degree 9), GFu8 can not be used. GFu16 must be used instead
//! type GF = gf2_lut::GFu8::<0xB>; // Or alternatively type GF = gf2::GFu8::<0xB>;
//! let a = GF::ONE;
//! let b = GF::new(2);
//! let c = GF {value: 3};
//!
//! let d = (a + b * c).inverse();
//! ```

#[cfg(test)]
use paste::paste;

#[cfg(test)]
use rand::Rng;

use core::fmt::{Debug, Display};
use core::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Sub, SubAssign};

pub mod gf2;
pub mod gf2_lut;

/// A trait used to indicate that a type can be used to represent the elements of a Galois Field.
pub trait GaloisField:
    Clone
    + Copy
    + Eq
    + PartialEq
    + Debug
    + Display
    + Add
    + AddAssign
    + Sub
    + SubAssign
    + Mul
    + MulAssign
    + Div
    + DivAssign
{
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

    /// Takes the inverse of an element in the field. Panics when inverting the additive identity (aka zero)
    fn inverse(&self) -> Self;

    /// Constructs a GF element using the underlying storage type
    fn new(value: Self::StorageType) -> Self;

    /// Used to check if the value stored is a valid element in the current field
    fn validate(&self) -> bool;
}

#[allow(dead_code)]
const fn calc_degree(x: u128) -> i16 {
    127 - (x.leading_zeros() as i16)
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
        };
    }

    macro_rules! inverse_multiplication_test {
        ($type:ty) => {
            for i in 1..GF::NUM_ELEM {
                let a = GF::new(i as $type);
                assert_eq!(a * a.inverse(), GF::ONE);
                assert_eq!(a / a, GF::ONE);
            }
        };
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
        };
    }

    macro_rules! field_test {
        ($($mod:tt: $type:ty: $poly:expr,)*) => {
        $(
            paste! {
                #[test]
                fn [<$mod _gf_ $poly>]() {
                    type GF = $mod::[<GF $type>]<$poly>;
                    addition_test!($type);
                    multiplication_test!($type);
                    distributive_test!($type);
                }

                #[test]
                #[should_panic]
                fn [<$mod _div_by_0_ $poly>]() {
                    type GF = $mod::[<GF $type>]<$poly>;
                    let _ = GF::ONE / GF::ZERO;
                }

                #[test]
                #[should_panic]
                fn [<$mod _inverse_of_0_ $poly>]() {
                    type GF = $mod::[<GF $type>]<$poly>;
                    let _ = GF::ZERO.inverse();
                }
            }
        )*
        }
    }

    field_test! {
        gf2_lut: u16: 0x3,
        gf2_lut: u16: 0x7,
        gf2_lut: u16: 0xb,
        gf2_lut: u16: 0x13,
        gf2_lut: u8: 0x25,
        gf2_lut: u8: 0x43,
        gf2_lut: u8: 0x83,
        gf2_lut: u8: 0x11d,

        gf2: u16: 0x3,
        gf2: u16: 0x7,
        gf2: u16: 0xb,
        gf2: u16: 0x13,
        gf2: u8: 0x25,
        gf2: u8: 0x43,
        gf2: u8: 0x83,
        gf2: u8: 0x11d,
    }

    macro_rules! associative_spot_test {
        ($type:ty, $a:tt, $b:tt, $c:tt, $op:tt) => {
            for i in 0..NUM_VALS {
                for j in 0..NUM_VALS {
                    for k in 0..NUM_VALS {
                        assert_eq!($a[i] $op ($b[j] $op $c[k]), ($a[i] $op $b[j]) $op $c[k]);
                    }
                }
            }
        }
    }

    macro_rules! commutative_spot_test {
        ($type:ty, $a:tt, $b:tt, $op:tt) => {
            for i in 0..NUM_VALS {
                for j in 0..NUM_VALS {
                    assert_eq!($a[i] $op $b[j], $b[j] $op $a[i]);
                }
            }
        }
    }

    macro_rules! identity_spot_test {
        ($type:ty, $a:tt, $op:tt, $identity:tt) => {
            for i in 0..NUM_VALS {
                assert_eq!($a[i] $op GF::$identity, $a[i]);
            }
        }
    }

    macro_rules! inverse_addition_spot_test {
        ($type:ty, $a:tt) => {
            for i in 0..NUM_VALS {
                assert_eq!($a[i] + $a[i], GF::ZERO);
            }
        };
    }

    macro_rules! inverse_multiplication_spot_test {
        ($type:ty, $a:tt) => {
            for i in 0..NUM_VALS {
                assert_eq!($a[i] * $a[i].inverse(), GF::ONE);
                assert_eq!($a[i] / $a[i], GF::ONE);
            }
        };
    }

    macro_rules! addition_spot_test {
        ($type:ty, $a:tt, $b:tt, $c:tt) => {
            associative_spot_test!($type, $a, $b, $c, +);
            commutative_spot_test!($type, $a, $b, +);
            identity_spot_test!($type, $a, +, ZERO);
            inverse_addition_spot_test!($type, $a);
        }
    }

    macro_rules! multiplication_spot_test {
        ($type:ty, $a:tt, $b:tt, $c:tt) => {
            associative_spot_test!($type, $a, $b, $c, *);
            commutative_spot_test!($type, $a, $b, *);
            identity_spot_test!($type, $a, *, ONE);
            inverse_multiplication_spot_test!($type, $a);
        }
    }

    macro_rules! distributive_spot_test {
        ($type:ty, $a:tt, $b:tt, $c:tt) => {
            for i in 0..NUM_VALS {
                for j in 0..NUM_VALS {
                    for k in 0..NUM_VALS {
                        assert_eq!($a[i] * ($b[j] + $c[k]), $a[i] * $b[j] + $a[i] * $c[k]);
                    }
                }
            }
        };
    }

    macro_rules! field_spot_test {
        ($($mod:tt: $type:ty: $poly:expr,)*) => {
        $(
            paste! {
                #[test]
                fn [<$mod _spot_gf_ $poly>]() {
                    type GF = $mod::[<GF $type>]<$poly>;

                    const NUM_VALS: usize = 30;
                    let mut a: [GF; NUM_VALS] = [GF::ZERO; NUM_VALS];
                    let mut b: [GF; NUM_VALS] = [GF::ZERO; NUM_VALS];
                    let mut c: [GF; NUM_VALS] = [GF::ZERO; NUM_VALS];
                    for i in 0..NUM_VALS {
                        a[i] = GF::new(rand::thread_rng().gen_range(1..GF::NUM_ELEM) as $type);
                        b[i] = GF::new(rand::thread_rng().gen_range(1..GF::NUM_ELEM) as $type);
                        c[i] = GF::new(rand::thread_rng().gen_range(1..GF::NUM_ELEM) as $type);
                    }

                    addition_spot_test!($type, a, b, c);
                    multiplication_spot_test!($type, a, b, c);
                    distributive_spot_test!($type, a, b, c);
                }
            }
        )*
        }
    }

    field_spot_test! {
        gf2_lut: u16: 0x211,
        gf2_lut: u16: 0x409,
        gf2_lut: u16: 0x805,
        gf2_lut: u16: 0x1053,

        gf2: u16: 0x211,
        gf2: u16: 0x409,
        gf2: u16: 0x805,
        gf2: u16: 0x1053,
        gf2: u32: 0x20009,
        gf2: u64: 0x2_0000_2001,
        gf2: u128: 0x2_0000_0000_0004_0001,
        gf2: u128: 0x8000_0000_0000_0000_0000_0000_0000_0003,
    }
}
