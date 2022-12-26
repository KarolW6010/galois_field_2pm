use core::fmt;
use core::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Sub, SubAssign};
use paste::paste;

mod clmul;
mod gf2_poly_div;

use crate::GaloisField;
use clmul::CarryLessMultiply;
use gf2_poly_div::GF2PolyDiv;

macro_rules! assign_operator_impl {
    ($($type:ty: $trait_name:ident: $trait_fn:ident: $op:tt,)*) => {
    $(
        paste! {
            impl<const POLY: u128> $trait_name for [<GF $type>]<POLY> {
                fn $trait_fn(&mut self, other: Self) {
                    *self = *self $op other;
                }
            }
        }
    )*
    }
}

macro_rules! setup_gf {
    ($($type:ty,)*) => {
    $(
        paste! {
            // Define the struct
            #[repr(transparent)]
            #[derive(Clone, Copy, Eq, PartialEq)]
            pub struct [<GF $type>]<const POLY: u128> {
                pub value: $type,
            }

            // Implement the traits
            impl<const POLY: u128> GaloisField for [<GF $type>]<POLY> {
                type StorageType = $type;

                const M: u128 = crate::calc_degree(POLY) as u128;
                const NUM_ELEM: u128 = 1 << Self::M;

                const ZERO: Self = Self {value: 0};
                const ONE: Self = Self {value: 1};

                fn inverse(&self) -> Self {
                    if *self == Self::ZERO {
                        panic!("Cannot take inverse of zero");
                    }
                    if *self == Self::ONE {
                        return Self::ONE;
                    }

                    let mut temp_t: Self::StorageType;
                    let mut t: Self::StorageType = 0;
                    let mut new_t: Self::StorageType = 1;

                    let mut quo = Self::ZERO;
                    let mut nt = Self::ZERO;

                    let mut new_remainder: Self::StorageType;
                    let mut remainder = self.value;
                    (quo.value, new_remainder) = $type::gf2_poly_div_poly(POLY, self.value);

                    nt.value = new_t;
                    temp_t = t;
                    t = new_t;
                    new_t = temp_t ^ ((quo * nt).value);

                    while new_remainder != 0 {
                        let temp_r = remainder;
                        remainder = new_remainder;
                        (quo.value, new_remainder) = $type::gf2_poly_div(temp_r, new_remainder);

                        nt.value = new_t;
                        temp_t = t;
                        t = new_t;
                        new_t = temp_t ^ ((quo * nt).value);
                    }

                    Self {
                        value: t
                    }
                }

                fn new(value: $type) -> Self {
                    Self {value: value}
                }

                fn validate(&self) -> bool {
                    (self.value as u128) >= Self::NUM_ELEM
                }
            }

            impl<const POLY: u128> fmt::Debug for [<GF $type>]<POLY> {
                fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
                    write!(f, "GF<{:#0X}>(value: 0x{:0width$x})", POLY, self.value, width = (Self::M as usize / 4))
                }
            }

            impl<const POLY: u128> fmt::Display for [<GF $type>]<POLY> {
                fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
                    write!(f, "0x{:0width$X}", self.value, width = (Self::M as usize / 4))
                }
            }

            impl<const POLY: u128> Add<[<GF $type>]<POLY>> for [<GF $type>]<POLY> {
                type Output = Self;

                fn add(self, other: Self) -> Self {
                    Self {
                        value: self.value ^ other.value,
                    }
                }
            }

            impl<const POLY: u128> Sub<[<GF $type>]<POLY>> for [<GF $type>]<POLY> {
                type Output = Self;

                fn sub(self, other: Self) -> Self {
                    Self {
                        value: self.value ^ other.value,
                    }
                }
            }

            impl<const POLY: u128> Mul<[<GF $type>]<POLY>> for [<GF $type>]<POLY> {
                type Output = Self;

                fn mul(self, other: Self) -> Self {
                    let hi = self.value.clmul_high(other.value);
                    let lo = self.value.clmul_low(other.value);

                    Self {
                        value: $type::gf2_poly_mod(hi, lo, POLY),
                    }
                }
            }

            impl<const POLY: u128> Div<[<GF $type>]<POLY>> for [<GF $type>]<POLY> {
                type Output = Self;

                fn div(self, other: Self) -> Self {
                    self * other.inverse()
                }
            }

            assign_operator_impl! {
                $type: AddAssign: add_assign: +,
                $type: SubAssign: sub_assign: -,
                $type: MulAssign: mul_assign: *,
                $type: DivAssign: div_assign: /,
            }
        }
    )*
    };
}

setup_gf! {
    u8,
    u16,
    u32,
    u64,
    u128,
}
