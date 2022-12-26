use core::fmt;
use core::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Sub, SubAssign};
use paste::paste;

use crate::GaloisField;

/// A trait used to indicate that the implementation of the Galois Field uses a look up table (LUT)
pub trait GaloisFieldLut: GaloisField {
    /// The primitive element of the field α. Stored as the value 2.
    const ALPHA: Self;

    /// Returns α<sup>power</sup>
    fn alpha_pow(power: isize) -> Self;

    /// For input α<sup>power</sup> returns power. For 0 returns -1.
    fn log_alpha(&self) -> isize;
}

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
                        panic!("Can not take inverse of zero");
                    } else {
                        return Self::alpha_pow(Self::DEGREE_MOD - self.log_alpha());
                    }
                }

                fn new(value: $type) -> Self {
                    Self {value: value}
                }

                fn validate(&self) -> bool {
                    (self.value as u128) >= Self::NUM_ELEM
                }
            }

            impl<const POLY: u128> GaloisFieldLut for [<GF $type>]<POLY> {
                const ALPHA: Self = Self {value: 2};

                fn alpha_pow(power: isize) -> Self {
                    let mut pow = power % Self::DEGREE_MOD;
                    pow += Self::DEGREE_MOD;
                    pow %= Self::DEGREE_MOD;
                    Self { value: Self::TABLES.exp_tbl[pow as usize] }
                }

                fn log_alpha(&self) -> isize {
                    Self::TABLES.log_tbl[self.value as usize]
                }
            }

            // Implement all the behind the scenes detail
            struct [<Tables $type:upper>] {
                exp_tbl: [$type; 1 << $type::BITS],
                log_tbl: [isize; 1 << $type::BITS],
            }

            const fn [<generate_lut_ $type:lower>](poly: u128) -> [<Tables $type:upper>] {
                let m: usize = crate::calc_degree(poly) as usize;
                let num_elems: usize = (1 << m) - 1;
                let feedback_mask: $type = 1 << (m - 1);

                let mut exp_tbl: [$type; 1 << $type::BITS] = [0; 1 << $type::BITS];
                let mut log_tbl: [isize; 1 << $type::BITS] = [0; 1 << $type::BITS];

                log_tbl[0] = -1;
                let mut value: $type = 1;
                let mut i: usize = 0;
                while i < num_elems {
                    exp_tbl[i] = value;
                    log_tbl[exp_tbl[i] as usize] = i as isize;

                    let leading_1: bool = value >= feedback_mask;
                    value <<= 1;
                    if leading_1 {
                        value ^= poly as $type;
                    }
                    i += 1;
                }

                [<Tables $type:upper>] {
                    exp_tbl: exp_tbl,
                    log_tbl: log_tbl,
                }
            }

            impl<const POLY: u128> [<GF $type>]<POLY> {
                const TABLES : [<Tables $type:upper>] = [<generate_lut_ $type:lower>](POLY);
                const DEGREE_MOD: isize = (Self::NUM_ELEM as isize) - 1;
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
                    if (self == Self::ZERO) || (other == Self::ZERO) {
                        return Self::ZERO;
                    } else {
                        return Self::alpha_pow(self.log_alpha() + other.log_alpha());
                    }
                }
            }

            impl<const POLY: u128> Div<[<GF $type>]<POLY>> for [<GF $type>]<POLY> {
                type Output = Self;

                fn div(self, other: Self) -> Self {
                    if other == Self::ZERO {
                        panic!("Divide by 0");
                    } else if self == Self::ZERO {
                        return Self::ZERO;
                    } else {
                        return Self::alpha_pow(self.log_alpha() - other.log_alpha());
                    }
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
    }
}

setup_gf! {
    u8,
    u16,
}

#[cfg(test)]
mod tests {
    use super::*;

    macro_rules! lut_specific_tests {
        ($($type:ty: $poly:expr,)*) => {
        $(
            paste! {
                #[test]
                fn [<alpha_pow_log_alpha_ $poly>]() {
                    type GF = [<GF $type>]<$poly>;
                    assert_eq!(GF::ZERO.log_alpha(), -1);

                    for i in 1..GF::DEGREE_MOD {
                        assert_eq!(GF::alpha_pow(i as isize).log_alpha(), i as isize);
                    }
                }
            }
        )*
        }
    }

    lut_specific_tests! {
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
