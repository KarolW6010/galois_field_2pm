use core::fmt;
use core::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Sub, SubAssign};
use paste::paste;

mod clmul;
use crate::GaloisField;
use clmul::{CarryLessMultiply, GF2PolyDiv, U256};

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

                    let mut temp_t: u128;
                    let mut t: u128 = 0;
                    let mut new_t: u128 = 1;

                    let mut temp_r: u128;
                    let mut r: u128 = POLY;
                    let mut new_r: u128 = self.value as u128;

                    let mut quotient: u128;
                    let mut quo = Self::ZERO;
                    let mut nt = Self::ZERO;

                    while new_r != 0 {
                        temp_r = r;
                        r = new_r;
                        (quotient, new_r) = u128::gf2_poly_div(temp_r, new_r);

                        quo.value = quotient as Self::StorageType;
                        nt.value = new_t as Self::StorageType;

                        temp_t = t;
                        t = new_t;
                        new_t = temp_t ^ ((quo * nt).value as u128);
                    }

                    Self {
                        value: t as Self::StorageType
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
                    if Self::M <= 4 {
                        return write!(f, "GF<{:#0X}>(value: 0x{:01X})", POLY, self.value);
                    } else if Self::M <= 8 {
                        return write!(f, "GF<{:#0X}>(value: 0x{:02X})", POLY, self.value);
                    } else if Self::M <= 12 {
                        return write!(f, "GF<{:#0X}>(value: 0x{:03X})", POLY, self.value);
                    } else {
                        return write!(f, "GF<{:#0X}>(value: 0x{:04X})", POLY, self.value);
                    }
                }
            }

            impl<const POLY: u128> fmt::Display for [<GF $type>]<POLY> {
                fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
                    if Self::M <= 4 {
                        return write!(f, "0x{:01X}", self.value);
                    } else if Self::M <= 8 {
                        return write!(f, "0x{:02X}", self.value);
                    } else if Self::M <= 12 {
                        return write!(f, "0x{:03X}", self.value);
                    } else {
                        return write!(f, "0x{:04X}", self.value);
                    }
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
                    let c = $type::clmul(self.value, other.value);

                    type Ctype = <$type as CarryLessMultiply>::OutType;

                    let (_, r) = Ctype::gf2_poly_div(c, POLY as Ctype);

                    Self {
                        value: r as $type,
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
}

#[repr(transparent)]
#[derive(Clone, Copy, Eq, PartialEq)]
pub struct GFu128<const POLY: u128> {
    pub value: u128,
}

impl<const POLY: u128> GaloisField for GFu128<POLY> {
    type StorageType = u128;

    const M: u128 = crate::calc_degree(POLY) as u128;
    const NUM_ELEM: u128 = 1 << Self::M;

    const ZERO: Self = Self { value: 0 };
    const ONE: Self = Self { value: 1 };

    fn inverse(&self) -> Self {
        if *self == Self::ZERO {
            panic!("Cannot take inverse of zero");
        }

        let mut temp_t: u128;
        let mut t: u128 = 0;
        let mut new_t: u128 = 1;

        let mut temp_r: u128;
        let mut r: u128 = POLY;
        let mut new_r: u128 = self.value as u128;

        let mut quotient: u128;
        let mut quo = Self::ZERO;
        let mut nt = Self::ZERO;

        while new_r != 0 {
            temp_r = r;
            r = new_r;
            (quotient, new_r) = u128::gf2_poly_div(temp_r, new_r);

            quo.value = quotient as Self::StorageType;
            nt.value = new_t as Self::StorageType;

            temp_t = t;
            t = new_t;
            new_t = temp_t ^ ((quo * nt).value as u128);
        }

        Self {
            value: t as Self::StorageType,
        }
    }

    fn new(value: u128) -> Self {
        Self { value: value }
    }

    fn validate(&self) -> bool {
        (self.value as u128) >= Self::NUM_ELEM
    }
}

impl<const POLY: u128> fmt::Debug for GFu128<POLY> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if Self::M <= 4 {
            return write!(f, "GF<{:#0X}>(value: 0x{:01X})", POLY, self.value);
        } else if Self::M <= 8 {
            return write!(f, "GF<{:#0X}>(value: 0x{:02X})", POLY, self.value);
        } else if Self::M <= 12 {
            return write!(f, "GF<{:#0X}>(value: 0x{:03X})", POLY, self.value);
        } else {
            return write!(f, "GF<{:#0X}>(value: 0x{:04X})", POLY, self.value);
        }
    }
}

impl<const POLY: u128> fmt::Display for GFu128<POLY> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if Self::M <= 4 {
            return write!(f, "0x{:01X}", self.value);
        } else if Self::M <= 8 {
            return write!(f, "0x{:02X}", self.value);
        } else if Self::M <= 12 {
            return write!(f, "0x{:03X}", self.value);
        } else {
            return write!(f, "0x{:04X}", self.value);
        }
    }
}

impl<const POLY: u128> Add<GFu128<POLY>> for GFu128<POLY> {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        Self {
            value: self.value ^ other.value,
        }
    }
}

impl<const POLY: u128> Sub<GFu128<POLY>> for GFu128<POLY> {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        Self {
            value: self.value ^ other.value,
        }
    }
}

impl<const POLY: u128> Mul<GFu128<POLY>> for GFu128<POLY> {
    type Output = Self;

    fn mul(self, other: Self) -> Self {
        let c = u128::clmul(self.value, other.value);

        type Ctype = <u128 as CarryLessMultiply>::OutType;

        let (_, r) = Ctype::gf2_poly_div(
            c,
            U256 {
                lower: POLY,
                upper: 0,
            },
        );

        Self { value: r.lower }
    }
}

impl<const POLY: u128> Div<GFu128<POLY>> for GFu128<POLY> {
    type Output = Self;

    fn div(self, other: Self) -> Self {
        self * other.inverse()
    }
}

assign_operator_impl! {
    u128: AddAssign: add_assign: +,
    u128: SubAssign: sub_assign: -,
    u128: MulAssign: mul_assign: *,
    u128: DivAssign: div_assign: /,
}
