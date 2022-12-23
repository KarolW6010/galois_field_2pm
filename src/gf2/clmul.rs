use paste::paste;
use core::ops::{BitAnd, BitOr, Not, BitXor, BitXorAssign, Shl, ShlAssign, Shr, ShrAssign};

use crate::calc_degree;

pub trait CarryLessMultiply {
    type InType;
    type OutType;

    fn clmul(a: Self::InType, b: Self::InType) -> Self::OutType;
}

macro_rules! clmul_impl{
    ($($in_type:ty: $out_type:ty,)*) => {
    $(
        impl CarryLessMultiply for $in_type {
            type InType = $in_type;
            type OutType = $out_type;

            fn clmul(a: Self::InType, b: Self::InType) -> Self::OutType {
                const NUM_BITS: usize = <$in_type>::BITS as usize;
                let mut c: $out_type = 0;
                for i in 0..NUM_BITS {
                    let set: bool = (b >> i) & 0x1 > 0;
                    if set {
                        c ^= (a as $out_type) << i;
                    }
                }
                c
            }
        }
    )*
    }
}

clmul_impl! {
    u8: u16,
    u16: u32,
    u32: u64,
    u64: u128,
}

#[derive(Debug, Clone, Copy)]
pub struct U256 {
    pub lower: u128,
    pub upper: u128,
}

impl PartialEq for U256 {
    fn eq(&self, other: &Self) -> bool {
        self.lower == other.lower && self.upper == other.upper
    }
}

impl Eq for U256 {}

impl U256 {
    fn new(val: u128, shift: u32) -> U256 {
        U256 {
            lower: val.checked_shl(shift).unwrap_or(0),
            upper: val.checked_shr(128 - shift).unwrap_or(0),
        }
    }
}

impl BitAnd<U256> for U256 {
    type Output = Self;

    fn bitand(self, rhs: U256) -> Self::Output {
        Self {
            lower: self.lower & rhs.lower,
            upper: self.upper & rhs.upper,
        }
    }
}

impl BitOr<U256> for U256 {
    type Output = Self;

    fn bitor(self, rhs: U256) -> Self::Output {
        Self {
            lower: self.lower | rhs.lower,
            upper: self.upper | rhs.upper,
        }
    }
}

impl Not for U256 {
    type Output = Self;

    fn not(self) -> Self::Output {
        Self {
            lower: !self.lower,
            upper: !self.upper,
        }
    }
}

impl BitXor<U256> for U256 {
    type Output = Self;

    fn bitxor(self, rhs: U256) -> Self::Output {
        Self {
            lower: self.lower ^ rhs.lower,
            upper: self.upper ^ rhs.upper,
        }
    }
}

impl BitXorAssign for U256 {
    fn bitxor_assign(&mut self, rhs: Self) {
        *self = *self ^ rhs
    }
}

impl Shl<usize> for U256 {
    type Output = Self;

    fn shl(self, shift: usize) -> Self::Output {
        match shift {
            0 => self,
            1..=127 => {
                let spill: u128 = self.lower >> (128 - shift);
                Self {
                    lower: self.lower << shift,
                    upper: (self.upper << shift) ^ spill,
                }
            }
            128..=255 => Self { lower: 0, upper: (self.lower << (shift - 128)) },
            _ => Self { lower: 0, upper: 0 }
        }
    }
}

impl ShlAssign<usize> for U256 {
    fn shl_assign(&mut self, rhs: usize) {
        *self = *self << rhs;
    }
}

impl Shr<usize> for U256 {
    type Output = Self;

    fn shr(self, shift: usize) -> Self::Output {
        match shift {
            0 => self,
            1..=127 => {
                let spill: u128 = self.upper << (128 - shift);
                Self {
                    lower: (self.lower >> shift) ^ spill,
                    upper: self.upper >> shift,
                }
            }
            128..=255 => Self { lower: (self.upper >> (shift - 128)), upper: 0 },
            _ => Self { lower: 0, upper: 0 }
        }
    }
}

impl ShrAssign<usize> for U256 {
    fn shr_assign(&mut self, rhs: usize) {
        *self = *self >> rhs;
    }
}

impl CarryLessMultiply for u128 {
    type InType = u128;
    type OutType = U256;

    fn clmul(a: Self::InType, b: Self::InType) -> Self::OutType {
        const NUM_BITS: usize = u128::BITS as usize;
        let mut c = U256 { lower: 0, upper: 0 };

        for i in 0..NUM_BITS {
            let set: bool = (b >> i) & 0x1 > 0;
            if set {
                c ^= U256::new(a, i as u32);
            }
        }
        c
    }
}

pub trait GF2PolyDiv {
    type Elem;

    // Perfroms dividend / divisor and returns (quotient, remainder)
    fn gf2_poly_div(dividend: Self::Elem, divisor: Self::Elem) -> (Self::Elem, Self::Elem);
}

macro_rules! gf2_poly_div_impl {
    ($($type:ty,)*) => {
    $(
        paste! {
            impl GF2PolyDiv for $type {
                type Elem = $type;

                fn gf2_poly_div(dividend: Self::Elem, divisor: Self::Elem) -> (Self::Elem, Self::Elem) {
                    if divisor == 0 {
                        panic!("Can not divide by zero");
                    }
                    if divisor == 1 {
                        return (dividend, 0);
                    }

                    let mut q: Self::Elem = 0;  // Quotient
                    let mut r: Self::Elem = 0;  // Remainder

                    let deg_divisor = calc_degree(divisor as u128); // Degree of divisor
                    let d = deg_divisor - 1;                        // Position of last register
                    let mask: Self::Elem = divisor & !((1 as Self::Elem) << deg_divisor);   // The divisor with the leading x^n term zeroed out

                    for cur_deg in (0..=calc_degree(dividend as u128)).rev() {
                        let bit: bool = ((r >> d) & 0x1) > 0;
                        if bit {
                            q = q << 1 | 0x1;
                        } else {
                            q <<= 1;
                        }
                        r = (r << 1) | ((dividend >> cur_deg) & 0x1);
                        if bit {
                            r ^= mask;
                        }
                    }

                    r &= ((1 as Self::Elem) << deg_divisor) - 1;

                    (q, r)
                }
            }
        }
    )*
    };
}

gf2_poly_div_impl! {
    u8,
    u16,
    u32,
    u64,
    u128,
}

fn calc_degree_256(x: U256) -> i32 {
    if x.upper > 0 {
        return (calc_degree(x.upper) as i32) + 128;
    }
    calc_degree(x.lower) as i32
}

impl GF2PolyDiv for U256 {
    type Elem = U256;
    
    fn gf2_poly_div(dividend: Self::Elem, divisor: Self::Elem) -> (Self::Elem, Self::Elem) {
        const ZERO: U256 = U256{upper: 0, lower: 0};
        const ONE: U256 = U256{upper: 0, lower: 1};
        if divisor == ZERO {
            panic!("Can not divide by zero");
        }
        if divisor == ONE {
            return (dividend, ZERO);
        }

        let mut q: Self::Elem = ZERO;  // Quotient
        let mut r: Self::Elem = ZERO;  // Remainder

        let deg_divisor = calc_degree_256(divisor);     // Degree of divisor
        let d = deg_divisor - 1;                        // Position of last register
        let mask: Self::Elem = divisor & !(U256::new(1, deg_divisor as u32));   // The divisor with the leading x^n term zeroed out

        for cur_deg in (0..=calc_degree_256(dividend)).rev() {
            let bit: bool = ((r >> d as usize).lower & 0x1) > 0;
            
            if bit {
                q = q << 1;
                q.lower |= 0x1;
            } else {
                q <<= 1;
            }

            r = (r << 1) | ((dividend >> (cur_deg as usize)) & ONE);
            if bit {
                r ^= mask;
            }
        }

        let mut mask: U256 = U256 { lower: 0, upper: 0 };
        if deg_divisor < 128 {
            mask.lower = (1 << deg_divisor) - 1;
        } else {
            mask.lower = u128::MAX;
            mask.upper = (1 << (deg_divisor - 128)) - 1;
        }

        r = r & mask;

        (q, r)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn clmul_check() {
        assert_eq!(u8::clmul(0xF, 0x81), 0x078F);
        assert_eq!(u16::clmul(0xF, 0x8001), 0x7800F);
        assert_eq!(
            u128::clmul(0xF, 0x8000_0000_0000_0000_0000_0000_0000_0001),
            U256 {
                lower: 0x8000_0000_0000_0000_0000_0000_0000_000F,
                upper: 0x7
            }
        );
    }
}
