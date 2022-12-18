use std::ops::{BitXor, BitXorAssign};

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
