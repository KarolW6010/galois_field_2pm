use std::cmp;
use std::fmt;
use std::ops::{Shr, BitXor, Add, AddAssign, Sub, SubAssign, Mul, MulAssign};

const fn calc_degree(x: u128) -> i8 {
    if x > 0 {
        return calc_degree(x >> 1) + 1;
    } else {
        return -1;
    }
}

#[derive(Clone, Copy)]
pub struct GF<T: crate::GFStorageType, const POLY: u128> {
    pub value: T,
}

impl<T: crate::GFStorageType, const POLY: u128> GF<T, POLY> {
    const M: i8 = calc_degree(POLY);
}

impl<T: crate::GFStorageType, const POLY: u128> fmt::Debug for GF<T, POLY> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "GF<{:#0X}>(value: {:#0X})", POLY, self.value)
    }
}

impl<T: crate::GFStorageType, const POLY: u128> fmt::Display for GF<T, POLY> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "({:#0X})", self.value)
    }
}

impl<T: crate::GFStorageType, const POLY: u128> cmp::PartialEq for GF<T, POLY> {
    fn eq(&self, other: &Self) -> bool {
        self.value == other.value
    }
}

impl<T: crate::GFStorageType, const POLY: u128> cmp::Eq for GF<T, POLY> {}

impl<T: crate::GFStorageType + BitXor<Output = T>, const POLY: u128> Add<GF<T, POLY>> for GF<T, POLY> {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        Self {
            value: self.value ^ other.value,
        }
    }
}

impl<T: crate::GFStorageType + BitXor<Output = T>, const POLY: u128> Sub<GF<T, POLY>> for GF<T, POLY> {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        Self {
            value: self.value ^ other.value,
        }
    }
}

//impl<T: crate::GFStorageType + BitXor<Output = T> + Shr<Output = T>, const POLY: u128> Mul<GF<T, POLY>> for GF<T, POLY> {
//    type Output = Self;
//
//    fn mul(self, other: Self) -> Self {
//        let mut res: T = self.value ^ other.value;
//        let mut count: usize = 1;
//        while (count as usize) < (GF::<T, POLY>::M as usize) {
//            res ^= self.value ^ (other.value.shr(count));
//            count += 1;
//        }
//
//        Self {
//            value: res,
//        }
//    }
//}

impl<T: crate::GFStorageType + BitXor<Output = T>, const POLY: u128> AddAssign for GF<T, POLY> {
    fn add_assign(&mut self, other: Self) {
        *self = *self + other;
    }
}

impl<T: crate::GFStorageType + BitXor<Output = T>, const POLY: u128> SubAssign for GF<T, POLY> {
    fn sub_assign(&mut self, other: Self) {
        *self = *self - other;
    }
}

impl<T: crate::GFStorageType + BitXor<Output = T>, const POLY: u128> MulAssign for GF<T, POLY> {
    fn mul_assign(&mut self, other: Self) {
        *self = *self * other;
    }
}


macro_rules! implement_traits_for_types {
    ($($type:ty,)*) => {
    $(
        impl<const POLY: u128> Mul<GF<$type, POLY>> for GF<$type, POLY> {
            type Output = Self;
            
            fn mul(self, other: Self) -> Self {
                let mut res: $type = 0;
                for i in 0..GF::<$type, POLY>::M {
                    res ^= self.value ^ (other.value.shr(i));
                }

                Self { 
                    value: res,
                }
            }
        }

        impl<const POLY: u128> crate::GaloisField<GF<$type, POLY>> for GF<$type, POLY> {
            type Elem = GF<$type, POLY>;
            fn zero() -> Self::Elem { Self::Elem { value: 0 } }
            fn one() -> Self::Elem { Self::Elem { value: 1 } }
        }
    )*
    }
}

implement_traits_for_types! {
    u8,
    u16,
    u32,
    u64,
    u128,
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
