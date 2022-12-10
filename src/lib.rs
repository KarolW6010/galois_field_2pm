use std::fmt::{Debug, Display, UpperHex};
use std::ops::{Add, AddAssign,
               Sub, SubAssign,
               Mul, MulAssign,
               Div, DivAssign,
               BitXor, BitXorAssign,
               Shr};

pub mod gf2;

pub trait GFStorageType: Clone + Copy + Sized + Shr + BitXor + BitXorAssign + Eq + PartialEq + Debug + Display + UpperHex {}
impl GFStorageType for u8 {}
impl GFStorageType for u16 {}
impl GFStorageType for u32 {}
impl GFStorageType for u64 {}
impl GFStorageType for u128 {}

pub trait GaloisField<T: Eq + PartialEq + Clone + Copy +
                         Add + AddAssign + Sub + SubAssign +
                         Mul + MulAssign> {
    type Elem;
    fn zero() -> Self::Elem;
    fn one() -> Self::Elem;
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::gf2;

    #[test]
    fn compile_test() {
        type GF = gf2::GF<u8, 0x7>;
        //let gf0 = GF{value: 1};
        let gf0 = GF::zero();
        let mut gf1 = GF{value: 2};
        gf1 += gf1;
        let gf2 = gf0 + gf1;
        assert_eq!(gf2, gf1);
    }
}
