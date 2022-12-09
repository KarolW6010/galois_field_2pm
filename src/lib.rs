use std::ops::{Add, AddAssign,
               Sub, SubAssign,
               Mul, MulAssign,
               Div, DivAssign,
               BitXor, BitXorAssign,
               BitAnd};

trait GFStorageType {}
trait GFType {}
trait GaloisField<T: GFType> {
    type Elem;
    fn zero() -> Self::Elem;
    fn one() -> Self::Elem;
}

struct GF<T: GFStorageType, const POLY: u128> {
    value: T,
}

impl<T: GFStorageType + BitXor<Output = T>, const POLY: u128> Add<GF<T, POLY>> for GF<T, POLY> {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        Self {
            value: self.value ^ other.value,
        }
    }
}

impl<T: GFStorageType, const POLY: u128> GFType for GF<T, POLY> {}

macro_rules! implement_traits_for_types {
    ($($type:ty,)*) => {
    $(
        impl GFStorageType for $type {}

        impl<const POLY: u128> GaloisField<GF<$type, POLY>> for GF<$type, POLY> {
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
    fn compile_test() {
        let gf0 = GF::<u32, 0x7>{value: 1};
        let gf0 = GF::<u32, 7>::zero();
        let gf1 = GF::<u32, 0x7>{value: 2};
        let gf2 = gf0 + gf1;
        assert_eq!(gf2.value, 2);
    }
}
