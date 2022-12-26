pub trait CarryLessMultiply {
    fn clmul_low(&self, rhs: Self) -> Self;
    fn clmul_high(&self, rhs: Self) -> Self;
}

macro_rules! clmul_impl {
    ($($type:ty,)*) => {
    $(
        impl CarryLessMultiply for $type {
            fn clmul_low(&self, rhs: Self) -> Self {
                const NUM_BITS: u32 = <$type>::BITS;
                let mut output: Self = 0;

                for i in 0..NUM_BITS {
                    if ((rhs >> i) & 0x1) > 0 {
                        output ^= self << i;
                    }
                }

                output
            }

            fn clmul_high(&self, rhs: Self) -> Self {
                const NUM_BITS: u32 = <$type>::BITS;
                let mut output: Self = 0;

                for i in 1..NUM_BITS {
                    if ((rhs >> i) & 0x1) > 0 {
                        output ^= self >> (NUM_BITS - i);
                    }
                }

                output
            }
        }
    )*
    };
}

clmul_impl! {
    u8,
    u16,
    u32,
    u64,
    u128,
}
