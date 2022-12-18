mod clmul;
use clmul::{CarryLessMultiply, U256};

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
