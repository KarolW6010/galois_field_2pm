pub trait GF2PolyDiv {
    type Elem;

    // This is used for inverting an element in GF(2^M)
    fn gf2_poly_div(dividend: Self::Elem, divisor: Self::Elem) -> (Self::Elem, Self::Elem);

    // This is used for inverting an element in GF(2^M). Used to divide POLY by lower degree value.
    // Divisor must be > 1
    fn gf2_poly_div_poly(dividend: u128, divisor: Self::Elem) -> (Self::Elem, Self::Elem);

    // Take the modulus of (hi | lo) by POLY
    fn gf2_poly_mod(hi: Self::Elem, lo: Self::Elem, poly: u128) -> Self::Elem;
}

macro_rules! gf2_poly_div_impl {
    ($($type:ty,)*) => {
    $(
        impl GF2PolyDiv for $type {
            type Elem = $type;

            fn gf2_poly_div(dividend: Self::Elem, divisor: Self::Elem) -> (Self::Elem, Self::Elem) {
                if divisor == 1 {
                    return (dividend, 0);
                }

                let mut quotient: Self::Elem = 0;
                let mut remainder: Self::Elem = 0;

                let deg_divisor = crate::calc_degree(divisor as u128);
                let d = deg_divisor - 1;
                let mask: Self::Elem = divisor & !(1 << deg_divisor);

                for cur_deg in (0..=crate::calc_degree(dividend as u128)).rev() {
                    let bit: bool = ((remainder >> d) & 0x1) > 0;
                    if bit {
                        quotient = (quotient << 1) | 0x1;
                    } else {
                        quotient <<= 1;
                    }

                    remainder = (remainder << 1) | ((dividend >> cur_deg) & 0x1);
                    if bit {
                        remainder ^= mask;
                    }
                }

                let bit_limit: Self::Elem = (1 << deg_divisor) - 1;
                remainder &= bit_limit;
                (quotient, remainder)
            }

            fn gf2_poly_div_poly(dividend: u128, divisor: Self::Elem) -> (Self::Elem, Self::Elem) {
                let mut quotient: Self::Elem = 0;
                let mut remainder: Self::Elem = 0;

                let deg_divisor = crate::calc_degree(divisor as u128);
                let d = deg_divisor - 1;
                let mask: Self::Elem = divisor & !(1 << deg_divisor);

                for cur_deg in (0..=crate::calc_degree(dividend as u128)).rev() {
                    let bit: bool = ((remainder >> d) & 0x1) > 0;
                    if bit {
                        quotient = (quotient << 1) | 0x1;
                    } else {
                        quotient <<= 1;
                    }

                    remainder = (remainder << 1) | (((dividend >> cur_deg) as Self::Elem) & 0x1) ;
                    if bit {
                        remainder ^= mask;
                    }
                }

                let bit_limit: Self::Elem = (1 << deg_divisor) - 1;
                remainder &= bit_limit;
                (quotient, remainder)
            }

            fn gf2_poly_mod(hi: Self::Elem, lo: Self::Elem, poly: u128) -> Self::Elem {
                let deg_divisor: usize = crate::calc_degree(poly) as usize;
                let d = deg_divisor - 1;
                let mask: Self::Elem = (poly & !(1u128 << deg_divisor)) as Self::Elem;

                let mut remainder: Self::Elem = 0;

                let hi_deg = crate::calc_degree(hi as u128);
                let lo_deg = if hi_deg < 0 { crate::calc_degree(lo as u128) } else { (Self::Elem::BITS - 1) as i16};

                // Start with the hi portion
                for cur_deg in (0..=hi_deg).rev() {
                    let bit: bool = ((remainder >> d) & 0x1) > 0;
                    remainder = (remainder << 1) | ((hi >> cur_deg) & 0x1);
                    if bit {
                        remainder ^= mask;
                    }
                }

                // Now use the lo portion
                for cur_deg in (0..=lo_deg).rev() {
                    let bit: bool = ((remainder >> d) & 0x1) > 0;
                    remainder = (remainder << 1) | ((lo >> cur_deg) & 0x1);
                    if bit {
                        remainder ^= mask;
                    }
                }

                let bit_limit: Self::Elem = ((1u128 << deg_divisor) - 1) as Self::Elem;
                remainder &= bit_limit;

                remainder
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
