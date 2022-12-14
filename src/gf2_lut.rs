use paste::paste;
use std::fmt;
use std::ops::{Add, Sub, Mul, Div, AddAssign, SubAssign, MulAssign, DivAssign};

#[allow(dead_code)]
struct TablesU8 {
    exp_tbl: [u8; 1 << u8::BITS],
    log_tbl: [isize; 1 << u8::BITS],
}

#[allow(dead_code)]
const fn calc_degree(x: u128) -> i8 {
    if x > 0 {
        return calc_degree(x >> 1) + 1;
    } else {
        return -1;
    }
}

#[allow(dead_code)]
const fn generate_lut(poly: u128) -> TablesU8 {
    let m: usize = calc_degree(poly) as usize;
    let num_elems: usize = (1 << m) - 1;
    let feedback_mask: u8 = 1 << (m - 1);

    let mut exp_tbl: [u8; 1 << u8::BITS] = [0; 1 << u8::BITS];
    let mut log_tbl: [isize; 1 << u8::BITS] = [0; 1 << u8::BITS];
    
    log_tbl[0] = -1;
    let mut value: u8 = 1;
    let mut i: usize = 0;
    while i < num_elems {
        exp_tbl[i] = value;
        log_tbl[exp_tbl[i] as usize] = i as isize;
    
        let leading_1: bool = value >= feedback_mask;
        value <<= 1;
        if leading_1 {
            value ^= poly as u8;
        }
        i += 1;
    }
    
    TablesU8 {
        exp_tbl: exp_tbl,
        log_tbl: log_tbl,
    }
}

#[derive(Clone, Copy)]
pub struct GFu8<const POLY: u128> {
    value: u8,
}

#[allow(dead_code)]
const fn zero<const POLY: u128>() -> GFu8::<POLY> {
    GFu8::<POLY> {value: 0}
}

#[allow(dead_code)]
const fn one<const POLY: u128>() -> GFu8::<POLY> {
    GFu8::<POLY> {value: 1}
}

#[allow(dead_code)]
const fn alpha<const POLY: u128>() -> GFu8::<POLY> {
    GFu8::<POLY> {value: 2}
}

impl<const POLY: u128> GFu8<POLY> {
    const TABLES: TablesU8 = generate_lut(POLY);
    const M: usize = calc_degree(POLY) as usize;
    const NUM_ELEM: u128 = 1 << Self::M;
    const DEGREE_MOD: isize = (Self::NUM_ELEM as isize) - 1;

    const ZERO: Self = zero::<POLY>();
    const ONE: Self = one::<POLY>();
    const ALPHA: GFu8<POLY> = alpha::<POLY>();

    fn validate_value(&self) -> bool {
        (self.value as u128) >= Self::NUM_ELEM
    }

    fn new(value: u8) -> Self {
        Self {value: value}
    }

    fn new_validate(value: u8) -> Result<GFu8<POLY>, String> {
        let x = Self {value: value};
        if x.validate_value() {
            return Ok(x);
        } else {
            return Err(format!("GF(2^{}) represents values in the range [0x00 - 0x{:02X}]. The value 0x{:02X} is outside this range.", Self::M, Self::DEGREE_MOD, value));
        }
    }

    fn alpha_pow(mut power: isize) -> Self {
        power %= Self::DEGREE_MOD;
        power += Self::DEGREE_MOD;
        power %= Self::DEGREE_MOD;
        Self { value: Self::TABLES.exp_tbl[power as usize] }
    }

    fn log_alpha(&self) -> isize {
        Self::TABLES.log_tbl[self.value as usize]
    }

    fn inverse(&self) -> Self {
        if *self == Self::ZERO {
            panic!("Can not take inverse of zero");
        } else {
            return Self::alpha_pow(Self::DEGREE_MOD - Self::log_alpha(self));
        }
    }

    #[allow(dead_code)]
    fn validate() {
        if u8::BITS < (Self::M as u32) {
            panic!("GFu8 is not able to represent all elements in GF(2^{}) generated by polynomial {:#0X}", Self::M, POLY);
        }

        // Check if the polynomial is irreducible
        if POLY == 0x3 {
            return;
        }

        if POLY % 2 == 0 {
            panic!("{:#0X} is not irreducible. 0 is a root.", POLY);
        } else if POLY.count_ones() % 2 == 0 {
            panic!("{:#0X} is not irreducible. 1 is a root.", POLY);
        } else {
            for i in 2..Self::NUM_ELEM {
                if Self::TABLES.exp_tbl[i as usize] == 1 {
                    panic!("{:#0X} is irreducible but not primitive. Can not use LUT version.", POLY);
                }
            }
        }
    }
}

impl<const POLY: u128> PartialEq for GFu8<POLY> {
    fn eq(&self, other: &Self) -> bool {
        self.value == other.value
    }
}

impl<const POLY: u128> fmt::Debug for GFu8<POLY> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "GF<{:#0X}>(value: 0x{:02X})", POLY, self.value)
    }
}

impl<const POLY: u128> fmt::Display for GFu8<POLY> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "0x{:02X}", self.value)
    }
}

impl<const POLY: u128> Eq for GFu8<POLY> {}

impl<const POLY: u128> Add<GFu8<POLY>> for GFu8<POLY> {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        Self {
            value: self.value ^ other.value,
        }
    }
}

impl<const POLY: u128> Sub<GFu8<POLY>> for GFu8<POLY> {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        Self {
            value: self.value ^ other.value,
        }
    }
}

impl<const POLY: u128> Mul<GFu8<POLY>> for GFu8<POLY> {
    type Output = Self;

    fn mul(self, other: Self) -> Self {
        if (self == Self::ZERO) || (other == Self::ZERO) {
            return Self::ZERO;
        } else {
            return Self::alpha_pow(self.log_alpha() + other.log_alpha());
        }
    }
}

impl<const POLY: u128> Div<GFu8<POLY>> for GFu8<POLY> {
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

macro_rules! assign_operator_impl {
    ($($trait_name:ident: $trait_fn:ident: $op:tt,)*) => {
    $(
        impl<const POLY: u128> $trait_name for GFu8<POLY> {
            fn $trait_fn(&mut self, other: Self) {
                *self = *self $op other;
            }
        }
    )*
    }
}

assign_operator_impl! {
    AddAssign: add_assign: +,
    SubAssign: sub_assign: -,
    MulAssign: mul_assign: *,
    DivAssign: div_assign: /,
}

impl<const POLY: u128> crate::GaloisField for GFu8<POLY> {}

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

    macro_rules! associative_test {
        ($type:ty, $op:tt) => {
            for i in 0..GF::NUM_ELEM {
                for j in 0..GF::NUM_ELEM {
                    for k in 0..GF::NUM_ELEM {
                        let a = GF::new(i as $type);
                        let b = GF::new(j as $type);
                        let c = GF::new(k as $type);
                        assert_eq!(a $op (b $op c), (a $op b) $op c);
                    }
                }
            }
        }
    }

    macro_rules! commutative_test {
        ($type:ty, $op:tt) => {
            for i in 0..GF::NUM_ELEM {
                for j in 0..GF::NUM_ELEM {
                    let a = GF::new(i as $type);
                    let b = GF::new(j as $type);
                    assert_eq!(a $op b, b $op a);
                }
            }
        }
    }

    macro_rules! identity_test {
        ($type:ty, $op:tt, $identity:tt) => {
            for i in 0..GF::NUM_ELEM {
                let a = GF::new(i as $type);
                assert_eq!(a $op GF::$identity, a);
            }
        }
    }

    macro_rules! inverse_addition_test {
        ($type:ty) => {
            for i in 0..GF::NUM_ELEM {
                let a = GF::new(i as $type);
                assert_eq!(a + a, GF::ZERO);
            }
        }
    }

    macro_rules! inverse_multiplication_test {
        ($type:ty) => {
            for i in 1..GF::NUM_ELEM {
                let a = GF::new(i as $type);
                assert_eq!(a * a.inverse(), GF::ONE);
                assert_eq!(a / a, GF::ONE);
            }
        }
    }

    macro_rules! addition_test {
        ($type:ty) => {
            associative_test!($type, +);
            commutative_test!($type, +);
            identity_test!($type, +, ZERO);
            inverse_addition_test!($type);
        }
    }

    macro_rules! multiplication_test {
        ($type:ty) => {
            associative_test!($type, *);
            commutative_test!($type, *);
            identity_test!($type, *, ONE);
            inverse_multiplication_test!($type);
        }
    }

    macro_rules! distributive_test {
        ($type:ty) => {
            for i in 0..GF::NUM_ELEM {
                for j in 0..GF::NUM_ELEM {
                    for k in 0..GF::NUM_ELEM {
                        let a = GF::new(i as $type);
                        let b = GF::new(j as $type);
                        let c = GF::new(k as $type);
                        assert_eq!(a * (b + c), a * b + a * c);
                    }
                }
            }
        }
    }

    macro_rules! field_test_u8 {
        ($($name:ident: $poly:expr,)*) => {
        $(
            #[test]
            fn $name() {
                type GF = GFu8<$poly>;
                GF::validate();
                addition_test!(u8);
                multiplication_test!(u8);
                distributive_test!(u8);
            }

            paste! {
                #[test]
                #[should_panic]
                fn [<div_by_0_ $name>]() {
                    type GF = GFu8<$poly>;
                    let _ = GF::ONE / GF::ZERO;
                }

                #[test]
                #[should_panic]
                fn [<inverse_of_0_ $name>]() {
                    type GF = GFu8<$poly>;
                    let _ = GF::ZERO.inverse();
                }
            }
        )*
        }
    }

    field_test_u8! {
        gf_1: 0x3,
        gf_2: 0x7,
        gf_3: 0xB,
        gf_4: 0x13,
        gf_5: 0x25,
        gf_6: 0x43,
        gf_7: 0x83,
        gf_8: 0x11D,
    }

    macro_rules! check_not_primitive {
        ($($poly:expr,)*) => {
        $(
            paste! {
                #[test]
                #[should_panic]
                fn [<not_primitive_ $poly>]() {
                    GFu8::<$poly>::validate();
                }
            }
        )*
        }
    }

    check_not_primitive! {
        0x4,    // Zero is a root
        0x5,    // One is a root
        0x15,   // Irreducible but not primitive    
    }
}

