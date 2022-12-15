use paste::paste;
use std::fmt;
use std::ops::{Add, Sub, Mul, Div, AddAssign, SubAssign, MulAssign, DivAssign};

#[allow(dead_code)]
const fn calc_degree(x: u128) -> i8 {
    if x > 0 {
        return calc_degree(x >> 1) + 1;
    } else {
        return -1;
    }
}

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
            struct [<Tables $type:upper>] {
                exp_tbl: [$type; 1 << $type::BITS],
                log_tbl: [isize; 1 << $type::BITS],
            }

            const fn [<generate_lut_ $type:lower>](poly: u128) -> [<Tables $type:upper>] {
                let m: usize = calc_degree(poly) as usize;
                let num_elems: usize = (1 << m) - 1;
                let feedback_mask: $type = 1 << (m - 1);
            
                let mut exp_tbl: [$type; 1 << $type::BITS] = [0; 1 << $type::BITS];
                let mut log_tbl: [isize; 1 << $type::BITS] = [0; 1 << $type::BITS];
                
                log_tbl[0] = -1;
                let mut value: $type = 1;
                let mut i: usize = 0;
                while i < num_elems {
                    exp_tbl[i] = value;
                    log_tbl[exp_tbl[i] as usize] = i as isize;
                
                    let leading_1: bool = value >= feedback_mask;
                    value <<= 1;
                    if leading_1 {
                        value ^= poly as $type;
                    }
                    i += 1;
                }
                
                [<Tables $type:upper>] {
                    exp_tbl: exp_tbl,
                    log_tbl: log_tbl,
                }
            }

            #[derive(Clone, Copy)]
            pub struct [<GF $type>]<const POLY: u128> {
                value: $type,
            }

            #[allow(dead_code)]
            const fn [<zero_ $type>]<const POLY: u128>() -> [<GF $type>]::<POLY> {
                [<GF $type>]::<POLY> {value: 0}
            }

            #[allow(dead_code)]
            const fn [<one_ $type>]<const POLY: u128>() -> [<GF $type>]::<POLY> {
                [<GF $type>]::<POLY> {value: 1}
            }

            #[allow(dead_code)]
            const fn [<alpha_ $type>]<const POLY: u128>() -> [<GF $type>]::<POLY> {
                [<GF $type>]::<POLY> {value: 2}
            }
        
            impl<const POLY: u128> [<GF $type>]<POLY> {
                const TABLES: [<Tables $type:upper>] = [<generate_lut_ $type:lower>](POLY);
                const M: usize = calc_degree(POLY) as usize;
                const NUM_ELEM: u128 = 1 << Self::M;
                const DEGREE_MOD: isize = (Self::NUM_ELEM as isize) - 1;

                const ZERO: Self = [<zero_ $type>]::<POLY>();
                const ONE: Self = [<one_ $type>]::<POLY>();
                const ALPHA: Self = [<alpha_ $type>]::<POLY>();

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

                fn new(value: $type) -> Self {
                    Self {value: value}
                }

                fn validate() {
                    // Check if the polynomial is irreducible
                    if POLY == 0x3 {
                        return;
                    }
            
                    if POLY % 2 == 0 {
                        panic!("{:#0X} is not irreducible. 0 is a root.", POLY);
                    }
                    if POLY.count_ones() % 2 == 0 {
                        panic!("{:#0X} is not irreducible. 1 is a root.", POLY);
                    }

                    // Polynomial is irreducible, check if primitive
                    for i in 2..Self::NUM_ELEM {
                        if Self::TABLES.exp_tbl[i as usize] == 1 {
                            panic!("{:#0X} is irreducible but not primitive. Can not use LUT version.", POLY);
                        }
                    }
                }
            }

            impl<const POLY: u128> PartialEq for [<GF $type>]<POLY> {
                fn eq(&self, other: &Self) -> bool {
                    self.value == other.value
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

            impl<const POLY: u128> Eq for [<GF $type>]<POLY> {}

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
                    if (self == Self::ZERO) || (other == Self::ZERO) {
                        return Self::ZERO;
                    } else {
                        return Self::alpha_pow(self.log_alpha() + other.log_alpha());
                    }
                }
            }
            
            impl<const POLY: u128> Div<[<GF $type>]<POLY>> for [<GF $type>]<POLY> {
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

            assign_operator_impl! {
                $type: AddAssign: add_assign: +,
                $type: SubAssign: sub_assign: -,
                $type: MulAssign: mul_assign: *,
                $type: DivAssign: div_assign: /,
            }

            impl<const POLY: u128> crate::GaloisField for [<GF $type>]<POLY> {}
        }
    )*
    }
}

setup_gf!{
    u8,
    u16,
}

impl<const POLY: u128> GFu8<POLY> {
/*
    fn validate_value(&self) -> bool {
        (self.value as u128) >= Self::NUM_ELEM
    }

    fn new_validate(value: u8) -> Result<GFu8<POLY>, String> {
        let x = Self {value: value};
        if x.validate_value() {
            return Ok(x);
        } else {
            return Err(format!("GF(2^{}) represents values in the range [0x00 - 0x{:02X}]. The value 0x{:02X} is outside this range.", Self::M, Self::DEGREE_MOD, value));
        }
    }
*/
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

    macro_rules! field_test {
        ($($type:ty: $poly:expr,)*) => {
        $(
            paste! {
                #[test]
                fn [<gf_ $poly>]() {
                    type GF = [<GF $type>]<$poly>;
                    GF::validate();
                    addition_test!($type);
                    multiplication_test!($type);
                    distributive_test!($type);
                }

                #[test]
                #[should_panic]
                fn [<div_by_0_ $poly>]() {
                    type GF = [<GF $type>]<$poly>;
                    let _ = GF::ONE / GF::ZERO;
                }

                #[test]
                #[should_panic]
                fn [<inverse_of_0_ $poly>]() {
                    type GF = [<GF $type>]<$poly>;
                    let _ = GF::ZERO.inverse();
                }
            }
        )*
        }
    }

    field_test! {
        u16: 0x3,
        u16: 0x7,
        u16: 0xB,
        u16: 0x13,
        u8: 0x25,
        u8: 0x43,
        u8: 0x83,
        u8: 0x163,
    }

    macro_rules! check_not_primitive {
        ($($type:ty: $poly:expr,)*) => {
        $(
            paste! {
                #[test]
                #[should_panic]
                fn [<not_primitive_ $poly>]() {
                    [<GF $type>]::<$poly>::validate();
                }
            }
        )*
        }
    }

    check_not_primitive! {
        u8: 0x4,    // Zero is a root
        u8: 0x5,    // One is a root
        u8: 0x15,   // Irreducible but not primitive    
        u16: 0x203, // Irreducible but not primitive
    }
}

