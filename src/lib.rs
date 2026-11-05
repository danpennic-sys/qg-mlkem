//! # qg-mlkem
//!
//! Sealed, dependency-free, `no_std` ML-KEM ([FIPS 203]) with ML-KEM-512/768/1024.
//! Parameter sets, NTT tables, and reduction constants are produced by the in-repo compiler
//! `mlkemc` during the build; every derived constant is exhaustively verified there.
//!
//! * `#![forbid(unsafe_code)]`, no dependencies, no allocation.
//! * Secret-dependent paths are branch-free; decapsulation uses implicit rejection.
//! * [`Kem`] is sealed: only the three FIPS 203 sets exist.
//!
//! **Status:** draft, cross-checked only against `tools/ref.py`, not NIST ACVP. Run ACVP vectors and a
//! timing review before use. CNSA 2.0 requires ML-KEM-1024.
//!
//! [FIPS 203]: https://csrc.nist.gov/pubs/fips/203/final
#![no_std]
#![allow(clippy::needless_range_loop, clippy::too_many_arguments)]
#![forbid(unsafe_code)]
#![deny(missing_debug_implementations)]

#[cfg(feature = "std")]
extern crate std;

mod ct {
    //! Constant-time helpers (safe Rust; `black_box` as optimisation barrier, best effort).
    use core::hint::black_box;
    pub(crate) fn eq(a: &[u8], b: &[u8]) -> u32 {
        let mut d = 0u8;
        for (x, y) in a.iter().zip(b) {
            d |= x ^ y;
        }
        let x = black_box(d) as u32;
        ((x | x.wrapping_neg()) >> 31) ^ 1
    }
    pub(crate) fn select(out: &mut [u8], a: &[u8], b: &[u8], flag: u32) {
        let m = 0u8.wrapping_sub(black_box(flag) as u8);
        for i in 0..out.len() {
            out[i] = (a[i] & m) | (b[i] & !m);
        }
    }
    pub(crate) fn scrub(b: &mut [u8]) {
        b.fill(0);
        black_box(&*b);
    }
    pub(crate) fn scrub16(b: &mut [u16]) {
        b.fill(0);
        black_box(&*b);
    }
}
mod engine;
#[allow(clippy::all, dead_code)]
mod generated {
    include!(concat!(env!("OUT_DIR"), "/generated.rs"));
}
mod keccak;

pub use generated::{ml_kem_1024, ml_kem_512, ml_kem_768, MlKem1024, MlKem512, MlKem768};

mod sealed {
    pub trait Sealed {}
}

/// Source of cryptographically secure random bytes supplied by the caller.
pub trait Rng {
    fn fill_bytes(&mut self, dest: &mut [u8]);
}

/// Errors from key parsing.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Error {
    InvalidLength,
    InvalidEncapsulationKey,
    InvalidDecapsulationKey,
}

impl core::fmt::Display for Error {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.write_str(match self {
            Error::InvalidLength => "invalid length",
            Error::InvalidEncapsulationKey => "encapsulation key failed FIPS 203 modulus check",
            Error::InvalidDecapsulationKey => "decapsulation key failed hash check",
        })
    }
}
#[cfg(feature = "std")]
impl std::error::Error for Error {}

/// Sealed KEM interface implemented by [`MlKem512`], [`MlKem768`], [`MlKem1024`].
pub trait Kem: sealed::Sealed + Sized {
    const NAME: &'static str;
    const EK_LEN: usize;
    const DK_LEN: usize;
    const CT_LEN: usize;
    const SS_LEN: usize;
    type EncapsulationKey;
    type DecapsulationKey;
    type Ciphertext;
    type SharedSecret;
    fn generate<R: Rng>(rng: &mut R) -> (Self::EncapsulationKey, Self::DecapsulationKey);
    fn encapsulate<R: Rng>(
        ek: &Self::EncapsulationKey,
        rng: &mut R,
    ) -> (Self::Ciphertext, Self::SharedSecret);
    fn decapsulate(dk: &Self::DecapsulationKey, ct: &Self::Ciphertext) -> Self::SharedSecret;
}

/// OS randomness via `/dev/urandom` (unix, `std` feature). Panics if the read fails.
#[cfg(feature = "std")]
#[derive(Debug, Default)]
pub struct UrandomRng;
#[cfg(feature = "std")]
impl Rng for UrandomRng {
    fn fill_bytes(&mut self, dest: &mut [u8]) {
        use std::io::Read;
        std::fs::File::open("/dev/urandom")
            .and_then(|mut f| f.read_exact(dest))
            .expect("/dev/urandom unavailable");
    }
}
