include!("vectors.inc");
use qg_mlkem::{
    ml_kem_1024, ml_kem_512, ml_kem_768, Error, Kem, MlKem1024, MlKem512, MlKem768, Rng,
};
use std::fmt::Write;

struct TestRng(u64);
impl Rng for TestRng {
    fn fill_bytes(&mut self, dest: &mut [u8]) {
        for b in dest.iter_mut() {
            self.0 = self.0.wrapping_add(0x9e3779b97f4a7c15);
            let mut z = self.0;
            z = (z ^ (z >> 30)).wrapping_mul(0xbf58476d1ce4e5b9);
            z = (z ^ (z >> 27)).wrapping_mul(0x94d049bb133111eb);
            *b = (z ^ (z >> 31)) as u8;
        }
    }
}
fn hex(b: &[u8]) -> String {
    let mut s = String::with_capacity(b.len() * 2);
    for x in b {
        write!(s, "{x:02x}").unwrap();
    }
    s
}
fn range(start: u8, n: usize) -> Vec<u8> {
    (0..n as u8).map(|i| start + i).collect()
}

macro_rules! suite {
    ($m:ident, $vec:expr, $name:ident) => {
        mod $name {
            use super::*;
            use $m::*;

            #[test]
            fn sizes() {
                let (ek, dk) = keygen_from_seed(&[0u8; SEED_LEN]);
                let (ct, ss) = ek.encapsulate_from_seed(&[0u8; 32]);
                assert_eq!(
                    (
                        ek.as_bytes().len(),
                        dk.as_bytes().len(),
                        ct.as_bytes().len(),
                        ss.as_bytes().len()
                    ),
                    (EK_LEN, DK_LEN, CT_LEN, 32)
                );
            }

            #[test]
            fn roundtrip_many() {
                let mut rng = TestRng(7);
                for _ in 0..16 {
                    let (ek, dk) = keygen(&mut rng);
                    let (ct, ss) = ek.encapsulate(&mut rng);
                    assert!(dk.decapsulate(&ct) == ss);
                    assert_eq!(dk.encapsulation_key(), ek);
                }
            }

            #[test]
            fn cross_check_vs_python_reference() {
                let seed: [u8; 64] = range(0, 64).try_into().unwrap();
                let m: [u8; 32] = range(64, 32).try_into().unwrap();
                let v = &$vec;
                let (ek, dk) = keygen_from_seed(&seed);
                let (ct, ss) = ek.encapsulate_from_seed(&m);
                let k = K;
                assert_eq!(hex(&ek.as_bytes()[..32]), v.ek_head);
                assert_eq!(hex(&ek.as_bytes()[EK_LEN - 32..]), v.rho);
                assert_eq!(hex(&dk.as_bytes()[768 * k + 32..768 * k + 64]), v.dk_hash);
                assert_eq!(hex(&ct.as_bytes()[..32]), v.ct_head);
                assert_eq!(hex(&ct.as_bytes()[CT_LEN - 32..]), v.ct_tail);
                assert_eq!(hex(ss.as_bytes()), v.ss);
                // implicit rejection: flip ct[0] bit 0
                let mut bad = *ct.as_bytes();
                bad[0] ^= 1;
                let r = dk.decapsulate(&Ciphertext::from_bytes(&bad));
                assert_eq!(hex(r.as_bytes()), v.ss_reject);
                assert!(r != ss);
            }

            #[test]
            fn rejection_is_deterministic_and_input_dependent() {
                let mut rng = TestRng(99);
                let (ek, dk) = keygen(&mut rng);
                let (ct, _) = ek.encapsulate(&mut rng);
                for pos in [0usize, 1, CT_LEN / 2, CT_LEN - 1] {
                    let mut bad = *ct.as_bytes();
                    bad[pos] ^= 0x80;
                    let c = Ciphertext::from_bytes(&bad);
                    assert!(dk.decapsulate(&c) == dk.decapsulate(&c));
                }
                let mut b1 = *ct.as_bytes();
                b1[0] ^= 1;
                let mut b2 = *ct.as_bytes();
                b2[0] ^= 2;
                assert!(
                    dk.decapsulate(&Ciphertext::from_bytes(&b1))
                        != dk.decapsulate(&Ciphertext::from_bytes(&b2))
                );
            }

            #[test]
            fn ek_modulus_check() {
                let (ek, _) = keygen_from_seed(&[1u8; SEED_LEN]);
                assert!(EncapsulationKey::from_bytes(ek.as_bytes()).is_ok());
                let mut bad = *ek.as_bytes();
                bad[0] = 0xff;
                bad[1] |= 0x0f; // first 12-bit field = 4095 >= q
                assert_eq!(
                    EncapsulationKey::from_bytes(&bad),
                    Err(Error::InvalidEncapsulationKey)
                );
                assert_eq!(
                    EncapsulationKey::from_slice(&bad[1..]),
                    Err(Error::InvalidLength)
                );
            }

            #[test]
            fn dk_hash_check() {
                let (_, dk) = keygen_from_seed(&[2u8; SEED_LEN]);
                assert!(DecapsulationKey::from_bytes(dk.as_bytes()).is_ok());
                let mut bad = *dk.as_bytes();
                bad[384 * K + 5] ^= 1; // corrupt embedded ek
                assert!(matches!(
                    DecapsulationKey::from_bytes(&bad),
                    Err(Error::InvalidDecapsulationKey)
                ));
            }
        }
    };
}
suite!(ml_kem_512, V512, s512);
suite!(ml_kem_768, V768, s768);
suite!(ml_kem_1024, V1024, s1024);

#[test]
fn sealed_trait_is_generic() {
    fn rt<T: Kem>() -> usize
    where
        T::SharedSecret: PartialEq,
    {
        let mut rng = TestRng(1);
        let (ek, dk) = T::generate(&mut rng);
        let (ct, ss) = T::encapsulate(&ek, &mut rng);
        assert!(T::decapsulate(&dk, &ct) == ss);
        T::CT_LEN
    }
    assert_eq!(
        (rt::<MlKem512>(), rt::<MlKem768>(), rt::<MlKem1024>()),
        (768, 1088, 1568)
    );
    assert_eq!(
        (MlKem512::NAME, MlKem768::NAME, MlKem1024::NAME),
        ("ML-KEM-512", "ML-KEM-768", "ML-KEM-1024")
    );
}
