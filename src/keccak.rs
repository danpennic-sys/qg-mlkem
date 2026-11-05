//! Keccak-f[1600], SHA3-256/512 and SHAKE128/256 (FIPS 202). Byte-oriented, small, auditable.
use core::hint::black_box;

const RC: [u64; 24] = [
    0x0000000000000001,
    0x0000000000008082,
    0x800000000000808a,
    0x8000000080008000,
    0x000000000000808b,
    0x0000000080000001,
    0x8000000080008081,
    0x8000000000008009,
    0x000000000000008a,
    0x0000000000000088,
    0x0000000080008009,
    0x000000008000000a,
    0x000000008000808b,
    0x800000000000008b,
    0x8000000000008089,
    0x8000000000008003,
    0x8000000000008002,
    0x8000000000000080,
    0x000000000000800a,
    0x800000008000000a,
    0x8000000080008081,
    0x8000000000008080,
    0x0000000080000001,
    0x8000000080008008,
];
const ROTC: [u32; 24] = [
    1, 3, 6, 10, 15, 21, 28, 36, 45, 55, 2, 14, 27, 41, 56, 8, 25, 43, 62, 18, 39, 61, 20, 44,
];
const PILN: [usize; 24] = [
    10, 7, 11, 17, 18, 3, 5, 16, 8, 21, 24, 4, 15, 23, 19, 13, 12, 2, 20, 14, 22, 9, 6, 1,
];

fn f1600(st: &mut [u64; 25]) {
    for rc in RC {
        let mut bc = [0u64; 5];
        for i in 0..5 {
            bc[i] = st[i] ^ st[i + 5] ^ st[i + 10] ^ st[i + 15] ^ st[i + 20];
        }
        for i in 0..5 {
            let t = bc[(i + 4) % 5] ^ bc[(i + 1) % 5].rotate_left(1);
            for j in (0..25).step_by(5) {
                st[j + i] ^= t;
            }
        }
        let mut t = st[1];
        for i in 0..24 {
            let j = PILN[i];
            let b = st[j];
            st[j] = t.rotate_left(ROTC[i]);
            t = b;
        }
        for j in (0..25).step_by(5) {
            bc.copy_from_slice(&st[j..j + 5]);
            for i in 0..5 {
                st[j + i] ^= !bc[(i + 1) % 5] & bc[(i + 2) % 5];
            }
        }
        st[0] ^= rc;
    }
}

pub(crate) struct Sponge {
    st: [u64; 25],
    rate: usize,
    pos: usize,
    ds: u8,
    squeezing: bool,
}

impl Sponge {
    pub fn new(rate: usize, ds: u8) -> Self {
        Sponge {
            st: [0; 25],
            rate,
            pos: 0,
            ds,
            squeezing: false,
        }
    }
    pub fn absorb(&mut self, data: &[u8]) {
        debug_assert!(!self.squeezing);
        for &b in data {
            self.st[self.pos / 8] ^= (b as u64) << (8 * (self.pos % 8));
            self.pos += 1;
            if self.pos == self.rate {
                f1600(&mut self.st);
                self.pos = 0;
            }
        }
    }
    fn finalize(&mut self) {
        self.st[self.pos / 8] ^= (self.ds as u64) << (8 * (self.pos % 8));
        let last = self.rate - 1;
        self.st[last / 8] ^= 0x80u64 << (8 * (last % 8));
        f1600(&mut self.st);
        self.pos = 0;
        self.squeezing = true;
    }
    pub fn squeeze(&mut self, out: &mut [u8]) {
        if !self.squeezing {
            self.finalize();
        }
        for o in out.iter_mut() {
            if self.pos == self.rate {
                f1600(&mut self.st);
                self.pos = 0;
            }
            *o = (self.st[self.pos / 8] >> (8 * (self.pos % 8))) as u8;
            self.pos += 1;
        }
    }
}
impl Drop for Sponge {
    fn drop(&mut self) {
        self.st = [0; 25];
        black_box(&self.st);
    }
}

pub(crate) const SHAKE128_RATE: usize = 168;
pub(crate) fn shake128() -> Sponge {
    Sponge::new(SHAKE128_RATE, 0x1f)
}

fn hash<const L: usize>(rate: usize, ds: u8, parts: &[&[u8]]) -> [u8; L] {
    let mut s = Sponge::new(rate, ds);
    for p in parts {
        s.absorb(p);
    }
    let mut o = [0u8; L];
    s.squeeze(&mut o);
    o
}
/// H in FIPS 203.
pub(crate) fn sha3_256(parts: &[&[u8]]) -> [u8; 32] {
    hash::<32>(136, 0x06, parts)
}
/// G in FIPS 203.
pub(crate) fn sha3_512(parts: &[&[u8]]) -> [u8; 64] {
    hash::<64>(72, 0x06, parts)
}
/// PRF / J in FIPS 203.
pub(crate) fn shake256(parts: &[&[u8]], out: &mut [u8]) {
    let mut s = Sponge::new(136, 0x1f);
    for p in parts {
        s.absorb(p);
    }
    s.squeeze(out);
}

#[cfg(test)]
mod tests {
    use super::*;
    fn hex(b: &[u8]) -> alloc_free::Hex {
        alloc_free::Hex::from(b)
    }
    mod alloc_free {
        pub struct Hex(pub [u8; 128], pub usize);
        impl Hex {
            pub fn from(b: &[u8]) -> Hex {
                let mut h = Hex([0; 128], 0);
                for x in b {
                    for n in [x >> 4, x & 15] {
                        h.0[h.1] = b"0123456789abcdef"[n as usize];
                        h.1 += 1;
                    }
                }
                h
            }
            pub fn s(&self) -> &str {
                core::str::from_utf8(&self.0[..self.1]).unwrap()
            }
        }
    }
    #[test]
    fn sha3_256_empty() {
        assert_eq!(
            hex(&sha3_256(&[])).s(),
            "a7ffc6f8bf1ed76651c14756a061d662f580ff4de43b49fa82d80a4b80f8434a"
        );
    }
    #[test]
    fn sha3_512_empty() {
        assert_eq!(hex(&sha3_512(&[])).s(), "a69f73cca23a9ac5c8b567dc185a756e97c982164fe25859e0d1dcc1475c80a615b2123af1f5f94c11e3e9402c3ac558f500199d95b6d3e301758586281dcd26");
    }
    #[test]
    fn shake128_empty() {
        let mut o = [0u8; 32];
        let mut s = shake128();
        s.squeeze(&mut o);
        assert_eq!(
            hex(&o).s(),
            "7f9c2ba4e88f827d616045507605853ed73b8093f6efbc88eb1a6eacfa66ef26"
        );
    }
    #[test]
    fn shake256_empty() {
        let mut o = [0u8; 32];
        shake256(&[], &mut o);
        assert_eq!(
            hex(&o).s(),
            "46b9dd2b0ba88d13233b3feb743eeb243fcd52ea62b81b82b50c27646ed5762f"
        );
    }
    #[test]
    fn multi_block_and_split_absorb() {
        // 200 bytes crosses the SHA3-256 rate (136); splitting the input must not change the digest.
        let m = [0xa3u8; 200];
        assert_eq!(sha3_256(&[&m]), sha3_256(&[&m[..77], &m[77..]]));
    }
}
