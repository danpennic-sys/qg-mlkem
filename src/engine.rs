//! FIPS 203 algorithms, written in the spec's own structure. Constants come from mlkemc output.
use crate::ct;
use crate::generated::{BARRETT_M, DIV_MUL, DIV_SHIFT, GAMMAS, N_INV, Q, ZETAS};
use crate::keccak::{self, shake128, shake256};

pub(crate) const N: usize = 256;
const MAX_K: usize = 4;
const MAX_CT: usize = 1568;
type Poly = [u16; N];
const ZERO: Poly = [0; N];

pub(crate) struct Params {
    pub k: usize,
    pub eta1: usize,
    pub eta2: usize,
    pub du: usize,
    pub dv: usize,
}

// ---- modular arithmetic (all branch-free) ------------------------------------------------------
#[inline(always)]
fn csub(x: u32) -> u32 {
    let t = x.wrapping_sub(Q);
    t.wrapping_add(0u32.wrapping_sub(t >> 31) & Q)
}
#[inline(always)]
fn reduce(a: u32) -> u32 {
    let t = ((a as u64 * BARRETT_M) >> 32) as u32;
    csub(a - t * Q)
}
#[inline(always)]
fn add(a: u32, b: u32) -> u32 {
    csub(a + b)
}
#[inline(always)]
fn sub(a: u32, b: u32) -> u32 {
    csub(a + Q - b)
}

fn compress(x: u32, d: usize) -> u32 {
    let v = ((x as u64) << d) + (Q as u64 >> 1);
    (((v * DIV_MUL) >> DIV_SHIFT) as u32) & ((1u32 << d) - 1)
}
fn decompress(y: u32, d: usize) -> u32 {
    (y * Q + (1u32 << (d - 1))) >> d
}

// ---- NTT (Algorithms 9-12) ---------------------------------------------------------------------
fn ntt(f: &mut Poly) {
    let (mut i, mut len) = (1usize, 128usize);
    while len >= 2 {
        let mut start = 0;
        while start < N {
            let z = ZETAS[i] as u32;
            i += 1;
            for j in start..start + len {
                let t = reduce(z * f[j + len] as u32);
                let a = f[j] as u32;
                f[j + len] = sub(a, t) as u16;
                f[j] = add(a, t) as u16;
            }
            start += 2 * len;
        }
        len >>= 1;
    }
}
fn inv_ntt(f: &mut Poly) {
    let (mut i, mut len) = (127usize, 2usize);
    while len <= 128 {
        let mut start = 0;
        while start < N {
            let z = ZETAS[i] as u32;
            i -= 1;
            for j in start..start + len {
                let t = f[j] as u32;
                let u = f[j + len] as u32;
                f[j] = add(t, u) as u16;
                f[j + len] = reduce(z * sub(u, t)) as u16;
            }
            start += 2 * len;
        }
        len <<= 1;
    }
    for c in f.iter_mut() {
        *c = reduce(*c as u32 * N_INV) as u16;
    }
}
/// h += f * g in the NTT domain (MultiplyNTTs with accumulation).
fn mul_acc(h: &mut Poly, f: &Poly, g: &Poly) {
    for i in 0..128 {
        let (a0, a1) = (f[2 * i] as u32, f[2 * i + 1] as u32);
        let (b0, b1) = (g[2 * i] as u32, g[2 * i + 1] as u32);
        let c0 = reduce(a0 * b0 + reduce(a1 * b1) * GAMMAS[i] as u32);
        let c1 = reduce(a0 * b1 + a1 * b0);
        h[2 * i] = add(h[2 * i] as u32, c0) as u16;
        h[2 * i + 1] = add(h[2 * i + 1] as u32, c1) as u16;
    }
}

// ---- sampling (Algorithms 7, 8) and encoding (Algorithms 5, 6) ---------------------------------
fn sample_ntt(rho: &[u8], i: u8, j: u8, a: &mut Poly) {
    let mut x = shake128();
    x.absorb(rho);
    x.absorb(&[j, i]);
    let (mut c, mut b) = (0usize, [0u8; 3]);
    while c < N {
        x.squeeze(&mut b);
        let d1 = b[0] as u32 | ((b[1] as u32 & 15) << 8);
        let d2 = (b[1] as u32 >> 4) | ((b[2] as u32) << 4);
        if d1 < Q {
            a[c] = d1 as u16;
            c += 1;
        }
        if d2 < Q && c < N {
            a[c] = d2 as u16;
            c += 1;
        }
    }
}
fn sample_cbd(eta: usize, b: &[u8], f: &mut Poly) {
    for i in 0..N {
        let (mut x, mut y) = (0u32, 0u32);
        for j in 0..eta {
            let n = 2 * i * eta + j;
            let m = n + eta;
            x += ((b[n >> 3] >> (n & 7)) & 1) as u32;
            y += ((b[m >> 3] >> (m & 7)) & 1) as u32;
        }
        f[i] = csub(x + Q - y) as u16;
    }
}
fn prf(eta: usize, s: &[u8], n: u8, out: &mut [u8; 192]) {
    shake256(&[s, &[n]], &mut out[..64 * eta]);
}

fn byte_encode(f: &Poly, d: usize, out: &mut [u8]) {
    debug_assert_eq!(out.len(), 32 * d);
    ct::scrub(out);
    let mut pos = 0usize;
    for &c in f.iter() {
        for j in 0..d {
            out[pos >> 3] |= (((c >> j) & 1) as u8) << (pos & 7);
            pos += 1;
        }
    }
}
fn byte_decode(b: &[u8], d: usize, f: &mut Poly) {
    debug_assert_eq!(b.len(), 32 * d);
    let mut pos = 0usize;
    for c in f.iter_mut() {
        let mut v = 0u32;
        for j in 0..d {
            v |= (((b[pos >> 3] >> (pos & 7)) & 1) as u32) << j;
            pos += 1;
        }
        *c = (if d == 12 { csub(v) } else { v }) as u16;
    }
}
fn compress_encode(f: &Poly, d: usize, out: &mut [u8]) {
    let mut t = ZERO;
    for i in 0..N {
        t[i] = compress(f[i] as u32, d) as u16;
    }
    byte_encode(&t, d, out);
}
fn decode_decompress(b: &[u8], d: usize, f: &mut Poly) {
    byte_decode(b, d, f);
    for c in f.iter_mut() {
        *c = decompress(*c as u32, d) as u16;
    }
}

fn gen_matrix(rho: &[u8], k: usize) -> [[Poly; MAX_K]; MAX_K] {
    let mut a = [[ZERO; MAX_K]; MAX_K];
    for i in 0..k {
        for j in 0..k {
            sample_ntt(rho, i as u8, j as u8, &mut a[i][j]);
        }
    }
    a
}

// ---- K-PKE (Algorithms 13-15) ------------------------------------------------------------------
fn kpke_keygen(p: &Params, d: &[u8; 32], ek: &mut [u8], dk_pke: &mut [u8]) {
    let k = p.k;
    let mut g_in = [0u8; 33];
    g_in[..32].copy_from_slice(d);
    g_in[32] = k as u8;
    let mut g = keccak::sha3_512(&[&g_in]);
    let (rho, sigma) = g.split_at(32);
    let a = gen_matrix(rho, k);
    let (mut s, mut e) = ([ZERO; MAX_K], [ZERO; MAX_K]);
    let (mut n, mut buf) = (0u8, [0u8; 192]);
    for i in 0..k {
        prf(p.eta1, sigma, n, &mut buf);
        n += 1;
        sample_cbd(p.eta1, &buf[..64 * p.eta1], &mut s[i]);
        ntt(&mut s[i]);
    }
    for i in 0..k {
        prf(p.eta1, sigma, n, &mut buf);
        n += 1;
        sample_cbd(p.eta1, &buf[..64 * p.eta1], &mut e[i]);
        ntt(&mut e[i]);
    }
    for i in 0..k {
        let mut t = ZERO;
        for j in 0..k {
            mul_acc(&mut t, &a[i][j], &s[j]);
        }
        for c in 0..N {
            t[c] = add(t[c] as u32, e[i][c] as u32) as u16;
        }
        byte_encode(&t, 12, &mut ek[384 * i..384 * (i + 1)]);
        byte_encode(&s[i], 12, &mut dk_pke[384 * i..384 * (i + 1)]);
    }
    ek[384 * k..384 * k + 32].copy_from_slice(rho);
    ct::scrub(&mut g_in);
    ct::scrub(&mut g);
    ct::scrub(&mut buf);
    for i in 0..MAX_K {
        ct::scrub16(&mut s[i]);
        ct::scrub16(&mut e[i]);
    }
}

fn kpke_encrypt(p: &Params, ek: &[u8], m: &[u8], r: &[u8], out: &mut [u8]) {
    let k = p.k;
    let mut t_hat = [ZERO; MAX_K];
    for i in 0..k {
        byte_decode(&ek[384 * i..384 * (i + 1)], 12, &mut t_hat[i]);
    }
    let a = gen_matrix(&ek[384 * k..384 * k + 32], k);
    let (mut y, mut e1) = ([ZERO; MAX_K], [ZERO; MAX_K]);
    let (mut n, mut buf) = (0u8, [0u8; 192]);
    for i in 0..k {
        prf(p.eta1, r, n, &mut buf);
        n += 1;
        sample_cbd(p.eta1, &buf[..64 * p.eta1], &mut y[i]);
        ntt(&mut y[i]);
    }
    for i in 0..k {
        prf(p.eta2, r, n, &mut buf);
        n += 1;
        sample_cbd(p.eta2, &buf[..64 * p.eta2], &mut e1[i]);
    }
    let mut e2 = ZERO;
    prf(p.eta2, r, n, &mut buf);
    sample_cbd(p.eta2, &buf[..64 * p.eta2], &mut e2);
    for i in 0..k {
        let mut u = ZERO;
        for j in 0..k {
            mul_acc(&mut u, &a[j][i], &y[j]);
        }
        inv_ntt(&mut u);
        for c in 0..N {
            u[c] = add(u[c] as u32, e1[i][c] as u32) as u16;
        }
        compress_encode(&u, p.du, &mut out[32 * p.du * i..32 * p.du * (i + 1)]);
    }
    let mut v = ZERO;
    for i in 0..k {
        mul_acc(&mut v, &t_hat[i], &y[i]);
    }
    inv_ntt(&mut v);
    let mut mu = ZERO;
    byte_decode(m, 1, &mut mu);
    for c in 0..N {
        let m_c = decompress(mu[c] as u32, 1);
        v[c] = add(add(v[c] as u32, e2[c] as u32), m_c) as u16;
    }
    compress_encode(&v, p.dv, &mut out[32 * p.du * k..32 * (p.du * k + p.dv)]);
    ct::scrub(&mut buf);
    for i in 0..MAX_K {
        ct::scrub16(&mut y[i]);
        ct::scrub16(&mut e1[i]);
    }
    ct::scrub16(&mut e2);
}

fn kpke_decrypt(p: &Params, dk_pke: &[u8], c: &[u8], m: &mut [u8; 32]) {
    let k = p.k;
    let mut w = ZERO;
    let (mut up, mut s_hat) = (ZERO, ZERO);
    for i in 0..k {
        decode_decompress(&c[32 * p.du * i..32 * p.du * (i + 1)], p.du, &mut up);
        ntt(&mut up);
        byte_decode(&dk_pke[384 * i..384 * (i + 1)], 12, &mut s_hat);
        mul_acc(&mut w, &s_hat, &up);
    }
    inv_ntt(&mut w);
    let mut v = ZERO;
    decode_decompress(&c[32 * p.du * k..32 * (p.du * k + p.dv)], p.dv, &mut v);
    for i in 0..N {
        w[i] = sub(v[i] as u32, w[i] as u32) as u16;
    }
    let mut bits = ZERO;
    for i in 0..N {
        bits[i] = compress(w[i] as u32, 1) as u16;
    }
    byte_encode(&bits, 1, m);
    ct::scrub16(&mut w);
    ct::scrub16(&mut s_hat);
    ct::scrub16(&mut bits);
}

// ---- ML-KEM (Algorithms 16-18) -----------------------------------------------------------------
pub(crate) fn keygen_internal(
    p: &Params,
    d: &[u8; 32],
    z: &[u8; 32],
    ek: &mut [u8],
    dk: &mut [u8],
) {
    let k = p.k;
    kpke_keygen(p, d, ek, &mut dk[..384 * k]);
    dk[384 * k..768 * k + 32].copy_from_slice(ek);
    let h = keccak::sha3_256(&[ek]);
    dk[768 * k + 32..768 * k + 64].copy_from_slice(&h);
    dk[768 * k + 64..768 * k + 96].copy_from_slice(z);
}

pub(crate) fn encaps_internal(
    p: &Params,
    ek: &[u8],
    m: &[u8; 32],
    out_ct: &mut [u8],
    ss: &mut [u8; 32],
) {
    let h = keccak::sha3_256(&[ek]);
    let mut g = keccak::sha3_512(&[m, &h]);
    ss.copy_from_slice(&g[..32]);
    kpke_encrypt(p, ek, m, &g[32..], out_ct);
    ct::scrub(&mut g);
}

pub(crate) fn decaps_internal(p: &Params, dk: &[u8], c: &[u8], ss: &mut [u8; 32]) {
    let k = p.k;
    let dk_pke = &dk[..384 * k];
    let ek_pke = &dk[384 * k..768 * k + 32];
    let h = &dk[768 * k + 32..768 * k + 64];
    let z = &dk[768 * k + 64..768 * k + 96];
    let mut m = [0u8; 32];
    kpke_decrypt(p, dk_pke, c, &mut m);
    let mut g = keccak::sha3_512(&[&m, h]);
    let mut k_bar = [0u8; 32];
    shake256(&[z, c], &mut k_bar);
    let mut c2 = [0u8; MAX_CT];
    let c2 = &mut c2[..c.len()];
    kpke_encrypt(p, ek_pke, &m, &g[32..], c2);
    let ok = ct::eq(c, c2);
    ct::select(ss, &g[..32], &k_bar, ok);
    ct::scrub(&mut m);
    ct::scrub(&mut g);
    ct::scrub(&mut k_bar);
}

/// FIPS 203 §7.2 input check: every 12-bit field of t-hat must already be reduced mod q.
pub(crate) fn ek_check(p: &Params, ek: &[u8]) -> bool {
    let mut f = ZERO;
    let mut re = [0u8; 384];
    for i in 0..p.k {
        let part = &ek[384 * i..384 * (i + 1)];
        byte_decode(part, 12, &mut f);
        byte_encode(&f, 12, &mut re);
        if re[..] != *part {
            return false;
        }
    }
    true
}
/// FIPS 203 §7.3 input check: H(ek) embedded in dk must match the embedded ek.
pub(crate) fn dk_check(p: &Params, dk: &[u8]) -> bool {
    let k = p.k;
    keccak::sha3_256(&[&dk[384 * k..768 * k + 32]])[..] == dk[768 * k + 32..768 * k + 64]
}
