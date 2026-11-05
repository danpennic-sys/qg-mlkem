#!/usr/bin/env python3
"""Independent FIPS 203 reference (spec-style, slow, hashlib only).
Used to (a) cross-check the Rust crate and (b) generate tests/vectors.rs.
NOT NIST ACVP output: run the official ACVP KATs before any real use."""
import hashlib, sys
Q = 3329; ZETA = 17
def brv7(i): return int(f"{i:07b}"[::-1], 2)
ZETAS = [pow(ZETA, brv7(i), Q) for i in range(128)]
G = lambda b: hashlib.sha3_512(b).digest()
H = lambda b: hashlib.sha3_256(b).digest()
J = lambda b: hashlib.shake_256(b).digest(32)
PRF = lambda eta, s, b: hashlib.shake_256(s + bytes([b])).digest(64 * eta)

def byte_encode(f, d):
    bits = []
    for c in f: bits += [(c >> j) & 1 for j in range(d)]
    return bytes(sum(bits[8*i+j] << j for j in range(8)) for i in range(len(bits)//8))
def byte_decode(b, d):
    bits = [(b[i//8] >> (i % 8)) & 1 for i in range(8*len(b))]
    m = Q if d == 12 else 1 << d
    return [sum(bits[i*d+j] << j for j in range(d)) % m for i in range(256)]
def compress(x, d): return ((x << d) * 2 + Q) // (2 * Q) % (1 << d)
def decompress(y, d): return (Q * y * 2 + (1 << d)) // (1 << (d + 1))

def sample_ntt(b):
    stream = hashlib.shake_128(b).digest(3 * 1024); a = []; i = 0
    while len(a) < 256:
        c0, c1, c2 = stream[i:i+3]; i += 3
        d1 = c0 + 256 * (c1 % 16); d2 = c1 // 16 + 16 * c2
        if d1 < Q: a.append(d1)
        if d2 < Q and len(a) < 256: a.append(d2)
    return a
def cbd(b, eta):
    bits = [(b[i//8] >> (i % 8)) & 1 for i in range(8*len(b))]
    return [(sum(bits[2*i*eta+j] for j in range(eta)) - sum(bits[2*i*eta+eta+j] for j in range(eta))) % Q for i in range(256)]
def ntt(f):
    f = f[:]; i = 1; ln = 128
    while ln >= 2:
        for s in range(0, 256, 2*ln):
            z = ZETAS[i]; i += 1
            for j in range(s, s+ln):
                t = z * f[j+ln] % Q; f[j+ln] = (f[j] - t) % Q; f[j] = (f[j] + t) % Q
        ln //= 2
    return f
def intt(f):
    f = f[:]; i = 127; ln = 2
    while ln <= 128:
        for s in range(0, 256, 2*ln):
            z = ZETAS[i]; i -= 1
            for j in range(s, s+ln):
                t = f[j]; f[j] = (t + f[j+ln]) % Q; f[j+ln] = z * (f[j+ln] - t) % Q
        ln *= 2
    return [x * 3303 % Q for x in f]
def mul(f, g):
    h = [0]*256
    for i in range(128):
        gam = pow(ZETA, 2*brv7(i)+1, Q)
        a0, a1, b0, b1 = f[2*i], f[2*i+1], g[2*i], g[2*i+1]
        h[2*i] = (a0*b0 + a1*b1*gam) % Q; h[2*i+1] = (a0*b1 + a1*b0) % Q
    return h
def vadd(a, b): return [(x+y) % Q for x, y in zip(a, b)]

P = {512: (2,3,2,10,4), 768: (3,2,2,10,4), 1024: (4,2,2,11,5)}
def matrix(rho, k): return [[sample_ntt(rho + bytes([j, i])) for j in range(k)] for i in range(k)]

def kpke_keygen(d, p):
    k, e1, e2, du, dv = p
    g = G(d + bytes([k])); rho, sigma = g[:32], g[32:]; A = matrix(rho, k); n = 0
    s = []; e = []
    for _ in range(k): s.append(ntt(cbd(PRF(e1, sigma, n), e1))); n += 1
    for _ in range(k): e.append(ntt(cbd(PRF(e1, sigma, n), e1))); n += 1
    t = []
    for i in range(k):
        acc = [0]*256
        for j in range(k): acc = vadd(acc, mul(A[i][j], s[j]))
        t.append(vadd(acc, e[i]))
    return b"".join(byte_encode(x, 12) for x in t) + rho, b"".join(byte_encode(x, 12) for x in s)
def kpke_enc(ek, m, r, p):
    k, e1, e2, du, dv = p
    t = [byte_decode(ek[384*i:384*i+384], 12) for i in range(k)]; rho = ek[384*k:]
    A = matrix(rho, k); n = 0; y = []; E1 = []
    for _ in range(k): y.append(ntt(cbd(PRF(e1, r, n), e1))); n += 1
    for _ in range(k): E1.append(cbd(PRF(e2, r, n), e2)); n += 1
    E2 = cbd(PRF(e2, r, n), e2)
    u = []
    for i in range(k):
        acc = [0]*256
        for j in range(k): acc = vadd(acc, mul(A[j][i], y[j]))
        u.append(vadd(intt(acc), E1[i]))
    mu = [decompress(x, 1) for x in byte_decode(m, 1)]
    acc = [0]*256
    for i in range(k): acc = vadd(acc, mul(t[i], y[i]))
    v = vadd(vadd(intt(acc), E2), mu)
    return b"".join(byte_encode([compress(x, du) for x in ui], du) for ui in u) + byte_encode([compress(x, dv) for x in v], dv)
def kpke_dec(dk, c, p):
    k, e1, e2, du, dv = p
    u = [[decompress(x, du) for x in byte_decode(c[32*du*i:32*du*(i+1)], du)] for i in range(k)]
    v = [decompress(x, dv) for x in byte_decode(c[32*du*k:], dv)]
    s = [byte_decode(dk[384*i:384*i+384], 12) for i in range(k)]
    acc = [0]*256
    for i in range(k): acc = vadd(acc, mul(s[i], ntt(u[i])))
    w = [(a - b) % Q for a, b in zip(v, intt(acc))]
    return byte_encode([compress(x, 1) for x in w], 1)
def keygen(d, z, p):
    ek, dkp = kpke_keygen(d, p); return ek, dkp + ek + H(ek) + z
def encaps(ek, m, p):
    g = G(m + H(ek)); return g[:32], kpke_enc(ek, m, g[32:], p)
def decaps(dk, c, p):
    k = p[0]; dkp = dk[:384*k]; ek = dk[384*k:768*k+32]; h = dk[768*k+32:768*k+64]; z = dk[768*k+64:]
    m = kpke_dec(dkp, c, p); g = G(m + h); K, r = g[:32], g[32:]
    return K if kpke_enc(ek, m, r, p) == c else J(z + c)

def vectors():
    out = {}
    for lvl, p in P.items():
        d = bytes(range(32)); z = bytes(range(32, 64)); m = bytes(range(64, 96))
        ek, dk = keygen(d, z, p); K, c = encaps(ek, m, p)
        assert decaps(dk, c, p) == K
        bad = bytearray(c); bad[0] ^= 1; Kr = decaps(dk, bytes(bad), p)
        assert Kr != K and Kr == J(z + bytes(bad))
        k = p[0]
        out[lvl] = dict(sizes=(len(ek), len(dk), len(c)), ek_head=ek[:32], rho=ek[-32:], dk_hash=dk[768*k+32:768*k+64],
                        ct_head=c[:32], ct_tail=c[-32:], ss=K, ss_reject=Kr)
    return out
if __name__ == "__main__":
    v = vectors()
    for lvl, d in v.items(): print(lvl, d["sizes"])
    if len(sys.argv) > 1:
        hx = lambda b: b.hex()
        with open(sys.argv[1], "w") as f:
            f.write("// @generated by tools/ref.py (independent spec-style reference; NOT NIST ACVP)\n")
            f.write("// seed d=00..1f, z=20..3f, m=40..5f; reject case flips ct[0] bit 0\n")
            f.write("pub struct V { pub ek_head: &'static str, pub rho: &'static str, pub dk_hash: &'static str, pub ct_head: &'static str, pub ct_tail: &'static str, pub ss: &'static str, pub ss_reject: &'static str }\n")
            for lvl, d in v.items():
                f.write(f"pub const V{lvl}: V = V {{ ek_head: \"{hx(d['ek_head'])}\", rho: \"{hx(d['rho'])}\", dk_hash: \"{hx(d['dk_hash'])}\", ct_head: \"{hx(d['ct_head'])}\", ct_tail: \"{hx(d['ct_tail'])}\", ss: \"{hx(d['ss'])}\", ss_reject: \"{hx(d['ss_reject'])}\" }};\n")
