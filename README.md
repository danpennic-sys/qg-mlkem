# qg-mlkem
Sealed no_std ML-KEM (FIPS 203) 512/768/1024, zero deps. `build.rs` is the compiler: it proves the NTT/Barrett/compress
constants exhaustively and emits one module per parameter set. `cargo test --release`; `python3 tools/ref.py tests/vectors.inc`
regenerates cross-check vectors. **Uncompiled draft** (no cargo in the authoring sandbox); vectors are from my own reference, not NIST ACVP.
