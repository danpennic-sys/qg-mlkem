# qg-mlkem
Sealed no_std ML-KEM (FIPS 203) 512/768/1024, zero deps. `build.rs` is the compiler: it proves the NTT/Barrett/compress
constants exhaustively and emits one module per parameter set. `cargo test --release`; `python3 tools/ref.py tests/vectors.inc`
regenerates cross-check vectors. Draft: vectors are from `tools/ref.py`, not NIST ACVP.

Licensed under either of Apache-2.0 or MIT, at your option.
