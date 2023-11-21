# VOLE in the head aimed for compatability with popular DSLs
This is an extremely fast NIZK prover capable of doing roughly one million constraints per second (benchmark pending after code is completed) on consumer hardware. It is not succinct. It uses the Quicksilver [] proving system with VOLE in the head [] as a linearly homorphic commitment scheme. It is optimized for the prime 21888242871839275222246405745257275088548364400416034343698204186575808495617 popular modern proving systems. We plan to support more finite fields.

# How to use
Supply R1CS and witness via CLI or wasm

# How witness is formatted
Geared to be compatible with arbitrary witnesses and R1CS. Witness and their corresponding R1CS lengths are padded with zeroes to be divisible by the number of small VOLEs which is currently 1024. An extra 1024 entries are added to the VOLE so the universal hash in the subsapce VOLE correction does not reveal the linear combination of original VOLE inputs.

# 

## Keccak circuits in GF2^k, combined with ECDSA circuits over secp256k1's prime field 

