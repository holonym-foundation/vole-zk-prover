# VOLE in the head aimed for compatability with popular DSLs
This is an extremely efficient NIZK prover. It is currently doing about 300k constraints per second on consumer hardware on a 254-bit field. There is a tradeoff in that is not succinct. It uses the [Quicksilver](https://eprint.iacr.org/2021/076) proving system with [VOLE-in-the-head](https://eprint.iacr.org/2023/996) for the commitment scheme. It is optimized for the prime 21888242871839275222246405745257275088548364400416034343698204186575808495617 popular modern proving systems. We plan to support more finite fields.

# How to use
To obtain effeciency benefits of VitH with Quicksilver for a circom circuit, it's quite simple: simply pass the R1CS and witness as arguments to the prover and verifier. No verification key or proving key is necessary. For a rust example, take a look at the prover and verifier in `actors.rs`. Command line and WASM examples and interfaces do not exist, but pull requests with those are quite welcome. 


# How this is organized
- Linear code security (i.e. minimum distance) calculations are in `codeparams/`. They are slow for large codes.
- Quicksiver proving system is in `zkp/`
- Vector commitments used for creating the initial VOLE ("small VOLE" with delta chosen from a tiny set) are in `veccom/`. This is a relatively simple hash-based vector commitment.
- The small VOLE with delta chosen from a small set is in `smallvole/`. Here, the set delta can be chosen from has the smallest cardinality possible, i.e. 2
- The small VOLEs are then "stacked" then transformed into a subspace VOLE in `subspacevole/`. In subspace VOLE, the small VOLE values become interconnected in that they are parts of codewords. If a prover cheats in one small VOLE by guessing delta from the tiny set, it won't be part of the codeword anymore. So it ensures he can't cheat without guessing $d$ deltas where $d$ is the minimum distance of the code.
- Fiat-Shamir Heuristic is generated in `challenges/` to render the proof noninteractive

# Known Issues
- When interpereting circom circuits, wire to labels map is currently assumed to be the identity map which could cause some circuits with different maps to fail
- Linear code (RMA code) security parameter is currently too low. The computation of how large the block size should be for a given security amount is expensive -- once this is computed, the security parameter will be increased. We do not expect performance to change in any practical way from changing this parameter.

