circom --r1cs --wasm test.circom
node test_js/generate_witness.js test_js/test.wasm input.json witness.wtns

