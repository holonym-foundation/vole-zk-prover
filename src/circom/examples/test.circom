pragma circom 2.1.5;
include "node_modules/circomlib/circuits/sha256/sha256.circom";
include "node_modules/circomlib/circuits/bitify.circom";

template TestCircuit(Sha256IterationsMinusOne) {  
   signal input testInputs[3];
   signal input publicInput;
   signal input anotherPublicInput;
   signal output digest[256];
   signal output unimportantOutput;

   signal sig <== publicInput * testInputs[0] + testInputs[1] + testInputs[2];
   unimportantOutput <== sig * anotherPublicInput;

   component n2b = Num2Bits(254);

   n2b.in <== testInputs[0] + testInputs[0] + testInputs[0] + testInputs[0];
   
   component initialHash = Sha256(256);
   initialHash.in[0] <== 0;
   initialHash.in[1] <== 0;
   var i;
   for(i=2; i<256; i++) {
      initialHash.in[i] <== n2b.out[i-2];
   }

   component newHashes[Sha256IterationsMinusOne]; 
   newHashes[0] = Sha256(256);
   newHashes[0].in <== initialHash.out;
   var j;
   for(j=1; j<Sha256IterationsMinusOne; j++) {
      newHashes[j] = Sha256(256);
      newHashes[j].in <== newHashes[j-1].out;
   }

   digest <== newHashes[Sha256IterationsMinusOne-1].out;
   
   // // Declaration of signals.  
   // signal input a[4];  
   // signal input b[4];  
   // signal output c[3];  
    
   // // Constraints.  
   // c[0] <== (a[0] + a[3]) * b[0];
   // c[1] <== (a[1] + a[2]) * b[1];
   // c[2] <== a[3] * b[2] + b[3];
}

component main { public [publicInput, anotherPublicInput] } = TestCircuit(9);
