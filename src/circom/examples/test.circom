pragma circom 2.1.5;
template TestCircuit () {  

   // Declaration of signals.  
   signal input a[4];  
   signal input b[4];  
   signal output c[3];  
    
   // Constraints.  
   c[0] <== (a[0] + a[3]) * b[0];
   c[1] <== (a[1] + a[2]) * b[1];
   c[2] <== a[3] * b[2] + b[3];
}

component main = TestCircuit();