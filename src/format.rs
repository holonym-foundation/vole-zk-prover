//! Reads and write proof formats.
use crate::{actors::actors::Proof, Fr, FrRepr, PF};
use ff::PrimeField;
use serde::{ser::{Serialize, Serializer, SerializeStruct}, de::{Deserialize, Visitor, Expected}};

// extern crate proc_macro;
// use proc_macro::TokenStream;


// pub trait SerializeEfficiently {
//     fn serialize_efficiently() -> Vec<u8>;
// }
// pub trait DeserializeEfficiently {
//     fn deserialize_efficiently(b: &Vec<u8>) -> Self;
// }

// #[proc_macro_derive(SerializeEfficiently)]
// pub fn derive_serialize_efficiently(_item: TokenStream) -> TokenStream {
//     "fn serialize_efficiently() -> Vec<u8> { vec![69; 69] }".parse().unwrap()
// }


impl<'a> Serialize for Fr {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
        where
            S: Serializer {
        serializer.serialize_bytes(&self.to_repr().0)
    }
}
impl<'de> Deserialize<'de> for Fr {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
        where
            D: serde::Deserializer<'de> {
        deserializer.deserialize_bytes(FrVisitor)
    }
}

struct FrVisitor;
impl<'de> Visitor<'de> for FrVisitor {
    type Value = Fr;
    fn expecting(&self, formatter: &mut std::fmt::Formatter) -> std::fmt::Result {
        formatter.write_str("Fr representation bytes (For BN254, 32 little-endian bytes)")
    }
    fn visit_bytes<E>(self, v: &[u8]) -> Result<Self::Value, E>
        where
            E: serde::de::Error, {
        let repr = FrRepr(v.try_into().map_err(|_|E::invalid_length(v.len(), &"32"))?);
        let f = Fr::from_repr(repr);
        if f.is_none().into() { 
            Err(E::invalid_value(serde::de::Unexpected::Bytes(v), &"valid representation of a field element"))
        } else { 
            Ok(f.unwrap())
        }
    }
}

#[cfg(test)]
mod test {
    use ff::Field;
    use rand::thread_rng;

    use crate::FVec;

    use super::*;
    #[test]
    fn serde_fr() {
        let zero = Fr::ZERO;
        let one = Fr::ONE;
        let negone = zero - one;
        let rand = Fr::random(&mut thread_rng());
        let test_cases = vec![zero, one, negone, rand]; 
        test_cases.iter().for_each(|x|{
            let s = bincode::serialize(&x).unwrap();
            println!("Binary serialized {:?}", &s);
            let d: Fr = bincode::deserialize(&s).unwrap();
            assert_eq!(*x, d);
        })
    }
    #[test]
    fn serde_fr_vec() {
        let zero = Fr::ZERO;
        let one = Fr::ONE;
        let negone = zero - one;
        let rand = Fr::random(&mut thread_rng());
        let test_cases = vec![zero, one, negone, rand]; 
        let v = FVec::<Fr>(test_cases);
        let s = bincode::serialize(&v).unwrap();
        let d: FVec<Fr> = bincode::deserialize(&s).unwrap();
        assert_eq!(v, d);
    }
}