use std::{io::{BufReader, BufWriter}, os::unix::net::UnixStream};

use ocelot::svole::{wykw::{self, LPN_EXTEND_SMALL}, SVoleSender};
use scuttlebutt::{AbstractChannel, Channel, field::F61p};

fn main() {
    println!("Hello, world!");
    let (sender, receiver) = UnixStream::pair().unwrap();
    // let sender_reader = BufReader::new();
    let mut sender_channel = Channel::new(
        BufReader::new(&sender),
        BufWriter::new(&receiver),
    );

    let rng = &mut rand::thread_rng();
    let sender= wykw::Sender::<F61p>::init(&mut sender_channel, rng, LPN_EXTEND_SMALL, LPN_EXTEND_SMALL).unwrap();
    sender.send(sender_channel, rng, out);
    // let mut receiver_channel = Channel::new(
    //     BufReader::new(&receiver),
    //     BufWriter::new(&sender),
    // );
    
    // sender_channel.write_bytes(&[0, 1, 2, 3, 4, 5, 6, 7, 8, 9]).unwrap();

    // let mut received: [u8; 10] = [0; 10];
    // receiver_channel.read_bytes(&mut received).unwrap();
    // println!("received {:?}", received);
    // let sender = wykw::Sender::init(channel, rng, lpn_setup, lpn_extend);
}
