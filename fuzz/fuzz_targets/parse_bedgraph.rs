#![no_main]

use bigwig::bedgraph::{Reader, Record};
use libfuzzer_sys::fuzz_target;

// Reads to EOF, alternating `read_record` and the iterator; errors are expected, panics are not.
fuzz_target!(|data: &[u8]| {
    let mut reader = Reader::new(data);
    let mut record = Record::default();
    loop {
        match reader.read_record(&mut record) {
            Ok(false) => break,
            Ok(true) | Err(_) => {}
        }
        if reader.next().is_none() {
            break;
        }
    }
    let _ = reader.track_lines();
});
