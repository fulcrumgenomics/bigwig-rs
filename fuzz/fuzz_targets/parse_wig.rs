#![no_main]

use bigwig::wig::{Event, Reader};
use libfuzzer_sys::fuzz_target;

// Reads to EOF, alternating `read_event` and the iterator; errors are expected, panics are not.
fuzz_target!(|data: &[u8]| {
    let mut reader = Reader::new(data);
    let mut event = Event::FixedValue(0.0);
    loop {
        match reader.read_event(&mut event) {
            Ok(false) => break,
            Ok(true) | Err(_) => {}
        }
        if reader.next().is_none() {
            break;
        }
    }
    let _ = reader.track_lines();
});
