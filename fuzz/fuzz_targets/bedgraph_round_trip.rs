#![no_main]

use bigwig::bedgraph::{Reader, Record, Writer};
use libfuzzer_sys::fuzz_target;

/// Equality that treats any two NaNs as equal, since `Record`'s `PartialEq` does not.
fn same_record(a: &Record, b: &Record) -> bool {
    a.chrom == b.chrom
        && a.start == b.start
        && a.end == b.end
        && (a.value.to_bits() == b.value.to_bits() || (a.value.is_nan() && b.value.is_nan()))
}

// Every record that parses must write and re-parse to an equal record.
fuzz_target!(|data: &[u8]| {
    for parsed in Reader::new(data).flatten() {
        let mut writer = Writer::new(Vec::new());
        writer
            .write_record(&parsed.chrom, parsed.start, parsed.end, parsed.value)
            .expect("a parsed record must be writable");
        let text = writer.into_inner().expect("writing to a Vec cannot fail");

        let reparsed: Vec<Record> = Reader::new(text.as_slice())
            .collect::<Result<_, _>>()
            .expect("written output must parse");
        assert_eq!(reparsed.len(), 1, "one record in, one record out");
        assert!(same_record(&parsed, &reparsed[0]), "{parsed:?} became {:?}", reparsed[0]);
    }
});
