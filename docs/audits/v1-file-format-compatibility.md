# GVC v1 file-format compatibility audit

This audit defines the compatibility contract for the GVC 1.0.x maintenance
line.

## Required invariants

GVC 1.0.x must preserve:

- byte-exact serialization for the maintained v1 golden structures;
- decoding of the retained v1 structural fixture corpus;
- parameter-set and access-unit framing;
- per-plane sorting/transpose/codec flags;
- row-bin-split optional payload ordering;
- phase payload ordering and phase convention;
- mixed ploidy and haploid reconstruction;
- random-access behavior for maintained metadata sidecars.

## Malformed-stream rejection

The maintained automated suite requires deterministic rejection of:

- every byte truncation of the retained structural fixtures;
- declared data-unit length mismatches;
- access-unit lengths smaller than the fixed header;
- unknown data-unit/content identifiers where applicable;
- access units referencing an unavailable parameter set;
- trailing bytes after a complete v1 structural stream;
- invalid, duplicate, truncated, and trailing permutation payload data;
- inconsistent native decoder dimensions and invalid ploidy inputs.

Malformed input may raise a format-specific ValueError/EOFError (or an
equivalent explicitly tested parse failure), but it must not be silently
accepted as a valid v1 stream.

## Native/reference parity

Native Cython and standalone libgvc helpers are regression-checked against
Python reference behavior for:

- query-column expansion;
- row-bin-split decoding;
- permutation decoding;
- uniform phase reconstruction.

The native path must reject malformed inputs at least as strictly as the
maintained Python reference behavior.

## Historical compatibility

The pinned LUH historical VCF gate verifies that current 1.0.x encoding and
decoding behavior remains compatible with the reviewed upstream fixture.

Any historical .gvc artifact available outside the repository remains an
offline release gate: decode it with the candidate build and retain the
evidence with the release record.

## Release rule

A 1.0.x patch release must not intentionally change the v1 serialized format.
Any proposed change that alters golden bytes, framing, identifier semantics, or
payload ordering requires a new compatibility decision rather than a patch-only
release.
