# GVC 1.0 release-readiness checklist

This checklist defines the conditions for creating the `release/1.0`
maintenance branch and preparing `1.0.1rc1`.

## Automated gates

The release candidate must satisfy all of the following on the candidate commit:

- [ ] Python 3.8 reviewed dependency stack passes.
- [ ] Python 3.8 unconstrained resolver probe passes.
- [ ] Python 3.9-3.14 compatibility matrix passes.
- [ ] Wheel and sdist build successfully.
- [ ] Wheel and sdist install into isolated environments and import native
      extensions from outside the source checkout.
- [ ] CLI smoke test passes from the installed artifact.
- [ ] Real JBIG-KIT encode/decode gate passes.
- [ ] Pinned LUH historical VCF gate passes with verified upstream blob ID.
- [ ] Multiprocessing fork/spawn correctness and failure tests pass.
- [ ] Random-access tests pass.
- [ ] Byte-exact v1 golden-format tests pass.
- [ ] Benchmark smoke report generation passes.

## Controlled/offline gates

Before tagging `v1.0.1`, record evidence for:

- [ ] Full historical LUH fixture processed, not only the bounded CI block.
- [ ] Production-scale real-JBIG multiprocessing run with realistic worker
      counts and block sizes.
- [ ] Controlled-machine benchmark JSON captured for sequential and parallel
      configurations.
- [ ] Candidate benchmark compared to the retained baseline with an explicitly
      chosen regression budget.
- [ ] Peak memory inspected on the controlled machine.
- [ ] Abrupt termination / scheduler cancellation leaves no final partial
      artifact and no child processes.
- [ ] Random access exercised on a production-scale metadata sidecar.
- [ ] Non-Linux validation performed if Windows/macOS are release targets.
- [ ] Any historical .gvc artifact available to the project has been decoded
      successfully by the candidate.

## Version and branch policy

The `release/1.0` line must retain:

- Python `>=3.8`;
- the v1 serialized format;
- backward-compatible bug/security/compatibility fixes only;
- the reviewed Python 3.8 dependency ceiling;
- golden-format fixtures and historical compatibility gates.

Patch releases on `release/1.0` must not raise the Python floor.

## RC transition

Only after the automated gates and required offline gates are accepted:

1. create/freeze `release/1.0`;
2. change `1.0.1.dev0` to `1.0.1rc1`;
3. build wheel + sdist from a clean checkout;
4. install and test both artifacts;
5. capture final release benchmark JSON;
6. update changelog/release notes;
7. tag `v1.0.1rc1`;
8. after RC acceptance, change to `1.0.1` and tag `v1.0.1`.

Do not start 1.1-only modernization on the release branch.


## Controlled release-validation runner

The expensive release gates are intentionally separate from ordinary PR CI.
Run them locally with:

```bash
python scripts/run_release_validation.py \
  tmp/historical/test_block01.vcf.gz \
  --historical-max-blocks 0 \
  --include-sorting \
  --workers 0 1 2 4 \
  --repetitions 3 \
  --benchmark-output tmp/release-validation/benchmark.json \
  --evidence-output tmp/release-validation/evidence.json
```

If a retained controlled-machine baseline exists, add:

```bash
  --baseline /path/to/baseline.json \
  --max-regression-percent <accepted-budget>
```

The same orchestration is available as the manually dispatched
`Release validation` GitHub Actions workflow. It is never triggered by push
or pull request, so the production-scale historical/benchmark workload does
not consume CI automatically. The workflow uploads the benchmark and evidence
JSON as a temporary Actions artifact.

The evidence JSON records which controlled gates actually ran. It does not
replace human review of the benchmark environment, peak-memory results,
non-Linux validation, or any historical `.gvc` artifact supplied externally.
