# GVC 1.0.1 release-candidate checklist

This checklist is intentionally separate from normal CI. Do not create the
release tag until every item below refers to the same commit.

## Promotion

1. Start from a fully green `feat/ci-versioning` / release-preparation commit.
2. Preserve `requires-python = ">=3.8"` and the v1 serialization contract.
3. Change `gvc/_version.py` from `1.0.1.dev0` to `1.0.1rc1`.
4. Move relevant Unreleased changelog entries under a `1.0.1rc1` heading.
5. Run:

   ```bash
   python scripts/check_release_readiness.py --rc
   python scripts/check_release_tag.py v1.0.1rc1
   python scripts/ci.py all
   ```

## Controlled release evidence

Run the manual Release validation workflow with:

- `historical_max_blocks=0` for the full LUH fixture;
- sorting enabled for the final RC validation;
- controlled benchmark worker counts appropriate to the machine;
- `rc_mode=true`;
- a retained benchmark baseline and explicit regression budget when a trusted
  same-machine baseline exists.

Retain:

- `evidence.json`;
- `artifacts.json`;
- wheel and source distribution;
- benchmark JSON;
- historical-validation logs.

The artifact evidence records SHA-256 hashes. The artifacts eventually
published must be the exact files represented by those hashes; do not rebuild
them after approval.

## Tagging

Only after evidence review:

```bash
python scripts/check_release_tag.py v1.0.1rc1
git tag -a v1.0.1rc1 -m "GVC 1.0.1rc1"
```

The tag must point to the exact commit whose CI and release-validation evidence
were reviewed.

## Final 1.0.1

After RC validation and any necessary fixes, repeat the same process with
`1.0.1` and tag `v1.0.1`. Do not move the Python floor or v1 file-format
contract in a 1.0.x patch release.
