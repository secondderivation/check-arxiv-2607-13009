# Model-Independent Indication for a Localized Anomaly in the Late-Time Expansion History

[2607.13009v1](https://arxiv.org/abs/2607.13009v1) · Review 0.2.0-rc2

**Partly reproduced.** Exact pressure, distance and NEC identities reproduce; saved numerical agreement is conditional on our implementation.

**Limit:** Headline significance, location and full author-pipeline fidelity remain unverified; the numerical reconstruction underwent informed repair.

## Checked scope

- Exact pressure, distance and NEC identities reproduce; saved numerical agreement is conditional on our implementation.

Not checked: headline significance and location; H0 and h\*r\_d numerical fidelity; full robustness fits; printed mock validation table; global calibration of authors statistic

## Review record

| Claim | Method | Outcome / adequacy | Evidence |
|---|---|---|---|
| A1: Exact pressure/distance/NEC identities — Exact identities reproduce under stated flat-FLRW/GR premises. | independent\_derivation | held:reproduced / pass | [reference/algebra.json](reference/algebra.json) |
| N1: Specified implementation and pointwise NEC summaries — Numerical agreement is conditional on this implementation and informed repair, not author-pipeline fidelity. | rerun | unresolved / inconclusive | [reference/I2_RA_20261009.json](reference/I2_RA_20261009.json), [reference/I2_RA_731009.json](reference/I2_RA_731009.json), [reference/I2_REIRA_20261009.json](reference/I2_REIRA_20261009.json), [reference/I2_REIRA_731009.json](reference/I2_REIRA_731009.json), [reference/I4_RA_20261009.json](reference/I4_RA_20261009.json), [reference/I4_RA_731009.json](reference/I4_RA_731009.json), [reference/I4_REIRA_20261009.json](reference/I4_REIRA_20261009.json), [reference/I4_REIRA_731009.json](reference/I4_REIRA_731009.json) |

[Protocol and provenance](PROTOCOL.md). This review record is retrospective.

## Review provenance

Prepared with AI assistance. Historical Claude first pass and Codex independent review/critique remain distinct; this assembly reads saved evidence and is neither blind B nor new P2. No scientific code was run in this assembly. Before the replay command below, extract legacy-original.zip into a new directory and read its README.md and requirements.txt for installation, dependency, data and scope instructions.

Tools: Claude (historical A); Codex (historical B/C and current informed assembly). Human involvement: Owner selected preparation scope; no publication decision is recorded for this batch..

## Replay

~~~sh
mkdir replay-work && python3 -m zipfile -e legacy-original.zip replay-work && cd replay-work && python3 reproduce.py --all
~~~

See original source-payload README and recorded P2; not rerun during assembly.. Unpacks the unchanged original candidate into a new work directory and runs its documented scoped replay. Package validation checked the artifact and inherited the earlier replay after verifying the relevant files were unchanged; it did not rerun the scientific calculation.
