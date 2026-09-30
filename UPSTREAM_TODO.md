# Upstream TODO — measurement defects to resolve before the next production run

Opened 2026-09-29. Items are added here only when explicitly selected after
discussion. Candidates found during the 29 September audit are parked in
`scratchpad/audit_findings.md` and are deliberately NOT listed here.

Legend: `[ ]` open · `[x]` resolved · **RECAL** = fixing it invalidates prior runs.

---

## [ ] 1. Remove the sequencing mechanism; persist diagnosis times  **RECAL**

Agreed before this audit. `sequencing_rate = 0.05` yields ~20 sequences per
simulation out of ~3,500 infections; in SOT and edge_case the median
simulation sequences no long shedder at all. `seq_rate` never enters a
propensity, so the mechanism is separable from the dynamics.

**Do:** drop the in-simulation sequencing draw; persist `t_diagnosis` per
individual so any sequencing ratio can be reconstructed post hoc from the
complete intra-host record.

**Why first:** it removes a sampling layer that every downstream diagnostic
has to reason around, and it is the only item here that changes what data
exists rather than what a number means.

---
