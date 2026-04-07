# Changes

## Allow two loss events per region per lineage (required for CN=0)

COMPASS has a hard prior penalty (−1,000,000) for trees where any region is affected by more
than one CNA event along a root-to-leaf lineage path (`rec_check_max_one_event_per_region_per_lineage`
in Tree.cpp). This prevented CN=0, which requires two consecutive losses on the same region
(CN 2→1→0) in ancestor-descendant nodes.

### Tree.cpp

- **Line 512** (`rec_check_max_one_event_per_region_per_lineage`, function begins at line 504): Changed the per-lineage
  event limit from 1 to 2. This specifically enables the CN 2→1→0 double-loss scenario.
  This also enables CN=4 (two gains, CN 2→3→4; see below). All other two-event combinations
  (gain+loss, etc.) are net-neutral and rejected by the MCMC.

---


## Allow copy numbers 0 and 4 (homozygous deletion and high-level amplification)

COMPASS previously constrained copy numbers to {1, 2, 3}. The following changes allow regions to reach copy number 0 or 4.

### Node.cpp

- **Line 171**: Relaxed the CNA validity constraint: lower bound from `>= 1` to `>= 0` (enables CN=0), upper bound from `<= 3` to `<= 4` (enables CN=4). Allele-level checks below it still prevent losing an allele that doesn't exist. CN=0 naturally blocks further CNAs on a deleted region (both alleles are 0). CN=4 is a terminal high-amplification state; further gains are blocked by the validity check.

### Scores.cpp

- **Lines 83–85** (`compute_SNV_loglikelihoods`): Added handling for `c_ref=0, c_alt=0`, which occurs when a locus is in a CN=0 region. Previously, the code fell into `if (c_ref==0) c_alt=1`, giving an alt frequency of ~1 (incorrect). Now `c_ref` is set to 1 instead, so the locus is treated as homozygous reference (alt frequency ≈ sequencing error rate), reflecting that all observed reads are background noise. Also fixed a pre-existing typo: `c_ref==1` (no-op comparison) corrected to `c_ref=1` in the `else if (c_alt==0)` branch.

- **Line 197** (`compute_CNA_loglikelihoods`): Added a `+1e-6` floor to `expected_read_count_region` to avoid `log(0)` when a region has CN=0 (which gives `region_proportion=0`). This allows a small number of noise reads to fall on a deleted region without producing NaN or `-inf` in the likelihood.

- **Lines 208, 231** (`get_dropoutref_counts_genotype`, `get_dropoutalt_counts_genotype`): Applied the same `c_ref=0, c_alt=0` fix as in `compute_SNV_loglikelihoods`.

---

## Synthetic tests for CN=0 (`test/run_tests.py`)

Added `test/run_tests.py`, a self-contained test script that generates hand-crafted synthetic data and verifies CN=0 inference end-to-end. The script runs two tests:

- **Test 1 — smoke**: COMPASS exits zero, all expected output files are present, and none contain NaN or inf.
- **Test 2 — CN=0 recovery**: the best tree contains a node with two `Loss Region0` events (CN 2→1→0), the copynumbers TSV confirms CN=0, and ≥70% of the deleted cells are assigned to that node.

Design notes (CN=0): 110 cells (60 normal, 50 deleted); `--filterregions 0` keeps the zero-read region; `--cnacost 20 --lohcost 20` lowers the CNA discovery barrier.

- **Test 3 — CN=4 recovery**: the best tree contains a node with two `Gain Region0` events (CN 2→3→4), the copynumbers TSV confirms CN=4, and ≥70% of the amplified cells are assigned to that node.

Design notes (CN=4): 360 cells (60 normal, 300 amplified). Both gains are on *different alleles* — (Region0, +1, allele=0) and (Region0, +1, allele=1) — so they produce distinct CNA event tuples and coexist in a single node. SNV0 is a het mutation in amplified cells (50:50 reads); at CN=4 with one alt allele per allele pair the expected alt frequency is exactly 0.5, which matches the data and rules out the CNLOH-at-root alternative. Four chains × 5000 steps with `--cnacost 20 --lohcost 20`.
