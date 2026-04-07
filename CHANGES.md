# Changes

## Allow two loss events per region per lineage (required for CN=0)

COMPASS has a hard prior penalty (−1,000,000) for trees where any region is affected by more
than one CNA event along a root-to-leaf lineage path (`rec_check_max_one_event_per_region_per_lineage`
in Tree.cpp). This prevented CN=0, which requires two consecutive losses on the same region
(CN 2→1→0) in ancestor-descendant nodes.

### Tree.cpp

- **Line 512** (`rec_check_max_one_event_per_region_per_lineage`, function begins at line 504): Changed the per-lineage
  event limit from 1 to 2. This specifically enables the CN 2→1→0 double-loss scenario.
  All other two-event combinations (gain+gain, loss+gain, etc.) are harmless: additional gains
  beyond CN=3 are already rejected by `update_genotype`'s validity check, and net-neutral
  combinations have no likelihood benefit and will be rejected by the MCMC.

---


## Allow copy number 0 (homozygous deletion) events

COMPASS previously constrained copy numbers to {1, 2, 3}. The following changes allow regions to reach copy number 0.

### Node.cpp

- **Line 171**: Relaxed the CNA validity constraint from `cn_regions[region]+gain_loss >= 1` to `>= 0`. This is the single gate that blocked homozygous deletions. Allele-level checks below it still prevent losing an allele that doesn't exist, and naturally block further CNAs on an already-deleted region (since both n_ref and n_alt are 0).

### Scores.cpp

- **Lines 83–85** (`compute_SNV_loglikelihoods`): Added handling for `c_ref=0, c_alt=0`, which occurs when a locus is in a CN=0 region. Previously, the code fell into `if (c_ref==0) c_alt=1`, giving an alt frequency of ~1 (incorrect). Now `c_ref` is set to 1 instead, so the locus is treated as homozygous reference (alt frequency ≈ sequencing error rate), reflecting that all observed reads are background noise. Also fixed a pre-existing typo: `c_ref==1` (no-op comparison) corrected to `c_ref=1` in the `else if (c_alt==0)` branch.

- **Line 197** (`compute_CNA_loglikelihoods`): Added a `+1e-6` floor to `expected_read_count_region` to avoid `log(0)` when a region has CN=0 (which gives `region_proportion=0`). This allows a small number of noise reads to fall on a deleted region without producing NaN or `-inf` in the likelihood.

- **Lines 208, 231** (`get_dropoutref_counts_genotype`, `get_dropoutalt_counts_genotype`): Applied the same `c_ref=0, c_alt=0` fix as in `compute_SNV_loglikelihoods`.

---

## Synthetic tests for CN=0 (`test/run_tests.py`)

Added `test/run_tests.py`, a self-contained test script that generates hand-crafted synthetic data and verifies CN=0 inference end-to-end. The script runs two tests:

- **Test 1 — smoke**: COMPASS exits zero, all expected output files are present, and none contain NaN or inf.
- **Test 2 — CN=0 recovery**: the best tree contains a node with two `Loss Region0` events (CN 2→1→0), the copynumbers TSV confirms CN=0, and ≥70% of the deleted cells are assigned to that node.

Design notes: 110 cells (60 normal, 50 deleted) to satisfy COMPASS's hard thresholds for phase-2 CNA inference; `--filterregions 0` keeps the zero-read region; `--cnacost 20 --lohcost 20` lowers the CNA discovery barrier for this small dataset.
