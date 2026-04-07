# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## About

This is a fork of [COMPASS](https://github.com/cbg-ethz/COMPASS), a C++ tool for joint copy number alteration (CNA) and SNV phylogeny reconstruction from MissionBio Tapestri amplicon single-cell sequencing data. See `CHANGES.md` for a log of modifications from the original and `README_original.md` for full usage documentation.

## Build

Requires GCC with OpenMP. On the cluster:

```bash
ml stack/.2024-05-silent gcc/13.2.0
make
```

Clean object files: `make clean`

## Running

```bash
./COMPASS -i [sample_name] -o [output_name] --nchains 4 --chainlength 5000 --CNA 1
```

- `-i`: prefix for input files (`[sample_name]_variants.csv` and `[sample_name]_regions.csv`)
- `-o`: output prefix (produces `_tree.gv`, `_tree.json`, `_cellAssignments.tsv`, etc.)
- `--CNA 1`: enable CNA inference (requires `_regions.csv`)
- `--sex female|male`: chromosome X ploidy (default: female)

Key tuning parameters: `--cnacost` (default 85), `--lohcost` (default 85), `--nodecost` (default 1), `--nchains`, `--chainlength`, `--burnin`.

## Architecture

The codebase is organized around four main classes, with global state declared in `COMPASS.cpp` and shared via `extern`:

- **`Structures.h`** — Plain data structs: `Cell` (per-cell read counts), `Data` (locus/region metadata), `Params` (model hyperparameters). These are globals (`n_cells`, `n_loci`, `n_regions`, `cells`, `data`, `parameters`) accessed via `extern` throughout.

- **`Scores`** (`Scores.h/cpp`)  — Stateless likelihood calculator with aggressive caching (`std::map`-based). Computes beta-binomial SNV log-likelihoods (`compute_SNV_loglikelihoods`) and negative-binomial CNA log-likelihoods (`compute_CNA_loglikelihoods`). A single `Scores*` instance is shared across the tree via pointer.

- **`Node`** (`Node.h/cpp`)  — Represents one clone in the phylogeny. Stores somatic mutations (locus indices), CNA events (region, type ±1/0, affected alleles), and the resulting per-locus allele counts (`n_ref_allele`, `n_alt_allele`) and per-region copy numbers (`cn_regions`). Genotype is propagated from parent via `update_genotype()`. CNA validity is enforced here (the `>= 0` floor for CN=0 support is at Node.cpp:171).

- **`Tree`** (`Tree.h/cpp`)  — Owns a vector of `Node*` plus a parent/children index structure. Manages cell attachment likelihoods (`cells_attach_loglik`), runs EM to update node probabilities and dropout rates, computes the full score (`log_prior_score + log_likelihood`), and implements all MCMC proposal moves (prune-reattach, add/delete node, move/add/remove CNA, etc.). The lineage CNA event limit check is `rec_check_max_one_event_per_region_per_lineage` (Tree.cpp:~504).

- **`Inference`** (`Inference.h/cpp`)  — Runs parallel MCMC chains (via OpenMP) with simulated annealing. Each chain holds a current tree, proposed tree, and best-seen tree. Calls `Tree::mcmc()` which iterates: propose move → score → Metropolis accept/reject.

- **`input`** (`input.h/cpp`)  — Parses `_variants.csv` and `_regions.csv` into the global `cells` and `data` structures.

### Scoring pipeline (per MCMC step)

1. A proposal move modifies `t_prime` (copy of current tree `t`).
2. `Node::update_genotype()` propagates allele counts down the affected subtree.
3. `Node::compute_attachment_scores_parent()` uses `Scores` to compute per-cell log-likelihoods incrementally (only changed nodes).
4. `Tree::compute_likelihood()` runs EM: E-step marginalizes cell attachments, M-step updates node probs and dropout rates.
5. `Tree::compute_prior_score()` applies penalties for node count, CNA events, LOH, and the lineage event limit.
6. Metropolis criterion accepts/rejects based on `log_score` difference and current temperature.

## Fork-specific changes

All modifications are documented in `CHANGES.md`. In brief:
- CN=0 (homozygous deletion) is now supported via two consecutive loss events.
- The per-lineage CNA event limit was raised from 1 to 2 (Tree.cpp:504–508).
- `Scores.cpp` handles `c_ref=0, c_alt=0` edge cases and adds a `+1e-6` floor to avoid `log(0)`.
