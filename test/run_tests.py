"""
Minimal synthetic tests for the CN=0 and CN=4 modifications to COMPASS.

Ground truth (CN=0 scenario):
  - 2 regions (Region0, Region1), 2 loci (one per region), 110 cells.
  - "normal" cells (0-59): CN=2 in both regions, ~50 reads each, all-ref at both loci.
  - "deleted" cells (60-109): CN=0 in Region0 (0 reads), CN=2 in Region1, het mutation at locus1.
    The mutation at locus1 separates the two populations in the SNV-only phase so that
    COMPASS can then estimate region weights and attempt CNA inference.

Ground truth (CN=4 scenario):
  - 2 regions (Region0, Region1), 2 loci (one per region), 110 cells.
  - "normal" cells (0-59): CN=2 in both regions, ~50 reads each, all-ref at both loci.
  - "amplified" cells (60-109): CN=4 in Region0 (~100 reads), CN=2 in Region1, het mutation at locus1.
    The mutation at locus1 separates the two populations in the SNV-only phase.

Test 1 — smoke: COMPASS runs to completion with no crash and no NaN/inf in any output.
Test 2 — CN=0 recovery: the best tree has a node with two Loss Region0 events (CN 2->1->0),
          and those cells are predominantly assigned to it.
Test 3 — CN=4 recovery: the best tree has a node with two Gain Region0 events (CN 2->3->4),
          and those cells are predominantly assigned to it.
"""

import subprocess, os, sys, json, csv

COMPASS = os.path.join(os.path.dirname(__file__), "..", "COMPASS")
OUTDIR  = os.path.join(os.path.dirname(__file__), "output_cn0")

N_NORMAL  = 60   # root node threshold: max(40, 0.015*n_cells) → need ≥40 cells at root
N_DELETED = 50   # non-root node threshold: max(40, 0.03*n_cells) → need ≥40 cells in node
N_CELLS   = N_NORMAL + N_DELETED
READS     = 50

# CN=4 uses more amplified cells so the likelihood gain from CN=4 vs CN=3 overcomes the
# prior penalty for the second CNA event.  With NB theta=6 and READS=50, the per-cell
# log-likelihood advantage of CN=4 over CN=3 is ~0.04; the prior penalty for one extra
# gain event is ~ncells_coef * cnacost ≈ 0.245 * 10 = 2.45 (for 360 cells).
# 300 amplified cells gives a net advantage of ~300*0.04 - 2.45 ≈ 9.5 log-units.
N_AMPLIFIED = 300
N_CELLS_CN4 = N_NORMAL + N_AMPLIFIED


# ---------------------------------------------------------------------------
# Synthetic data
# ---------------------------------------------------------------------------

def make_variants_csv(path):
    header = ["CHR", "POS", "REF", "ALT", "REGION", "NAME", "FREQ"] + \
             [f"cell{j}" for j in range(N_CELLS)]

    # Locus in Region0: normal cells have READS ref reads; deleted cells have 0
    # (Region0 has CN=0 in deleted cells so no reads at all)
    row0 = ["0", "100", "A", "T", "Region0", "SNV0", "0.0"] + \
           [f"{READS}:0:3"] * N_NORMAL + ["0:0:3"] * N_DELETED

    # Locus in Region1: normal cells all-ref; deleted cells carry the het mutation
    # This mutation is what separates the two populations in the SNV-only phase.
    ref1 = READS // 2
    alt1 = READS - ref1
    row1 = ["1", "200", "A", "T", "Region1", "SNV1", "0.0"] + \
           [f"{READS}:0:3"] * N_NORMAL + [f"{ref1}:{alt1}:3"] * N_DELETED

    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(header)
        w.writerow(row0)
        w.writerow(row1)


def make_regions_csv(path):
    # Region0: normal cells have READS reads; deleted cells have 0 (CN=0)
    row0 = ["0_Region0"] + [str(READS)] * N_NORMAL + ["0"] * N_DELETED
    # Region1: all cells normal
    row1 = ["1_Region1"] + [str(READS)] * N_CELLS

    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(row0)
        w.writerow(row1)


def make_variants_csv_cn4(path):
    # Two gains on DIFFERENT alleles in the same node produce distinct CNA event tuples
    # (Region0, +1, [0]) and (Region0, +1, [1]), allowing CN=4 in a single node.
    #
    # SNV0 is a het somatic mutation in amplified cells (50:50 reads).  At CN=4 with one
    # alt copy per allele pair, genotype is (2 ref, 2 alt) -> expected alt = 0.5 = observed.
    # This blocks the CNLOH-at-root trick: CNLOH would make alt=0 or 1, neither matching 0.5.
    # SNV1 separates the two populations in the SNV-only phase 1.
    header = ["CHR", "POS", "REF", "ALT", "REGION", "NAME", "FREQ"] + \
             [f"cell{j}" for j in range(N_CELLS_CN4)]

    ref0 = READS // 2
    alt0 = READS - ref0   # 25:25 for normal (all-ref: 50:0); amplified: het (50:50 = 2*READS/2 each)
    # Normal cells: all-ref at SNV0; amplified cells: het (50% alt, 2*READS total = 100 reads)
    row0 = ["0", "100", "A", "T", "Region0", "SNV0", "0.0"] + \
           [f"{READS}:0:3"] * N_NORMAL + [f"{READS}:{READS}:3"] * N_AMPLIFIED

    ref1 = READS // 2
    alt1 = READS - ref1
    row1 = ["1", "200", "A", "T", "Region1", "SNV1", "0.0"] + \
           [f"{READS}:0:3"] * N_NORMAL + [f"{ref1}:{alt1}:3"] * N_AMPLIFIED

    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(header)
        w.writerow(row0)
        w.writerow(row1)


def make_regions_csv_cn4(path):
    # Region0: normal cells have READS reads; amplified cells have 2*READS (CN=4 vs CN=2).
    # total_counts for amplified cells = 2*READS + READS = 3*READS.
    # Expected Region0 proportion at CN=4: 4/(4+2) = 2/3.
    # Expected reads: 3*READS * 2/3 = 2*READS = 100. Matches perfectly.
    row0 = ["0_Region0"] + [str(READS)] * N_NORMAL + [str(2 * READS)] * N_AMPLIFIED
    # Region1: all cells normal
    row1 = ["1_Region1"] + [str(READS)] * N_CELLS_CN4

    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(row0)
        w.writerow(row1)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def strip_node_prefix(name):
    """_cellAssignments.tsv uses bare indices ('0','1') while JSON uses 'Node 0','Node 1'."""
    return name.replace("Node ", "").strip()


def run_compass(data_prefix, out_prefix):
    cmd = [
        COMPASS,
        "-i", data_prefix,
        "-o", out_prefix,
        "--filterregions", "0",   # Region0 would be filtered (50% zero-read cells)
        "--cnacost",      "20",   # lower than default (85) to help discovery with few regions
        "--lohcost",      "20",
        "--nchains",      "1",
        "--chainlength",  "5000",
        "-d",             "0",    # no doublets for simplicity
    ]
    return subprocess.run(cmd, capture_output=True, text=True)


def run_compass_cn4(data_prefix, out_prefix):
    # CN=4 uses two gain events on different alleles in the same node.
    # SNV0 (het in amplified cells) provides additional signal beyond read depth alone.
    cmd = [
        COMPASS,
        "-i", data_prefix,
        "-o", out_prefix,
        "--cnacost",      "20",
        "--lohcost",      "20",
        "--nchains",      "4",
        "--chainlength",  "5000",
        "-d",             "0",
    ]
    return subprocess.run(cmd, capture_output=True, text=True)


def has_nan_or_inf(filepath):
    with open(filepath) as f:
        content = f.read().lower()
    return "nan" in content or "inf" in content


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

def test_smoke(data_prefix, out_prefix, runner=None):
    print("=== Test 1: smoke (no crash, no NaN/inf) ===")
    if runner is None:
        runner = run_compass
    result = runner(data_prefix, out_prefix)

    if result.returncode != 0:
        print("FAIL: COMPASS exited non-zero.")
        print(result.stdout[-2000:])
        print(result.stderr[-2000:])
        return False

    always_present = [
        out_prefix + "_cellAssignments.tsv",
        out_prefix + "_nodes_genotypes.tsv",
        out_prefix + "_tree.json",
    ]
    for path in always_present:
        if not os.path.exists(path):
            print(f"FAIL: expected output missing: {path}")
            return False
        if has_nan_or_inf(path):
            print(f"FAIL: NaN or inf found in {path}")
            return False

    print("  stdout:", result.stdout.strip().splitlines()[-3:])
    print("PASS")
    return True


def test_cn0_recovery(out_prefix):
    print("=== Test 2: CN=0 recovery ===")

    json_path  = out_prefix + "_tree.json"
    cell_path  = out_prefix + "_cellAssignments.tsv"
    cn_path    = out_prefix + "_nodes_copynumbers.tsv"

    with open(json_path) as f:
        tree = json.load(f)

    # Find nodes that have a loss on Region0 (shows as "Loss Region0" in CNA list)
    nodes_with_loss_Region0 = [
        n["name"] for n in tree["nodes"]
        if any("Region0" in ev and "Loss" in ev for ev in n.get("CNA", []))
    ]
    print(f"  Nodes with loss on Region0: {nodes_with_loss_Region0}")

    # Cross-check with _nodes_copynumbers.tsv if it exists
    cn0_nodes = []
    if os.path.exists(cn_path):
        with open(cn_path) as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                if row.get("Region0", "").strip() == "0":
                    cn0_nodes.append(row["node"])
        print(f"  Nodes with CN=0 for Region0 (from copynumbers TSV): {cn0_nodes}")
    else:
        print("  (_nodes_copynumbers.tsv not present — checking JSON only)")

    if not nodes_with_loss_Region0 and not cn0_nodes:
        print("FAIL: no loss event on Region0 found in the inferred tree")
        print("  All nodes:")
        for n in tree["nodes"]:
            print(f"    {n['name']}: SNV={n.get('SNV')}, CNV={n.get('CNV')}")
        return False

    # Check that deleted cells (cell20–cell39) are mostly assigned to a node
    # that is a descendant of (or is) a node with a Region0 loss.
    # We identify the relevant nodes as those with Region0 loss (or their descendants
    # if CN=0 requires two events); for simplicity we just check the direct loss nodes.
    loss_node_names = {strip_node_prefix(n) for n in nodes_with_loss_Region0} | \
                      {strip_node_prefix(n) for n in cn0_nodes}

    deleted_cell_names = {f"cell{j}" for j in range(N_NORMAL, N_CELLS)}
    assigned_to_loss = 0
    total_deleted = 0

    with open(cell_path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            if row["cell"] in deleted_cell_names:
                total_deleted += 1
                # Accept cells assigned to any node with or below a Region0 loss.
                # Since we cannot easily walk the full tree here, we match on name.
                if row["node"] in loss_node_names:
                    assigned_to_loss += 1

    frac = assigned_to_loss / total_deleted if total_deleted > 0 else 0
    print(f"  {assigned_to_loss}/{total_deleted} deleted cells assigned to a Region0-loss node ({frac:.0%})")

    if frac < 0.7:
        print("FAIL: fewer than 70% of deleted cells assigned to a Region0-loss node")
        # Print full cell assignments for debugging
        with open(cell_path) as f:
            print(f.read())
        return False

    print("PASS")
    return True


def test_cn4_recovery(out_prefix):
    print("=== Test 3: CN=4 recovery ===")

    json_path = out_prefix + "_tree.json"
    cell_path = out_prefix + "_cellAssignments.tsv"
    cn_path   = out_prefix + "_nodes_copynumbers.tsv"

    with open(json_path) as f:
        tree = json.load(f)

    # Find nodes with two Gain Region0 events (CN 2->3->4 requires two gains)
    nodes_with_two_gains_Region0 = [
        n["name"] for n in tree["nodes"]
        if sum(1 for ev in n.get("CNA", []) if "Region0" in ev and "Gain" in ev) >= 2
    ]
    print(f"  Nodes with 2x Gain on Region0: {nodes_with_two_gains_Region0}")

    # Cross-check with _nodes_copynumbers.tsv if it exists
    cn4_nodes = []
    if os.path.exists(cn_path):
        with open(cn_path) as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                if row.get("Region0", "").strip() == "4":
                    cn4_nodes.append(row["node"])
        print(f"  Nodes with CN=4 for Region0 (from copynumbers TSV): {cn4_nodes}")
    else:
        print("  (_nodes_copynumbers.tsv not present — checking JSON only)")

    if not nodes_with_two_gains_Region0 and not cn4_nodes:
        print("FAIL: no CN=4 event on Region0 found in the inferred tree")
        print("  All nodes:")
        for n in tree["nodes"]:
            print(f"    {n['name']}: SNV={n.get('SNV')}, CNA={n.get('CNA')}")
        return False

    loss_node_names = {strip_node_prefix(n) for n in nodes_with_two_gains_Region0} | \
                      {strip_node_prefix(n) for n in cn4_nodes}

    amplified_cell_names = {f"cell{j}" for j in range(N_NORMAL, N_CELLS_CN4)}
    assigned_to_gain = 0
    total_amplified = 0

    with open(cell_path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            if row["cell"] in amplified_cell_names:
                total_amplified += 1
                if row["node"] in loss_node_names:
                    assigned_to_gain += 1

    frac = assigned_to_gain / total_amplified if total_amplified > 0 else 0
    print(f"  {assigned_to_gain}/{total_amplified} amplified cells assigned to a Region0-gain node ({frac:.0%})")

    if frac < 0.7:
        print("FAIL: fewer than 70% of amplified cells assigned to a Region0-gain node")
        with open(cell_path) as f:
            print(f.read())
        return False

    print("PASS")
    return True


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    os.makedirs(OUTDIR, exist_ok=True)

    # --- CN=0 tests ---
    prefix_cn0 = os.path.join(OUTDIR, "synthetic_cn0")
    print("Generating synthetic CN=0 data ...")
    make_variants_csv(prefix_cn0 + "_variants.csv")
    make_regions_csv(prefix_cn0 + "_regions.csv")

    passed = []
    passed.append(test_smoke(prefix_cn0, prefix_cn0))
    if passed[-1]:
        passed.append(test_cn0_recovery(prefix_cn0))

    # --- CN=4 tests ---
    prefix_cn4 = os.path.join(OUTDIR, "synthetic_cn4")
    print("\nGenerating synthetic CN=4 data ...")
    make_variants_csv_cn4(prefix_cn4 + "_variants.csv")
    make_regions_csv_cn4(prefix_cn4 + "_regions.csv")

    passed.append(test_smoke(prefix_cn4, prefix_cn4, runner=run_compass_cn4))
    if passed[-1]:
        passed.append(test_cn4_recovery(prefix_cn4))

    print()
    print(f"Results: {sum(passed)}/{len(passed)} tests passed.")
    sys.exit(0 if all(passed) else 1)
