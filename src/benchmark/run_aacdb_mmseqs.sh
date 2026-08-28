#!/bin/bash
: "${VAAD_ROOT:?set VAAD_ROOT to your data root}"
# Usage: bash run_aacdb_mmseqs.sh AACDB_all_antigens.fasta benchmark_200_antigens.fasta
# Historical clustering rule, applied identically to both files.
set -e
MMSEQS=${VAAD_ROOT}/tools/mmseqs/bin/mmseqs
ALL=$1; BENCH=$2
PARAMS="--min-seq-id 0.7 -c 0.8 --cov-mode 0 --cluster-mode 0"

rm -rf mmseqs_aacdb_full mmseqs_bench_only
mkdir -p mmseqs_aacdb_full mmseqs_bench_only

echo "== mmseqs version =="; $MMSEQS version

echo "== clustering FULL AACDB =="
$MMSEQS easy-cluster "$ALL" mmseqs_aacdb_full/aacdb70 mmseqs_aacdb_full/tmp $PARAMS -v 1

echo "== clustering BENCHMARK ONLY (sanity check: expect 172) =="
$MMSEQS easy-cluster "$BENCH" mmseqs_bench_only/bench70 mmseqs_bench_only/tmp $PARAMS -v 1

python3 - "$ALL" "$BENCH" <<'PY'
import sys
all_fa, bench_fa = sys.argv[1], sys.argv[2]
def members(p):
    m={}
    for line in open(p):
        rep,mem=line.rstrip("\n").split("\t"); m[mem]=rep
    return m
full  = members("mmseqs_aacdb_full/aacdb70_cluster.tsv")
bench = members("mmseqs_bench_only/bench70_cluster.tsv")
panel = [l[1:].strip() for l in open(bench_fa) if l.startswith(">")]
n_all = sum(1 for l in open(all_fa) if l.startswith(">"))
N = len(set(full.values()))
B = len(set(bench.values()))
covered = {full[t] for t in panel if t in full}
K = len(covered)
print()
print(f"Full AACDB sequences clustered:              {n_all}")
print(f"Full AACDB clusters:                         {N}")
print(f"Benchmark clusters when clustered alone:     {B}")
print(f"Full-AACDB clusters represented by our 200:  {K}")
print(f"Cluster coverage: {K}/{N} = {100*K/N:.2f}%")
print(f"Complex coverage: {len(panel)}/{n_all} = {100*len(panel)/n_all:.2f}%")
# how much merging happened when the 200 are embedded in the full set
print(f"\nMerging check: our 200 occupy {B} clusters alone but only {K} within full AACDB "
      f"({B-K} benchmark-only clusters merged via additional AACDB sequences)")
PY
