#!/usr/bin/env python3
"""Extract antigen sequences for every AACDB complex using the HISTORICAL rule.

Historical rule (scan_candidates.py):
    n_ab    = len(parse_chains(target))            # 'XLY' -> [X,L]; else [X]
    ab_vals = sorted(set(chain_encoding))[:n_ab]
    ag      = residues whose chain-encoding is NOT an antibody chain
    ag_seq  = ag[:1200]                            # TRUNCATED at 1200
    header  = the complex id, e.g. 1IC7_HLY

Features (.pt) exist for only 4,238 of 7,695 complexes, so the antigen is read from the
PDB with parse_antigen_residues() -- the same parser the features were built from. The
--validate mode checks PDB-derived == feature-derived on the targets where both exist.
"""
import os, sys, csv
sys.path.insert(0, f"{VAAD_ROOT}/projects/mutation_sampling/lyra/scripts/nos_diffusion")
sys.path.insert(0, f"{VAAD_ROOT}/projects/mutation_sampling/lyra/src")
from cluster_steer import parse_chains, parse_antigen_residues
import os as _os
_VAAD_DEFAULT = _os.path.expanduser("~/vaad-data")
VAAD_ROOT = _os.environ.get("VAAD_ROOT", _VAAD_DEFAULT)

PDBS = f"{VAAD_ROOT}/datasets/AACDB/data_zip/complex_structure"
TRUNC = 1200

def ag_from_pdb(target):
    chains_str = target.split("_")[-1]
    antibody = parse_chains(target)
    antigen_chains = [c for c in chains_str if c not in antibody]
    if not antigen_chains:
        return None, "no_antigen_chain"
    ordered = parse_antigen_residues(f"{PDBS}/{target}.pdb", antigen_chains)
    seq = "".join(aa for _, _, _, aa in ordered)
    # keep "X" gap placeholders: the feature block retains them, so stripping would desync
    if len(seq) < 30:
        return None, f"tiny_antigen:{len(seq)}"
    return seq[:TRUNC], None

if __name__ == "__main__":
    mode = sys.argv[1]
    if mode == "validate":
        # compare against the historical ag_seq recorded in candidates.tsv
        n_ok = n_bad = 0; bad = []
        for r in csv.DictReader(open(f"{VAAD_ROOT}/scratch/genzoo/cand/candidates.tsv"), delimiter="\t"):
            got, err = ag_from_pdb(r["target"])
            if got is None:
                n_bad += 1; bad.append((r["target"], err)); continue
            if got == r["ag_seq"]: n_ok += 1
            else:
                n_bad += 1
                bad.append((r["target"], f"len {len(got)} vs {len(r['ag_seq'])}"))
        print(f"PDB-extraction == historical feature-extraction: {n_ok}/{n_ok+n_bad}")
        for t, why in bad[:10]: print("   MISMATCH", t, why)
    else:
        targets = sorted(f[:-4] for f in os.listdir(PDBS) if f.endswith(".pdb"))
        out, skipped = sys.argv[2], {}
        n = 0
        with open(out, "w") as fh:
            for t in targets:
                try: seq, err = ag_from_pdb(t)
                except Exception as e: seq, err = None, f"err:{type(e).__name__}"
                if seq is None:
                    skipped[err.split(":")[0]] = skipped.get(err.split(":")[0], 0) + 1
                    continue
                fh.write(f">{t}\n{seq}\n"); n += 1
        print(f"wrote {n}/{len(targets)} antigen sequences -> {out}")
        print("skipped:", skipped)
