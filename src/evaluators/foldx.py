#!/usr/bin/env python3
"""FoldX validation: mutation-DDG-grade binding energies for the cluster-steered designs.

Replaces the PRODIGY pipeline (which was insensitive to point mutations on a fixed
backbone). For each antibody variant (native WT + each design) x antigen variant
(WT + mutants):
  1. renumber the complex once  -> unique sequential resnums, NO insertion codes
                                   (so FoldX mutation strings are unambiguous), keep chains
  2. FoldX RepairPDB once        -> energy-minimized reference
  3. FoldX BuildModel per pair   -> apply CDR design + antigen point mutations (FoldX models
                                   the sidechains + repacks neighbours; sensitive to identity)
  4. FoldX AnalyseComplex        -> interaction energy between antibody (H,L) and antigen (N)
                                   = sum of the antibody-chain x antigen-chain pair energies

Prints Phase A (potency vs WT antigen, designs vs WT antibody) and Phase B (robustness:
mean + worst-case interaction energy across antigen mutants). More negative = tighter.

External tool: a FoldX 5 binary (free academic license). No ChimeraX/PDBFixer/PRODIGY.
Verify on one antibody first with --limit 1; FoldX command/output details are best-effort
and worth eyeballing on the first pair.

Usage
-----
    python3 scripts/nos_diffusion/validate_designs_foldx.py \
        --designs $RESULTS/validation/1A14_HLN_designs.json \
        --pdb /data/.../complex_structure/1A14_HLN.pdb \
        --workdir $RESULTS/validation_foldx/1A14_HLN \
        --foldx ${VAAD_ROOT}/tools/foldx/foldx_20251231 [--limit 1]
"""
from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import subprocess
import os as _os
_VAAD_DEFAULT = _os.path.expanduser("~/vaad-data")
VAAD_ROOT = _os.environ.get("VAAD_ROOT", _VAAD_DEFAULT)


_MSE_FIX = ("HETATM", "ATOM  ")


def renumber_pdb(pdb_in, pdb_out):
    """Write a clean PDB with unique sequential resnums per chain (insertion codes removed,
    MSE->MET, only protein ATOM records). Returns map (chain, orig_resnum:int, icode:str) ->
    new_resnum:int, matching the (chain,resnum,icode) keys the design JSON uses."""
    cmap, counters, out_lines = {}, {}, []
    for line in open(pdb_in, errors="ignore"):
        rec = line[:6]
        if rec == "HETATM" and line[17:20] == "MSE":
            line = line.replace("HETATM", "ATOM  ").replace("MSE", "MET")
            rec = "ATOM  "
        if rec != "ATOM  ":
            continue
        chain = line[21]
        field = line[22:27].strip()
        if not field:
            continue
        if field[-1].isalpha():
            resnum, icode = int(field[:-1]), field[-1].upper()
        else:
            resnum, icode = int(field), ""
        key = (chain, resnum, icode)
        if key not in cmap:
            counters[chain] = counters.get(chain, 0) + 1
            cmap[key] = counters[chain]
        new = cmap[key]
        out_lines.append(f"{line[:22]}{new:>4d} {line[27:]}")
    with open(pdb_out, "w") as fh:
        fh.write("".join(out_lines))
    return cmap


def _link_molecules(cwd, foldx_bin):
    """FoldX needs its molecules/ dir (and rotabase.txt for v4) findable from cwd; symlink them in."""
    src_dir = os.path.dirname(os.path.abspath(foldx_bin))
    for name in ("molecules", "rotabase.txt"):
        src, dst = os.path.join(src_dir, name), os.path.join(cwd, name)
        if os.path.exists(src) and not os.path.exists(dst):
            try:
                os.symlink(src, dst)
            except OSError:
                pass


def foldx(cmd_args, cwd, foldx_bin):
    """Run a FoldX command in cwd; return (stdout+stderr)."""
    _link_molecules(cwd, foldx_bin)
    proc = subprocess.run([foldx_bin] + cmd_args, cwd=cwd, capture_output=True, text=True)
    return (proc.stdout or "") + "\n" + (proc.stderr or "")


def repair(pdb_name, cwd, foldx_bin):
    foldx(["--command=RepairPDB", f"--pdb={pdb_name}", "--output-dir=."], cwd, foldx_bin)
    return pdb_name[:-4] + "_Repair.pdb"


def build_model(repaired, muts, cwd, foldx_bin, tag):
    """Apply mutations (list of FoldX strings) to `repaired` in cwd; return mutant pdb name."""
    listfile = f"individual_list_{tag}.txt"
    with open(os.path.join(cwd, listfile), "w") as fh:
        fh.write(",".join(muts) + ";\n")
    out = foldx(["--command=BuildModel", f"--pdb={repaired}", f"--mutant-file={listfile}",
                 "--numberOfRuns=1", "--output-dir=."], cwd, foldx_bin)
    mutant = repaired[:-4] + "_1.pdb"
    if not os.path.exists(os.path.join(cwd, mutant)):
        raise RuntimeError(f"BuildModel produced no mutant ({mutant}); tail:\n"
                           + "\n".join(out.strip().splitlines()[-8:]))
    return mutant


def analyse_complex(pdb_name, ab_chains, ag_chains, cwd, foldx_bin):
    """AnalyseComplex; sum interaction energy over antibody-chain x antigen-chain pairs."""
    chains = ",".join(ab_chains + ag_chains)
    foldx(["--command=AnalyseComplex", f"--pdb={pdb_name}",
           f"--analyseComplexChains={chains}", "--output-dir=."], cwd, foldx_bin)
    inter = os.path.join(cwd, f"Interaction_{pdb_name[:-4]}_AC.fxout")
    if not os.path.exists(inter):
        # FoldX sometimes names it without the trailing _AC
        cands = [f for f in os.listdir(cwd) if f.startswith("Interaction_") and f.endswith(".fxout")]
        if not cands:
            raise RuntimeError("AnalyseComplex produced no Interaction_*.fxout")
        inter = os.path.join(cwd, sorted(cands)[-1])
    return _sum_interface_energy(inter, set(ab_chains), set(ag_chains))


def _sum_interface_energy(path, ab, ag):
    """Sum 'Interaction Energy' rows whose two groups straddle the antibody/antigen split."""
    rows = [ln.rstrip("\n").split("\t") for ln in open(path) if ln.strip()]
    header = next((r for r in rows if "Interaction Energy" in r), None)
    if header is None:
        raise RuntimeError(f"no 'Interaction Energy' column in {os.path.basename(path)}")
    g1, g2 = header.index("Group1"), header.index("Group2")
    ie = header.index("Interaction Energy")
    total, found = 0.0, False
    for r in rows:
        if len(r) <= ie or r is header or not r[g1].strip():
            continue
        a, b = r[g1].strip(), r[g2].strip()
        if ({a, b} & ab) and ({a, b} & ag) and not ({a, b} <= ab) and not ({a, b} <= ag):
            total += float(r[ie]); found = True
    if not found:
        raise RuntimeError(f"no antibody-antigen chain pair found in {os.path.basename(path)}")
    return total


def _to_foldx(chain, new_resnum, wt1, mut1):
    return f"{wt1}{chain}{new_resnum}{mut1}"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--designs", required=True)
    ap.add_argument("--pdb", required=True)
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--foldx", required=True, help="path to the FoldX binary")
    ap.add_argument("--limit", type=int, default=0, help="only the first N antibodies (smoke test)")
    ap.add_argument("--antigens", nargs="+", default=None,
                    help="only score these antigen-variant names (e.g. WT) — fast targeted run")
    args = ap.parse_args()

    spec = json.load(open(args.designs))
    os.makedirs(args.workdir, exist_ok=True)
    ab_chains, ag_chains = spec["antibody_chains"], spec["antigen_chains"]
    cdr_res, native_cdr = spec["cdr_residues"], spec["native_cdr"]

    # 1. renumber once (strip insertion codes), 2. RepairPDB once
    renum = os.path.join(args.workdir, "complex_renum.pdb")
    cmap = renumber_pdb(args.pdb, renum)
    print(f"renumbered {len(cmap)} residues -> {renum}", flush=True)
    repaired = repair("complex_renum.pdb", args.workdir, args.foldx)
    if not os.path.exists(os.path.join(args.workdir, repaired)):
        raise SystemExit(f"RepairPDB failed (no {repaired}) — check the FoldX binary/license")

    def newnum(chain, resnum, icode):
        return cmap.get((chain, resnum, (icode or "").upper()))

    antibodies = [{"name": "WT", "cdr_seq": native_cdr}] + spec["antibody_designs"]
    if args.limit:
        antibodies = antibodies[:args.limit]
    antigens = spec["antigen_variants"]
    if args.antigens:
        keep = set(args.antigens)
        antigens = [a for a in antigens if a["name"] in keep]
        print(f"antigen filter: {len(antigens)} of {len(spec['antigen_variants'])} variants "
              f"({', '.join(a['name'] for a in antigens)})", flush=True)

    # Incremental checkpoint so a wall-clock timeout salvages everything scored so far and a
    # resubmit resumes instead of recomputing — big panels (e.g. 7 ab x 49 antigen = 343 FoldX
    # pairs) can exceed a single job's wall. The final _report() stdout is unchanged.
    ckpt = os.path.join(args.workdir, "E_checkpoint.json")
    E = json.load(open(ckpt)) if os.path.exists(ckpt) else {}

    def _save_ckpt():
        tmp = ckpt + ".tmp"
        with open(tmp, "w") as fh:
            json.dump(E, fh)
        os.replace(tmp, ckpt)

    # rebuild the (cdr_seq, antigen) cache from checkpointed scores so resumed runs skip them
    cache = {}
    for ab in antibodies:
        for agn, v in E.get(ab["name"], {}).items():
            cache[(ab["cdr_seq"], agn)] = v
    if cache:
        print(f"resume: {len(cache)} (antibody,antigen) scores loaded from {ckpt}", flush=True)

    for ab in antibodies:
        E.setdefault(ab["name"], {})
        cdr_m = [(r["chain"], newnum(r["chain"], r["resnum"], r["icode"]), w, a)
                 for r, a, w in zip(cdr_res, ab["cdr_seq"], native_cdr) if a != w]
        for ag in antigens:
            key = (ab["cdr_seq"], ag["name"])
            if key in cache:
                E[ab["name"]][ag["name"]] = cache[key]
                _save_ckpt()
                continue
            ag_m = [(m["chain"], newnum(m["chain"], m["resnum"], m["icode"]), m["wt"], m["mut"])
                    for m in ag["mutations"]]
            allm = cdr_m + ag_m
            tag = f"{ab['name']}__{re.sub(r'[^A-Za-z0-9]', '_', ag['name'])}"
            print(f"[{tag}] {len(allm)} mutations", flush=True)
            val = None
            if any(m[1] is None for m in allm):
                print("    skip: a mutation residue did not map to the renumbered PDB")
            else:
                pair_dir = os.path.join(args.workdir, tag)
                os.makedirs(pair_dir, exist_ok=True)
                shutil.copy(os.path.join(args.workdir, repaired), pair_dir)
                try:
                    if allm:
                        muts = [_to_foldx(c, n, wt, mut) for c, n, wt, mut in allm]
                        scored = build_model(repaired, muts, pair_dir, args.foldx, tag)
                    else:
                        scored = repaired
                    val = analyse_complex(scored, ab_chains, ag_chains, pair_dir, args.foldx)
                except Exception as e:  # noqa: BLE001
                    print(f"    {type(e).__name__}: {e}")
            print(f"    interaction E = {val if val is not None else 'n/a'}", flush=True)
            cache[key] = val
            E[ab["name"]][ag["name"]] = val
            _save_ckpt()

    _report(antibodies, antigens, E)


def _fmt(v):
    return f"{v:8.2f}" if isinstance(v, float) else f"{'n/a':>8}"


def _report(antibodies, antigens, E):
    wt_ag = "WT"
    mut_names = [a["name"] for a in antigens if a["name"] != wt_ag]

    print("\n=== Phase A: potency (FoldX interaction E vs WT antigen; more negative = tighter) ===")
    print(f"{'antibody':<14} {'E(WT ag)':>10} {'vs WT ab':>10}")
    wt_e = E.get("WT", {}).get(wt_ag)
    for ab in antibodies:
        v = E[ab["name"]].get(wt_ag)
        d = (v - wt_e) if isinstance(v, float) and isinstance(wt_e, float) else None
        ds = f"{d:+8.2f}" if d is not None else f"{'n/a':>8}"
        print(f"{ab['name']:<14} {_fmt(v):>10} {ds:>10}")

    if mut_names:
        print("\n=== Phase B: robustness across antigen variants (interaction E; worst = least negative) ===")
        print(f"{'antibody':<14} {'mean':>8} {'worst':>8} {'WT-ag':>8} {'drop':>8}")
        for ab in antibodies:
            vals = [E[ab["name"]].get(n) for n in mut_names]
            vals = [v for v in vals if isinstance(v, float)]
            wv = E[ab["name"]].get(wt_ag)
            if not vals:
                print(f"{ab['name']:<14} {'n/a':>8}")
                continue
            mean, worst = sum(vals) / len(vals), max(vals)
            drop = (worst - wv) if isinstance(wv, float) else None
            ds = f"{drop:+8.2f}" if drop is not None else f"{'n/a':>8}"
            print(f"{ab['name']:<14} {mean:8.2f} {worst:8.2f} {_fmt(wv):>8} {ds:>8}")
        print("\n(robustness claim: cluster designs have a better — more negative — WORST-CASE "
              "interaction E across mutants, and a smaller drop from WT antigen, than the WT antibody)")


if __name__ == "__main__":
    main()
