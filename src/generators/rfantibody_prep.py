#!/usr/bin/env python3
"""Prepare one AACDB complex for the RFantibody chain.

RFantibody consumes HLT-format PDBs (chains renamed H/L/T, CDR loops annotated as REMARKs),
and its converter identifies CDRs by CHOTHIA residue number -- but AACDB structures carry
author numbering. So this step:

  1. ANARCI-renumbers the antibody chains to Chothia (heavy and light independently),
  2. writes a Chothia-numbered copy of the complex,
  3. runs RFantibody's own chothia2HLT.py to produce the HLT framework file,
  4. writes an antigen-only PDB (RFdiffusion's `--target`),
  5. derives hotspot residues = antigen residues within `--contact-cutoff` A of the native
     CDR-H3, which is what tells RFdiffusion WHERE on the antigen to bind. Without hotspots it
     would dock anywhere on the surface and the comparison against a native antibody that binds
     a specific epitope would be meaningless.

Also emits `meta.json` with the native CDR-H3 length, so the design-loops spec can pin H3 to
the native length and keep the comparison length-matched.

Runs under diffab_env (which has ANARCI) with DiffAb's bin on PATH (hmmscan).
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import os as _os
_VAAD_DEFAULT = _os.path.expanduser("~/vaad-data")
VAAD_ROOT = _os.environ.get("VAAD_ROOT", _VAAD_DEFAULT)

AA3to1 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E",
          "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F",
          "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"}
RFAB = f"{VAAD_ROOT}/tools/RFantibody"


def read_atoms(pdb):
    rows = []
    with open(pdb) as fh:
        for line in fh:
            if line.startswith(("ATOM", "HETATM")):
                rows.append(line.rstrip("\n"))
    return rows


def residues_of(rows, chain):
    """Ordered unique (resnum, icode, resname) for one chain, ATOM records only."""
    seen, out = set(), []
    for line in rows:
        if not line.startswith("ATOM") or line[21] != chain:
            continue
        key = (int(line[22:26]), line[26].strip())
        if key in seen:
            continue
        seen.add(key)
        out.append((key[0], key[1], line[17:20].strip()))
    return out


def anarci_chothia(seq, chain_type):
    """-> list of (chothia_resnum, icode) aligned to the INPUT sequence positions, None where
    ANARCI did not assign a number (outside the variable domain)."""
    from anarci import run_anarci
    res = run_anarci([("q", seq)], scheme="chothia", allow=set([chain_type]))
    numbering = res[1][0]
    if not numbering:
        return None
    dom, start, end = numbering[0][0], numbering[0][1], numbering[0][2]
    out = [None] * len(seq)
    i = start
    for (num, ins), aa in dom:
        if aa == "-":
            continue
        if i <= end:
            out[i] = (num, ins.strip())
            i += 1
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--template", required=True)
    ap.add_argument("--pdb-dir", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--contact-cutoff", type=float, default=5.0)
    ap.add_argument("--max-hotspots", type=int, default=6)
    args = ap.parse_args()

    tmpl = json.load(open(args.template))
    target = tmpl["target"]
    native_cdr = tmpl["native_cdr"]
    ab = list(tmpl["antibody_chains"])
    ag = list(tmpl["antigen_chains"])
    heavy = tmpl["cdr_residues"][0]["chain"]
    light = next((c for c in ab if c != heavy), None)
    os.makedirs(args.outdir, exist_ok=True)

    pdb = os.path.join(args.pdb_dir, f"{target}.pdb")
    rows = read_atoms(pdb)

    # --- 1. Chothia renumbering of the antibody chains -----------------------
    remap = {}                                        # (chain, oldnum, oldicode) -> (new, icode)
    for ch, ctype in [(heavy, "H"), (light, "K")]:
        if ch is None:
            continue
        res = residues_of(rows, ch)
        seq = "".join(AA3to1.get(r[2], "X") for r in res)
        numb = anarci_chothia(seq, ctype)
        if numb is None and ctype == "K":             # lambda light chains
            numb = anarci_chothia(seq, "L")
        if numb is None:
            raise SystemExit(f"[{target}] ANARCI failed on chain {ch}")
        for (old, oic, _rn), new in zip(res, numb):
            if new is not None:
                remap[(ch, old, oic)] = new

    # RFdiffusion builds a rigid frame per residue from N/CA/C, and dies with
    # "Non-positive determinant ... in rotation matrix" on any residue whose backbone is
    # incomplete or zero-filled. Drop those residues here rather than letting the diffusion
    # crash 20 minutes in.
    bb = {}
    for line in rows:
        if not line.startswith("ATOM"):
            continue
        atom = line[12:16].strip()
        if atom in ("N", "CA", "C"):
            key = (line[21], int(line[22:26]), line[26].strip())
            xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
            if xyz != (0.0, 0.0, 0.0):
                bb.setdefault(key, set()).add(atom)
    complete = {k for k, v in bb.items() if len(v) == 3}
    n_dropped = 0

    cho_path = os.path.join(args.outdir, f"{target}_chothia.pdb")
    kept_ab = 0
    with open(cho_path, "w") as fh:
        for line in rows:
            if not line.startswith("ATOM"):
                continue
            ch = line[21]
            if (ch, int(line[22:26]), line[26].strip()) not in complete:
                n_dropped += 1
                continue
            if ch in (heavy, light):
                key = (ch, int(line[22:26]), line[26].strip())
                if key not in remap:
                    continue                          # outside the variable domain -> dropped
                num, ic = remap[key]
                fh.write(f"{line[:22]}{num:>4d}{(ic or ' '):1s}{line[27:]}\n")
                kept_ab += 1
            elif ch in ag:
                fh.write(line + "\n")
    if n_dropped:
        print(f"[{target}] dropped {n_dropped} atoms on residues with incomplete backbone",
              flush=True)
    if kept_ab == 0:
        raise SystemExit(f"[{target}] no antibody atoms survived Chothia renumbering")

    # --- 2. HLT framework via RFantibody's own converter ---------------------
    hlt_path = os.path.join(args.outdir, f"{target}_hlt.pdb")
    cmd = [sys.executable, os.path.join(RFAB, "scripts/util/chothia2HLT.py"), cho_path,
           "-H", heavy, "-T", ",".join(ag), "-o", hlt_path]
    if light:
        cmd += ["-L", light]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0 or not os.path.exists(hlt_path):
        sys.stderr.write(r.stdout[-1500:] + "\n" + r.stderr[-1500:] + "\n")
        raise SystemExit(f"[{target}] chothia2HLT failed")

    # --- 3. antigen-only target PDB -----------------------------------------
    tgt_path = os.path.join(args.outdir, f"{target}_antigen.pdb")
    with open(tgt_path, "w") as fh:
        for line in read_atoms(hlt_path):
            if line.startswith("ATOM") and line[21] == "T":
                fh.write(line + "\n")

    # --- 4. hotspots = antigen residues contacting the native CDR-H3 ---------
    hl = [l for l in read_atoms(hlt_path) if l.startswith("ATOM")]
    def xyz(l):
        return (float(l[30:38]), float(l[38:46]), float(l[46:54]))
    # The epitope is defined by the WHOLE Fv, not by CDR-H3 alone: the framework and the other
    # CDRs are kept native here and do most of the binding on many targets. (On 4NZT_HLM the
    # native CDR-H3's closest approach to the antigen is 4.9 A -- an H3-only cutoff finds a
    # single residue and would point RFdiffusion at essentially nothing.) H3 proximity is still
    # recorded as a diagnostic, since it says whether redesigning H3 can matter at all here.
    fv = [l for l in hl if l[21] in ("H", "L")]
    h3 = [l for l in hl if l[21] == "H" and 95 <= int(l[22:26]) <= 102]
    tg = [l for l in hl if l[21] == "T"]
    cut2 = args.contact_cutoff ** 2

    def contact_map(probe):
        out = {}
        for a in tg:
            ax, ay, az = xyz(a)
            for b in probe:
                bx, by, bz = xyz(b)
                d2 = (ax - bx) ** 2 + (ay - by) ** 2 + (az - bz) ** 2
                if d2 <= cut2:
                    k = int(a[22:26])
                    out[k] = min(out.get(k, 1e9), d2)
                    break
        return out

    contacts = contact_map(fv)
    h3_contacts = contact_map(h3)
    h3_min = min((min((((float(a[30:38]) - float(b[30:38])) ** 2 +
                        (float(a[38:46]) - float(b[38:46])) ** 2 +
                        (float(a[46:54]) - float(b[46:54])) ** 2) for b in h3), default=1e9)
                  for a in tg), default=1e9) ** 0.5 if h3 and tg else None
    hot = sorted(contacts, key=lambda k: contacts[k])[:args.max_hotspots]
    hotspots = ",".join(f"T{h}" for h in sorted(hot))

    meta = {"target": target, "heavy": heavy, "light": light, "antigen": ag,
            "native_cdr": native_cdr, "h3_len": len(native_cdr), "hotspots": hotspots,
            "n_contacts": len(contacts), "n_h3_contacts": len(h3_contacts),
            "h3_min_dist": (round(h3_min, 2) if h3_min is not None else None),
            "hlt": hlt_path, "antigen_pdb": tgt_path}
    with open(os.path.join(args.outdir, "meta.json"), "w") as fh:
        json.dump(meta, fh, indent=2)
    print(f"[{target}] HLT ok; H3 len {len(native_cdr)}; hotspots {hotspots or '(none)'}; "
          f"epitope {len(contacts)} res; H3 contacts {len(h3_contacts)} "
          f"(min dist {h3_min:.1f} A)" if h3_min is not None else "", flush=True)


if __name__ == "__main__":
    main()
