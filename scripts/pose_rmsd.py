#!/usr/bin/env python3
"""Ligand RMSD between two HADDOCK models after Kabsch superposition on chain A.

HADDOCK writes each run in its own coordinate frame, so ligand coordinates must
NOT be compared directly (project CLAUDE.md). Superpose on the common receptor
(chain A CA atoms), then measure ligand heavy-atom RMSD.

Run with the haddock3 venv python (numpy required):
  source /home/yakuri/opt/haddock3/.venv/bin/activate
"""
import gzip
import sys
import numpy as np


def read(path):
    op = gzip.open if path.endswith(".gz") else open
    out = []
    with op(path, "rt") as f:
        for l in f:
            if not l.startswith(("ATOM", "HETATM")):
                continue
            out.append(dict(name=l[12:16].strip(), resname=l[17:20].strip(),
                            chain=l[21], resnum=int(l[22:26]),
                            el=(l[76:78].strip() or l[12:16].strip()[0]).upper(),
                            xyz=np.array([float(l[30:38]), float(l[38:46]),
                                          float(l[46:54])])))
    return out


def kabsch(P, Q):
    """Return f(x) mapping frame of P onto frame of Q."""
    cp, cq = P.mean(0), Q.mean(0)
    H = (P - cp).T @ (Q - cq)
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1.0, 1.0, d])
    R = Vt.T @ D @ U.T
    return lambda x: (R @ (np.asarray(x) - cp)) + cq


def compare(a, b, label=""):
    A, B = read(a), read(b)
    ca_a = {x["resnum"]: x["xyz"] for x in A if x["chain"] == "A" and x["name"] == "CA"}
    ca_b = {x["resnum"]: x["xyz"] for x in B if x["chain"] == "A" and x["name"] == "CA"}
    common = sorted(set(ca_a) & set(ca_b))
    P = np.array([ca_a[r] for r in common])
    Q = np.array([ca_b[r] for r in common])
    tf = kabsch(P, Q)
    Pt = np.array([tf(p) for p in P])
    fit = float(np.sqrt(((Pt - Q) ** 2).sum(1).mean()))

    la = [x for x in A if x["chain"] == "B" and x["el"] != "H"]
    lb = {x["name"]: x["xyz"] for x in B if x["chain"] == "B" and x["el"] != "H"}
    pairs = [(x["xyz"], lb[x["name"]]) for x in la if x["name"] in lb]
    if not pairs:
        print("  no matching ligand atoms")
        return
    LA = np.array([tf(p) for p, _ in pairs])
    LB = np.array([q for _, q in pairs])
    rms = float(np.sqrt(((LA - LB) ** 2).sum(1).mean()))
    dcen = float(np.linalg.norm(LA.mean(0) - LB.mean(0)))
    if label:
        print(label)
    print("  receptor CA fit RMSD : %.3f A  (%d residues)" % (fit, len(common)))
    print("  ligand heavy RMSD    : %.3f A  (%d atoms, symmetry-unaware)"
          % (rms, len(pairs)))
    print("  ligand centroid shift: %.3f A" % dcen)


if __name__ == "__main__":
    print("A: %s\nB: %s" % (sys.argv[1], sys.argv[2]))
    compare(sys.argv[1], sys.argv[2])
