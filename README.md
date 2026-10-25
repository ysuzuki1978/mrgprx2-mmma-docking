# MRGPRX2 docking data — compound 48/80 and monomethyl methacrylate

Data, analysis code and Supplementary Tables for the Letter

> **Two routes into one pocket: compound 48/80 engages MRGPRX2 through charge and
> monomethyl methacrylate through van der Waals contact**

Each release of this repository is archived on Zenodo and gets its own DOI.
**Cite the DOI of the version you used**, not this repository URL.

## What is here

| | |
|---|---|
| `data/runs/` | every HADDOCK3 run the Letter cites (97 runs): the configuration, the per-cluster scores, the cluster assignment of all refined models, and the four models of the rank-1 cluster |
| `data/raw/` | the receptor, the variant receptors and the ligands as docked |
| `data/params/` | ACPYPE parameter and topology files for every ligand |
| `data/restraints/` | the restraint files of the restrained reference run |
| `data/tables/` | Supplementary Tables S1–S9 as TSV |
| `cache/` | `7vv6.cif` and `7vdm.cif`, read by Table S9 |
| `scripts/` | the analysis code |
| `protocol/` | the protocol written before the runs were submitted — **see the correction at the top of section 16 of `PREREG_TRIMER.md`** |

The runs are 18 conditions × 5 replicate seeds for the D184 substitution series,
plus 7 runs used only in the methods: four unrestrained replicates, the restrained
reference run, its Gasteiger-charge counterpart, and an antagonist control.

`cluster_<N>_model_<M>` files are named by cluster **rank**, not by cluster id.
Where the two differ, Table S1 gives both.

## Regenerating the numbers

With Python 3 and NumPy, from the top of this repository:

```sh
python3 scripts/make_ladder_table.py       # T9 = Table S3
python3 scripts/make_si_tables.py          # S1, S2, S4, S5, S7, S9
python3 scripts/make_s6_hbond_geometry.py  # S6
python3 scripts/make_s8_pose_rmsd.py       # S8   (needs NumPy)
python3 scripts/verify.py                  # checks every quoted number (46 items)
```

All nine tables are reproduced byte for byte, and `verify.py` exits non-zero if any
value quoted in the Letter disagrees with them.

`scripts/make_tables.py` regenerates the earlier tables T1–T8. Those draw on runs from
the first campaign that are not part of this deposit, so T1–T8 ship as generated;
`verify.py` reads them for the checks the two versions share.

HADDOCK3 (2025.11.0) is required only to repeat the docking itself.

## Notes

The HADDOCK3 configuration files contain the absolute paths of the machines the runs
were executed on. They are left exactly as they were used, because their value is that
they are the files that produced the runs.

`MANIFEST.sha256` lists a checksum for every file.

## Funding

Japan Society for the Promotion of Science (JSPS), Grant-in-Aid for Scientific Research,
KAKENHI grant number 21K16549.

## Licence

CC BY 4.0. See `LICENSE`.
