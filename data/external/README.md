# External reference data

Digitized reference curves used as overlays by the `btmw` paper-figure pipeline.
Each file is a two-column whitespace-separated text table (columns documented
below and in each file's header comment). These files are **not** produced by
`btmw` and should not be modified.

All files were copied from the original paper-plot directory tree on the CUHK
Central Cluster.

## Active files

| File | Columns | Used by | Reference |
|------|---------|---------|-----------|
| `hellwing2016_smf_r50c.txt` | `log10(mu)`, `log10(N(>mu))` | `btmw hmf`, `btmw hmf --compare-vr`, `btmw hmf --resolution-study` | [Hellwing+ 2016](https://doi.org/10.1093/mnras/stw214) — COCO best-fit power law for subhalos within R50c |
| `cautun2014_svf_r100c.txt` | `log10(nu)`, `log10(N(>nu))` | `btmw hvf`, `btmw hvf --compare-vr`, `btmw hvf --resolution-study` | [Cautun+ 2014](https://doi.org/10.1093/mnras/stu1829) — Cautun fit for subhalos within R100c |
| `springel2008_nprofile_r200c.txt` | `log10(r/kpc)`, `log10(n(r)/<n>)` | `btmw hrf`, `btmw hrf --compare-vr`, `btmw hrf --resolution-study` | [Springel+ 2008](https://doi.org/10.1111/j.1365-2966.2008.14066.x) — Aquarius A1 subhalo radial profile (R200c ≈ 245 kpc) |
| `lovell2014_cmf_r200c.txt` | `r [kpc]`, `M_sub(<r)/M_tot(<r)` | `btmw cmf` | [Lovell+ 2014](https://doi.org/10.1093/mnras/stt2431) — cumulative substructure mass fraction (R200c ≈ 256.1 kpc) |
| `grand2021_rmax_r200c.txt` | `Vmax [km/s]`, `Rmax [kpc]` | `btmw rvsv`, `btmw rvsv --resolution-study` | [Grand & White 2021](https://doi.org/10.1093/mnras/staa3993) — Subhalos within R200c (referred to as "Robert result" in figure legends) |

## Files parked in `trash/`

The `trash/` subdirectory contains digitized COCO M-Vmax and R-Vmax curves that
are **not loaded by the current Python code path**. They are kept temporarily so
a full figure-reproduction run can confirm they are safe to delete permanently.

Note: the black M-Vmax reference line in `btmw mvsv` is an **analytic**
BolshoiP+MDPL power-law hard-coded in `src/btmw/figures/mvsv.py`; it is not
read from `data/external/`.

If you discover a discrepancy between any active file and the original
publication, please open an issue.

---

## BibTeX references

```bibtex
@article{10.1093/mnras/stw214,
    author  = {Hellwing, Wojciech A. and Frenk, Carlos S. and Cautun, Marius
               and Bose, Sownak and Helly, John and Jenkins, Adrian
               and Sawala, Till and Cytowski, Maciej},
    title   = {The Copernicus Complexio: a high-resolution view of the small-scale Universe},
    journal = {Monthly Notices of the Royal Astronomical Society},
    volume  = {457},
    number  = {4},
    pages   = {3492--3509},
    year    = {2016},
    doi     = {10.1093/mnras/stw214},
}

@article{10.1093/mnras/stu1829,
    author  = {Cautun, Marius and Hellwing, Wojciech A. and van de Weygaert, Rien
               and Frenk, Carlos S. and Jones, Bernard J. T. and Sawala, Till},
    title   = {Subhalo statistics of galactic haloes: beyond the resolution limit},
    journal = {Monthly Notices of the Royal Astronomical Society},
    volume  = {445},
    number  = {2},
    pages   = {1820--1835},
    year    = {2014},
    doi     = {10.1093/mnras/stu1829},
}

@article{10.1111/j.1365-2966.2008.14066.x,
    author  = {Springel, V. and Wang, J. and Vogelsberger, M. and Ludlow, A.
               and Jenkins, A. and Helmi, A. and Navarro, J. F.
               and Frenk, C. S. and White, S. D. M.},
    title   = {The Aquarius Project: the subhaloes of galactic haloes},
    journal = {Monthly Notices of the Royal Astronomical Society},
    volume  = {391},
    number  = {4},
    pages   = {1685--1711},
    year    = {2008},
    doi     = {10.1111/j.1365-2966.2008.14066.x},
}

@article{10.1093/mnras/stt2431,
    author  = {Lovell, Mark R. and Frenk, Carlos S. and Eke, Vincent R.
               and Jenkins, Adrian and Gao, Liang and Theuns, Tom},
    title   = {The properties of warm dark matter haloes},
    journal = {Monthly Notices of the Royal Astronomical Society},
    volume  = {439},
    number  = {1},
    pages   = {300--317},
    year    = {2014},
    doi     = {10.1093/mnras/stt2431},
}

@ARTICLE{2021MNRAS.501.3558G,
    author  = {{Grand}, Robert J.~J. and {White}, Simon D.~M.},
    title   = {Baryonic effects on the detectability of annihilation radiation
               from dark matter subhaloes around the Milky Way},
    journal = {Monthly Notices of the Royal Astronomical Society},
    volume  = {501},
    number  = {3},
    pages   = {3558--3567},
    year    = {2021},
    doi     = {10.1093/mnras/staa3993},
}
```
