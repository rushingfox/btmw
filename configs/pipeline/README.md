# Pipeline configuration files

Configuration files used to run the nine zoom-in simulations analysed in this
paper and their halo-finding and post-processing pipeline.  Each file is
preserved verbatim from the corresponding simulation directory.

`fiducial_m12i/` — the three fiducial-resolution m12i runs (`m12i_cdmo`,
`m12i_btps_deep`, `m12i_btps_soft`).

`fiducial_m12f/` — the three fiducial-resolution m12f runs (`m12f_cdmo`,
`m12f_btps_deep`, `m12f_btps_soft`).

`highest_m12i/` — the three high-resolution m12i runs used for the resolution
study (`m12i_cdmo_highest`, `m12i_btps_deep_highest`,
`m12i_btps_soft_highest`).

Within each folder, the SWIFT parameter files of the three power-spectrum
models differ only in the output basename and initial-condition filename.
Relative to `fiducial_m12i/`, the SWIFT parameters of the other folders
differ as follows:

- `fiducial_m12f/`: `a_begin` is 0.01 (instead of 0.00990099),
  `mesh_side_length` is 458 (instead of 414), and `time_first` (only used by
  non-cosmological runs) is additionally listed.
- `highest_m12i/`: `mesh_side_length` is 826 and the DM softening length is
  halved (`comoving_DM_softening` and `max_physical_DM_softening` are
  0.00010375976 Mpc instead of 0.00020751953 Mpc).

`hbtplus.conf`, `soap_hbt.yml` and `soap_vr.yml` are kept once per folder
(from the `cdmo` run); the files of the other runs differ only in the
simulation name in their paths and, for m12f, in the final snapshot index
(49 instead of 50).  `soap_hbt.yml` and `soap_vr.yml` differ only in the
`HaloFinder` block (two lines); VELOCIraptor-based SOAP catalogues are only
used for m12i, so `soap_vr.yml` is only included there.
`velociraptor.cfg` is byte-identical across all nine runs (verified by md5),
so only one copy is kept.

| File | Tool | Purpose |
|------|------|---------|
| `swift_simulation_<sim>.yml` | [SWIFT](https://swift.strw.leidenuniv.nl) | SWIFT N-body simulation parameters (cosmology, softening, ICs, etc.) |
| `hbtplus.conf` | [HBT-HERONS](https://github.com/SWIFTSIM/HBT-HERONS) | HBT-HERONS subhalo finder configuration |
| `velociraptor.cfg` | [VELOCIraptor](https://github.com/pelahi/VELOCIraptor-STF) | VELOCIraptor halo/subhalo finder configuration |
| `soap_hbt.yml` | [SOAP](https://github.com/SWIFTSIM/SOAP) | SOAP post-processing config using HBT-HERONS catalogues |
| `soap_vr.yml` | [SOAP](https://github.com/SWIFTSIM/SOAP) | SOAP post-processing config using VELOCIraptor catalogues |
