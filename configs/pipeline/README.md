# Pipeline configuration files

Representative configuration files used to run the halo-finding and post-processing
pipeline for the `m12i_cdmo` simulation (the fiducial CDMO run).  All other
simulations share the same settings; only the snapshot file prefix and ICs path
differ in `swift_simulation.yml`.

| File | Tool | Purpose |
|------|------|---------|
| `swift_simulation.yml` | [SWIFT](https://swift.strw.leidenuniv.nl) | SWIFT N-body simulation parameters (cosmology, softening, ICs, etc.) |
| `hbtplus.conf` | [HBT-HERONS](https://github.com/SWIFTSIM/HBT-HERONS) | HBT-HERONS subhalo finder configuration |
| `velociraptor.cfg` | [VELOCIraptor](https://github.com/pelahi/VELOCIraptor-STF) | VELOCIraptor halo/subhalo finder configuration |
| `soap_hbt.yml` | [SOAP](https://github.com/SWIFTSIM/SOAP) | SOAP post-processing config using HBT-HERONS catalogues |
| `soap_vr.yml` | [SOAP](https://github.com/SWIFTSIM/SOAP) | SOAP post-processing config using VELOCIraptor catalogues |

`soap_hbt.yml` and `soap_vr.yml` differ only in the `HaloFinder` block
(two lines); both are included for completeness.
