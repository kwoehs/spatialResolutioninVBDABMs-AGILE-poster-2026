# spatialResolutioninVBDABMs-AGILE-poster-2026

## DengueABM — Spatial Resolution

An agent-based model (ABM) exploring the impact of **grid cell size on dengue fever transmission** in a small Central European community. The model is implemented in [GAMA](https://gama-platform.org/) and designed to be driven from [OpenMOLE](https://openmole.org/) for distributed parameter sweeps.

> **Not a prediction model.** This is a prototype ABM used to study how spatial resolution (raster cell size) shapes outbreak dynamics — one facet of the Modifiable Areal Unit Problem (MAUP) in vector-borne disease modelling.

## Purpose

- **Research question**: How does the spatial resolution of the environmental input (NDVI-derived risk raster) affect the emergent transmission patterns of dengue fever in an area where *Aedes albopictus* is newly established?
- **Context**: Part of a PhD dissertation on spatial aspects of vector-borne disease agent-based models (focus: dengue / *Aedes albopictus*) at Z_GIS, University of Salzburg.
- **Study area**: Elsbethen, a small community south of Salzburg, Austria. Chosen because *Aedes albopictus* was first registered in the Salzburg area in 2023 (AGES, 2024).
- **Based on**: `salzburgDengueABM-v1` (Wöhs, 2025, master thesis).

## Repository layout

```
spatialResolutioninVBDABMs-AGILE-poster-2026/
├── README.md                              ← this file
├── ODD.md                                 ← ODD protocol description (Grimm et al. 2020)
├── models/
│   ├── dengueABM_spatialResolution_v2.gaml    ← GAMA model
│   └── dengueABM_spatialResolution_v2.oms     ← OpenMOLE script for direct sampling
├── includes/                              ← input rasters and shapefiles
└── results/                               ← transmission CSV exports (created at runtime)
```

## How to run

### GAMA (single run / GUI)
1. Open `models/dengueABM_spatialResolution_v2.gaml` in GAMA.
2. Launch experiment `singleSim` for an interactive run with map, chart, and transmission-location display.
3. To swap resolutions, uncomment one of the `raster_path` options in the `global` block (or set it via the `raster_path` parameter in the GUI parameters panel).

### GAMA (local batch)
- Experiment `simulationBatch` runs 100 Monte Carlo replicates in parallel across 16 cores and produces a histogram of outcomes.

### OpenMOLE (distributed sweep)
See `models/dengueABM_spatialResolution_v2.oms`. The script performs a direct sampling across the six raster resolutions with replication per resolution, writing results to `results/resolution_sweep.csv`.

## v2 changelog (since master thesis v1)

### Mosquito behaviour
- **EIP (extrinsic incubation period)** added: mosquitoes do not recover (virus dies with them).
- **Lifespan**: mosquito agents die after 45 days.
- **Home range**: 250 m radius centred on birth location (Vavassori, Saddler & Müller 2019); permanent relocation on breakout (`proba_dispersal`).
- **Edge-effect fix**: removed `elsbethen_geom` clipping from the mosquito home range in `reflex move` / `reflex hunt` — unclipped circle prevents boundary truncation.
- **Probabilistic biting**: `proba_mosquito_bite` (~10 %) gates the per-step `make_mosquito_agent` trigger, roughly one bite per 3–4 days.
- **Gonotrophic cycle cooldown**: mosquitoes cannot feed again for `gonotrophic_cycle` days (default 3) after a bite; `reflex hunt` guarded by the cooldown, `hours_since_last_feed` reset in `reflex infect`, `my_human` cleared after bite.

### Human behaviour
- **Edge-effect fix in `reflex move`**: humans wander freely within the wider `elsbethen_geom` instead of being hard-clipped to `settlement_geom` (which caused clustering along settlement boundaries). Soft preference via `settlement_geom covers self.location` with probabilistic pull-back (`flip(0.6)`) when outside — allows occasional longer excursions while keeping humans centred on where they live.
- **Movement branches consolidated**: three duplicate SIR-dependent `wander` branches collapsed into one (health status did not actually change movement).

### Transmission logic
- **Transmission location tracking**: `list<point> transmission_locations` populated in `reflex infect`, displayed as a graphics layer and exported to CSV for post-hoc KDE in ArcGIS Pro.
- **`riskArea` bug fix**: flag now resets to `false` each step, correctly reflecting current position.
- **`nb_localInfections` fix**: removed circular no-op update; now incremented in `reflex infect` on each mosquito-to-human transmission.
- **`reflex infect` bug fixes**:
  - `healthstatus` comparison corrected from ambiguous `= "I" or "R"` to explicit `= "I" or self.healthstatus = "R"`.
  - `proba_getinfectious_human` flip moved *inside* the `if self.is_infected` block; previously independent flips produced impossible `is_infectious=T, is_infected=F` states and a mismatch between the diagram count and the transmission map.
- **`N` state dropped**: asymptomatic infections now go straight to `R` (previously they were stuck forever in `updateSIR`); `S` branch guarded with `is_recovered = false` so asymptomatic-R humans (`infectionDay=0`) are not reset to S.
- **End-of-simulation flag fix**: `end_simulation_noinfection` no longer unconditionally sets `noLocalInfections=true` / `localInfections=false`; both end reflexes partition runs based on actual `nb_localInfections` count.

### Output and tracking
- **Disease spread chart**: `nb_humans_infected` now tracks `is_infectious` (I state only).
- **`time_of_peak` + `peak_infected`**: `reflex track_peak` records the cycle and count of the epidemic peak; monitored in the experiment outputs.
- **CSV export**: transmission locations saved per run with `CRS_transform` to ETRS89 LAEA (EPSG:3035), written at the end of the simulation.
- **Standalone `transmission_map` display** for results output.

### Simulation control
- **Extended duration**: simulation runs for 540 cycles (180 days / 6 months) instead of 9 days; `end_season` and `end_simulation_noinfection` reflexes both save the CSV and flip `sim_over`.
- **Batch experiment**: `repeat 100, parallel 16` (16 cores) for Monte Carlo runs.

### OpenMOLE preparation
- **`raster_path` parameterized**: global string allows swapping 10 m / 50 m / 100 m / 500 m / 1000 m / single-cell resolutions without code changes.
- **Auto-sized grid**: grid declaration now uses `file: grid_file(raster_path)`; dimensions and `grid_value` auto-populated; manual matrix copy removed from `init`.
- **Minimal `openmole` experiment**: `type: gui`, headless (no displays); exposes `raster_path` as a parameter and monitors `Time of peak` / `Peak infected` / `Local infections` for OpenMOLE to read. `until: sim_over or cycle >= 540` allows early termination.
- **`.oms` script** created for `DirectSampling` across the 6 raster paths with `Replication`, results hooked to CSV.

## Known limitations and TODO

### Justification still pending
- `proba_dispersal` value (currently 0.05 per 8 h step).
- `gonotrophic_cycle` value (currently 3 days; literature range 2–5 depending on temperature).
- `proba_infection_human` value (currently 0.1).

## Citation and contact

- **Author**: Katharina Wöhs, Department of Geoinformatics Z_GIS, University of Salzburg, katharina.woehs@plus.ac.at
- **Based on**: Wöhs, K. (2025). *Salzburg-Dengue ABM*. Master thesis, University of Salzburg.

## Licence:  
CC-BY-SA 4.0
