# ODD Protocol: DengueABM — Spatial Resolution (v2)

Model description following the ODD protocol (Grimm et al. 2020). This is the **spatial-resolution variant (v2)** of the Salzburg-Dengue ABM developed as part of Wöhs (2025), extended to study the effect of input raster resolution on emergent outbreak dynamics.

---

## 1. Purpose and patterns

The model is a **theoretical exposition** (in the sense of Grimm et al. 2020) of how the spatial resolution of the environmental input affects dengue transmission dynamics in a small community where *Aedes albopictus* is newly established. It is not a prediction model.

**Research question**: How does the grid cell size of the risk-area raster (10 m / 50 m / 100 m / 500 m / 1000 m / non-spatial baseline) shape the number, timing, and location of local dengue transmissions in Elsbethen, Austria?

**Patterns used for evaluation**:
- Presence/absence of local transmission over a 6-month simulation.
- Epidemic peak (time and height) when an outbreak occurs.
- Spatial distribution of transmission events (for post-hoc KDE analysis).

---

## 2. Entities, state variables, and scales

### Entities

| Entity    | Represented   | Notes |
|-----------|---------------|-------|
| **Humans**     | Explicitly    | `nb_humans = 3252`, one seeded as the index case (infectious traveller). |
| **Mosquitoes** | Implicitly until a transmission event, then explicitly. | Created on demand when an infectious human steps into a risk cell. |
| **Environment (grid)** | Explicit static raster. | Values from Sen2Cube (semantically enriched NDVI); cell size varies with the experiment (see Scales). |

### Human state variables

`is_susceptible`, `is_infected`, `is_infectious`, `is_recovered` (bool); `healthstatus` (string: S / I / R); `infectionDay` (int, hours since infection); `riskArea` (bool, whether current cell is classified as risk); `action_radius` (float, per-agent speed).

### Mosquito state variables

`is_infected`, `is_infectious` (bool); `eipDayVector` (int, hours in EIP); `age` (int, hours since creation, max 45 days); `hours_since_last_feed` (int, gonotrophic-cycle counter); `birth_location` (point, centre of 250 m home range); `my_human` (humans, current target).

### Global state variables (system level)

`nb_humans_infected`, `nb_mosquitos_infectious`, `nb_localInfections`, `peak_infected`, `time_of_peak`, `transmission_locations` (list of points), `sim_over` (bool).

### Scales

| Dimension | Value |
|-----------|-------|
| Spatial extent | Cadastral community Elsbethen (≈ 3.47 × 3.23 km), south of Salzburg, Austria. |
| Spatial resolution (experiment variable) | 10 m (native Sen2Cube) / 50 m / 100 m / 500 m / 1000 m / single-cell non-spatial baseline. |
| Temporal resolution | 1 time step = 8 hours. |
| Temporal extent | 540 cycles = 180 days (6 months). Early termination when no infected humans or mosquitoes remain and no mosquito is in EIP. |

---

## 3. Process overview and scheduling

At each 8 h time step, the following processes run in order:

1. **Humans** — `reflex move`: check current cell for risk classification; if infectious and on a risk cell, trigger `make_mosquito_agent` with probability `proba_mosquito_bite`. Then wander within `elsbethen_geom` with a soft preference for `settlement_geom` (probabilistic pull-back when outside).
2. **Mosquitoes** — `reflex lifespan`: increment `age` and `hours_since_last_feed`; die after 45 days.
3. **Mosquitoes** — `reflex move`: if within home range (250 m), wander inside; otherwise disperse (with probability `proba_dispersal`) and update `birth_location`.
4. **Mosquitoes** — `reflex hunt` (only if `hours_since_last_feed / 24 ≥ gonotrophic_cycle`): acquire a human target within perception distance; otherwise wander inside home range.
5. **Mosquitoes** — `reflex infect` (only if `my_human ≠ nil`): attempt human-to-mosquito or mosquito-to-human transmission; record transmission location; reset cooldown.
6. **Humans** — `reflex updateSIR`: advance S → I → R health-state transitions based on `infectionDay`.
7. **Global** — `reflex track_peak`: update `time_of_peak` and `peak_infected` if the current infectious count exceeds the recorded peak.
8. **Global** — end-of-season / no-infection exit checks (`reflex end_season`, `reflex end_simulation_noinfection`): flip `sim_over`, write CSV of transmission locations.

All agents are updated once per cycle. The order among humans (and among mosquitoes) within a step is GAMA's default (random).

---

## 4. Design concepts

- **Emergence**: local outbreak or extinction emerges from individual-level interactions; neither is imposed. Whether an introduced infectious traveller produces secondary infections depends on the spatial overlap between human trajectories and risk cells at the given resolution.
- **Adaptation & objectives**: mosquitoes pursue blood meals (`reflex hunt`) as their only objective; humans wander without goal-directed behaviour, with a soft preference for the settlement polygon.
- **Sensing**: mosquitoes perceive humans within a perception distance; humans sense the grid-value of their current cell to decide whether they are in a risk area.
- **Interaction**: humans trigger mosquito creation when in a risk cell; mosquitoes interact with humans during feeding (bidirectional transmission).
- **Stochasticity**: all infection events, biting events, dispersal, and human/mosquito movement decisions are stochastic; parameter values from literature.
- **Observation**: per-step counts of infectious humans and mosquitoes; cumulative `nb_localInfections`; peak tracking; list of transmission locations for spatial post-processing (KDE in ArcGIS Pro). OpenMOLE batch experiment reads `Time of peak`, `Peak infected`, `Local infections` as outputs.

### What changed vs. master-thesis v1

- **EIP** for mosquitoes (v1 had none).
- **Mosquito lifespan** (45 days) and **home range** (250 m circle).
- **Gonotrophic cycle cooldown** (3 days between bites).
- **Asymptomatic N-state dropped**: asymptomatic infections now go directly to R.
- **Simulation extended** from 9 days (v1) to 180 days (v2).
- **`raster_path` parameterized** so the grid cell size can be swept via OpenMOLE.
- **Auto-sized grid** from the raster file (no manual matrix copy).
- **Transmission locations** are now recorded and exported for spatial analysis.

---

## 5. Initialisation

- Spatial environment is loaded from a shapefile (`elsbethen_shp` for the cadastral community, `elsbethen_settlement_shp` for the settlement polygon) and a raster (`raster_path`, variable across experiments) that defines the risk layer.
- `nb_humans = 3252` human agents are created, placed at random locations inside `settlement_geom`. One is selected as the infectious index case (`is_infected = true`, `is_infectious = true`, `healthstatus = "I"`).
- Mosquitoes are **not** created at initialisation; they are spawned on demand when an infectious human steps into a risk cell.
- `nb_localInfections = 0`, `localInfections = false`, `sim_over = false`, `peak_infected = 0`, `time_of_peak = 0`.

---

## 6. Input data

| Input | Description | Source |
|-------|-------------|--------|
| `raster_path` | Sen2Cube NDVI-derived risk layer. Default: `masked_map_greenness_id79229-Elsbethen.tiff` (10 m). Additional aggregated variants for 50 m / 100 m / 500 m / 1000 m / single-cell baseline in `/includes`. | Sudmanns et al. (2021) |
| `elsbethen_shp` | Cadastral community boundary (movement bound for all agents). | BEV (2023) |
| `elsbethen_settlement_shp` | Permanent settlement area (soft preference for humans). | Land Salzburg (2025) |
| `proba_infection_vector = 0.9` | Probability a mosquito gets infected when biting an infectious human (DENV-1, *Ae. albopictus*, 28 °C). | Bellone et al. (2020) |
| `proba_getinfectious_vector = 0.37` | Transmission efficiency: virus survives in mosquito to make it infectious. | Bellone et al. (2020) |
| `proba_infection_human = 0.1` | Probability a human gets infected when bitten by an infectious mosquito. | Model assumption |
| `proba_getinfectious_human = 0.5` | Probability an infected human becomes symptomatic/infectious. | Assumption after Asish et al. (2023) |
| `proba_mosquito_bite = 0.1` | Probability of mosquito creation per 8 h step when an infectious human is on a risk cell. | Model assumption (≈ 1 bite per 3–4 days) |
| `proba_dispersal = 0.05` | Probability of mosquito leaving its home range per 8 h step. | Model assumption (justification pending) |
| `dispersal_radius = 250 m` | Home range radius. | Vavassori, Saddler & Müller (2019) |
| `eip_mosquito = 8 days` | Extrinsic incubation period. | Dengue literature (8–12 d at optimal temperature) |
| `gonotrophic_cycle = 3 days` | Cooldown between mosquito blood meals. | Literature (2–5 d depending on temperature) |
| `humans_min_speed / humans_max_speed` | 0 / 200 m/h | Model assumption |
| `vector_min_speed / vector_max_speed` | 0 / 50 m/h | Model assumption |

---

## 7. Submodels

- **`make_mosquito_agent`** (humans): creates one mosquito at the human's location with `is_infected = flip(proba_infection_vector)`. The mosquito enters EIP; it does *not* become infectious immediately.
- **`reflex infect`** (mosquitoes): if the mosquito has a human target and is infectious and the target is susceptible, the target is infected with probability `proba_infection_human`; on success, `nb_localInfections` is incremented and the transmission location is recorded. With probability `proba_getinfectious_human` the new infection becomes symptomatic (→ I), otherwise it is asymptomatic and the human is marked R immediately (no stuck N state). After the attempt, `hours_since_last_feed` is reset and `my_human` is cleared.
- **`reflex updateSIR`** (humans): advances the S / I / R health-state transitions based on `infectionDay`. After the equivalent of 8 days in state I, the human recovers.
- **`reflex lifespan`** (mosquitoes): increments `age` and `hours_since_last_feed`; mosquitoes die after 45 days.
- **`reflex track_peak`** (global): records the cycle and count at which `nb_humans_infected` reaches its maximum.
- **End-of-simulation reflexes** (`end_season`, `end_simulation_noinfection`): exit the run at 540 cycles or when no infectious agents (and no mosquitoes in EIP) remain, write the transmission CSV, and flip `sim_over`.

---

## References

- Grimm, V., Railsback, S. F., Vincenot, C. E., Berger, U., et al. (2020). *The ODD Protocol for Describing Agent-Based and Other Simulation Models: A Second Update to Improve Clarity, Replication, and Structural Realism*. JASSS 23(2), 7.
- Vavassori, L., Saddler, A., & Müller, P. (2019). *Active dispersal of Aedes albopictus*. Parasites & Vectors 12, 583.
- Bellone, R., et al. (2020). *The role of temperature in shaping mosquito-borne viruses transmission*. Frontiers in Microbiology 11.
- Sudmanns, M., et al. (2021). *Semantic Earth Observation Data Cubes (Sen2Cube)*.
- Wöhs, K. (2025). *Salzburg-Dengue ABM*. Master thesis, University of Salzburg.
