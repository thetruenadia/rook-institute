# EEB-1 — Findings and evidence status
2026-11-14. This is a synthesis of the conversation, not a fresh verification of all archived experiments. Historical counts and outcomes below are reported, not independently re-executed here.

| Stage | Reported finding / status |
|---|---|
| Geography v0.2–0.10 | Terrain, runoff, depression filling, seasonal liquids and species bookkeeping explored; source/results not audited here. Do not use exact numerical predictions as verified. |
| E01–04 | Conditional reaction, diffusion and adsorption models; historical execution evidence incomplete. E04B withdrawn. |
| E05 | Heredity suspended; no demonstrated replication mechanism. |
| E06 | Balanced C2H2 + 3H2 → 2CH4; coupled transport can restrict potential chemical power. Reference energy 334 kJ/mol not local validated Gibbs energy. |
| E07 | Conditional H2 production/escape/transport balance; 0.1% atmospheric H2 not independently selectable from source and loss rates. |
| E08 | Partial saturation sometimes improves reactant access; weak catalysts and thick films can erase benefit. Geometries assumed. |
| E09 | Ideal capillary rise permits gas/liquid coexistence; isolated films may evaporate quickly. Actual film stability unresolved. |
| E10 | Ideal mixture can gain methane from atmosphere; replenishment may preserve dissolved acetylene; excludes nitrogen/nonideality. |
| E11 | Four related Peng–Robinson variants challenge 44:44:12 sea composition; revised ~40:53:7 candidate; approximate binary calibration only. Gas release quantities unvalidated. |
| E12 | Reported independent 21-point density benchmark failed (~18% mean error). Deep mixed liquids have larger sensible thermal inertia than shallow pools. |
| E13 | Reported held-out mean error fell from 19.5% to 4.1% with constant volume correction; small dataset, not independent new external validation; derivatives unvalidated. |
| E14 | Conditional methane-rich rain layering, sensitive to mixing; 100 m grid failed resolution screen. GitHub run directly inspected in prior conversation: two runs identical; 4/6 historical hashes matched, sweep.csv and profiles.csv did not. |
| E15 | Source committed, not executed; experiment programme paused by user. Do not report results. |

## Established external evidence versus EEB-1
Titan observations and laboratory work support hydrocarbon lakes, organic dunes, atmospheric photochemistry, nitrogen dissolution and possible gas release, molecular co-crystals and varied terrain. These establish analogues, not EEB-1 occurrences. See sources index.

## Analytical relations retained
g=GM/R²; escape speed=sqrt(2GM/R); capillary height=2γ cosθ/(Δρgr); sensible thermal change=Q/(ρcpd); conductive periodic skin depth=sqrt(κP/π); hydrostatic pressure increment≈ρgd. All require appropriate assumptions. Actual parameters and local conditions are not measured.

## Open issues
Atmospheric escape and methane replenishment; validated ternary phase equilibrium and density derivatives; geological material inventory/strength; wind climate and marine mixing; film stability and catalytic kinetics; nutrients and heredity. No biological mechanism verified.

## Primary archive
[E14 GitHub run](https://github.com/thetruenadia/rook-institute/actions/runs/43598742961), [source](Experiment-014/column.py), [scope](Experiment-014/GITHUB_REPRODUCTION.md). Workflow artifact retention is finite; preserve valuable raw data separately. [E15 source](Experiment-015/wind_mixing.py).
