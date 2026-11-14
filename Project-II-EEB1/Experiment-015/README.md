# EEB-1 Experiment 15: wind and waves
Purpose: screen whether prescribed wind could provide enough energy to homogenise a two-layer liquid, and estimate capillary-gravity wave dispersion, not growth or wave height.

E14 GitHub run 43598742961: two executions agreed byte-for-byte; historical timeseries, timestep_check, spatial_check and results hashes matched; sweep and profiles did not. Exact archival reproduction remains unresolved. E15 does not use these numerical profiles as verified physical inputs. Density contrasts, depths, drag and mixing coefficients are independent sensitivity assumptions. Gravity 1.33 m/s², liquid density 560 kg/m³, atmosphere density 6 kg/m³ are fictional reference choices. Surface tensions span pure-component reference values, not measured mixture properties.

Two-layer potential energy per unit area: g*h1*h2*(rho_bottom-rho_top)/2. Wind stress = rho_air*Cd*U10²; liquid-side friction velocity = sqrt(stress/rho_liquid); prescribed mixing power = m*stress*u_star. The coefficient m is uncertain and transferred only as an exploratory parameter from terrestrial work. Sufficient integrated energy is necessary within this budget, NOT sufficient to guarantee mixing. No actual wind circulation, entrainment, wave growth, heat/phase changes, lateral transport, or sediment erosion model. Dispersion uses the inviscid single-interface finite-depth formula ignoring gas inertia; minimum phase speed is not the wind threshold for wave growth.

External evidence:
- Hayes et al. (2013), Icarus 225, 403–412: https://www.sciencedirect.com/science/article/pii/S0019103513001644
- Reichl et al. (2022), appendix B: https://doi.org/10.1029/2021JC018140
- Detelich et al. (2026), Titan-specific wave model: https://doi.org/10.1029/2026JE009693

Run: python3 wind_mixing.py (standard library, no randomness). Outputs mixing.csv, waves.csv, results.json. GitHub workflow runs twice and compares hashes; logs and datasets uploaded as artifacts. Report execution only after reviewing logs. Negative and positive energy-screen results are conditional mathematical results, not physical validation. No new biological canon.
