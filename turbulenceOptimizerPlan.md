# Turbulence Topology Optimization Plan for `laminarOptimizer`

## 1) Goal and scope

Transition `laminarOptimizer` to a turbulent-flow topology-optimized channel design code capable of handling internal flows in the 1-3 m/s inlet velocity range using a robust, validated formulation and numerically stable continuation strategy.

## 2) Literature-backed design baseline

Use the turbulent literature available in:
- Sun, Y.; Hao, M.; Wang, Z. *Topology Optimization of Turbulent Flow Cooling Structures Based on the k- Model*, Entropy, 2023, 25, 1299, DOI: 10.3390/e25091299.
- Li, Yijun. *Design and Topology Optimization of Heat Sinks for the Cooling of Electronic Devices with Multiple Heat Sources*, Ph.D. thesis, University of Nantes/HAL, 2023 (HAL: tel-04059777).
- Kaikow, H. et al. *Thermal Management for Electrification in Aircraft Engines: Optimization of Coolant System*, GT 2022 conference paper, DOI: 10.1115/GT2022-82538.

Observed validated principles:
- Use variable-density topology parameterization with PDE filtering and projection continuation.
- Couple design variable to fluid properties in RANS momentum + turbulence transport equations.
- Optimize multi-objective thermal/hydraulic metrics, then re-evaluate selected designs in higher-fidelity CFD.
- Cover operating envelope by running separate single-velocity optimization cases.
- Keep turbulent boundary layer representation pragmatic; dense local channels should be checked against mesh and model limitations.

## 3) Target mathematical formulation

### 3.1 Governing flow model in optimization
- Use steady RANS in all optimization loops.
- Base model: `k-epsilon` or `k-omega SST` as practical options.
  - Prefer `k-epsilon` for continuity with the surveyed 2023 cold-plate formulation.
  - Use `k-omega SST` where near-wall behavior must be captured better during final validation.

### 3.2 Material interpolation
- Design density `rho_tilde` in [0,1] maps to flow capacity/permeability and effective conductivity terms.
- Use SIMP-like power law with filtering + projection:
  - `rho_eff = H(rho_tilde; beta, eta)` (continuation on beta for crispness)
  - Penalize momentum resistance and turbulence production/dissipation source terms in low-density regions.

### 3.3 Turbulence equation coupling strategy
- First-stage optimizer:
  - Keep turbulence equations solved in primal flow updates.
  - Use frozen-turbulence adjoint approximation initially:
    - Use `nu_eff` in adjoint momentum operator where feasible.
    - Do not fully differentiate turbulence model in first pass.
- Second-stage enhancement:
  - Differentiate turbulence source penalties only if residual stability and objective gain justify added complexity.

### 3.4 Objective functions
- Core combined objectives:
  - Minimize thermal metric (mean or max wall/solid temperature).
  - Minimize pressure-drop or pumping-power metric.
- Velocity-envelope coverage strategy:
  - Run one inlet velocity per optimization case.
  - Use a case set spanning `U_in in [1, 3] m/s` (for example: 1, 2, 3 m/s).
  - Compare selected designs across the case set in post-processing.

## 4) Numerical robustness and stability plan

1. Variable clipping and feasibility
   - Enforce `rho_tilde in [rho_min, 1]`, with `rho_min` between 1e-3 and 1e-2.

2. Solver stabilization
   - Start from laminar-converged design/flow fields as lower-speed seeds.
   - Increase turbulence nonlinearity gradually with robust under-relaxation.
   - Use bounded transport schemes for k, epsilon and turbulent-viscosity terms.

3. Filtering and projection continuation
   - Keep existing Helmholtz filtering and projection strategy.
   - Increase projection sharpness only in late iterations.
   - Tune filter radius and projection beta independently.

4. Sensitivity stabilization
   - Apply chain rule through filtered density and projected density.
   - Use move limits and periodic asymmetry checks to suppress checkerboarding.

5. Mesh and BC consistency
   - Ensure turbulence field BCs are consistent with initial and inlet conditions.
   - Match near-wall resolution to chosen turbulence model and check y+ behavior for selected test cases.

## 5) Concrete code implementation phases in repository

### Phase 1 — Enable turbulent forward model
- Update case dictionaries:
  - `laminarOptimizer/app/constant/turbulenceProperties`: activate RANS model and turbulence transport settings.
  - `laminarOptimizer/app/0/`: add required turbulence fields (k/epsilon or k-omega equivalent).
  - `laminarOptimizer/app/system/fvSchemes` and `fvSolution`: add robust discretization and solvers for turbulence equations.
- Confirm `readTransportProperties.H` supports turbulent viscosity setup used by OpenFOAM.
- Run a forward-only smoke test on fixed geometry for 1-3 m/s.

### Phase 2 — Couple topology loop to turbulence resistance
- In `Primal_U.H` and related flow setup:
  - Ensure porosity/permeability mapping feeds both laminar and turbulent resistance terms.
  - Call `turbulence->correct()` at stable cadence in the SIMPLE loop.
- In `sensitivity.H` and `filter_chainrule.H`:
  - Ensure adjoint forcing is consistent with effective transport terms (`nu_eff`).

### Phase 3 — Adjoint extension
- Keep current adjoint structure, but replace laminar-like constants by design-aware effective transport terms.
- Add an explicit frozen-turbulence adjoint mode switch:
  - Turbulence state is fixed while accumulated sensitivities are computed.
- Validate descent-direction quality on reduced meshes first.

### Phase 4 — Continuation and case sweep
- Keep objective accumulation single-case (one inlet velocity per run).
- Add run scripts/case templates for inlet-velocity sweep up to 3 m/s.
- Add continuation schedule:
  - SIMP exponent, projection beta, and filter radius schedules.
- Keep normalization and constraints consistent within each single case.

### Phase 5 — Validation and trust checks
- Stage A: optimizer-level checks
  - Compare against laminar baseline and Darcy-equivalent surrogate.
  - Verify mesh and filter-ratio convergence.
- Stage B: high-fidelity verification
  - Re-run selected designs with stricter convergence and optionally SST k-omega.
  - Report pressure drop, hotspot temperature, thermal resistance, Nusselt proxy, and robust performance at each speed.

### Phase 6 — Production hardening and documentation
- Add case scripts for 1/2/3 m/s turbulent benchmark set.
- Document setup in reproducible files:
  - mesh, BCs, turbulence model, relaxation factors, and stopping criteria.
- Add caveats:
  - boundary-layer errors in coarse meshes,
  - model validity in tiny channels,
  - objective noise near discontinuities.

## 6) File change map (primary)

### Primary candidate files
- `laminarOptimizer/src/readTransportProperties.H`
- `laminarOptimizer/src/Primal_U.H`
- `laminarOptimizer/src/Primal_T.H`
- `laminarOptimizer/src/sensitivity.H`
- `laminarOptimizer/src/filter_chainrule.H`
- `laminarOptimizer/src/filter_x.H`
- `laminarOptimizer/src/update.H`
- `laminarOptimizer/app/system/fvSolution`
- `laminarOptimizer/app/system/fvSchemes`
- `laminarOptimizer/app/constant/turbulenceProperties`
- `laminarOptimizer/app/constant/transportProperties`
- `laminarOptimizer/app/constant/thermalProperties`
- `laminarOptimizer/app/0/{p,U,T,k,epsilon,...}`

### Optional comparative folder
- `turbulenceOptimizer/` (template/reference) can be used selectively once flow control is stabilized.

## 7) Success criteria

- Convergent RANS topology-optimization runs at each 1, 2, 3 m/s case.
- Stable sensitivities (bounded, smooth through continuation, no optimizer stalling).
- Improved thermal-hydraulic objective versus laminar-analogue and Darcy baselines.
- High-fidelity re-evaluation confirms geometry-level performance trends.

## 8) Risks and mitigation

- Instability in coupled k-epsilon updates
  - Mitigate with stronger under-relaxation and fixed-turbulence phases.
- Non-smooth sensitivities from aggressive projection
  - Delay sharp projection and keep continuation conservative.
- Turbulence-model mismatch in tiny channels
  - Validate against alternate RANS model in Stage B.
- Overfit to a single Reynolds state
  - Run and compare separate cases over the inlet-velocity range.

## 9) Implementation checklist

1. Finalize turbulence model and BC/initialization conventions.
2. Activate RANS dictionaries and fields.
3. Validate forward RANS on fixed topology.
4. Integrate design coupling into momentum source and `nuEff` usage.
5. Keep objective/constraint evaluation single-case per run.
6. Add continuation and projection tuning inputs.
7. Add case sweep scripts/templates for 1-3 m/s runs.
8. Lock hyperparameters after 2-3 verification loops.
9. Run post-optimization high-fidelity re-evaluations.

## 10) Reference corpus for verification and cross-checking

### 10.1 Primary source papers

- Sun, Y.; Hao, M.; Wang, Z. *Topology Optimization of Turbulent Flow Cooling Structures Based on the k- Model*, Entropy, 2023, 25, 1299, DOI: 10.3390/e25091299.
- Li, Yijun. *Design and Topology Optimization of Heat Sinks for the Cooling of Electronic Devices with Multiple Heat Sources*, Ph.D. thesis, Universit? de Nantes, 2023 (HAL: tel-04059777).
- Kaikow, H. et al. *Thermal Management for Electrification in Aircraft Engines: Optimization of Coolant System*, GT 2022 (ASME), DOI: 10.1115/GT2022-82538.

### 10.2 Full references from *Topology Optimization of Turbulent Flow Cooling Structures Based on the k- Model* (Entropy, 2023)

- [1] Li, W.H. Investigation on Convective and Conjugate Heat Transfer Characteristics of Cooling Structures in Gas Turbine Thin-Wall Blade. Ph.D. Thesis, Tsinghua University, Beijing, China, 2018.
- [2] Wu, W.; Yao, R.; Wang, J.; Su, H.; Wu, X. Leading edge impingement cooling analysis with separators of a real gas turbine blade. Appl. Therm. Eng. 2022, 208, 118275.
- [3] Hassan, H.; Shafey, N.A. 3D study of convection-radiation heat transfer of electronic chip inside enclosure cooled by heat sink. Int. J. Therm. Sci. 2021, 159, 106585.
- [4] Mukesh Kumar, P.C.; Arun Kumar, C.M. Numerical study on heat transfer performance using Al2O3/water nanouids in six circular channel heat sink for electronic chip. Mater. Today Proc. 2020, 21(1), 194201.
- [5] Hao, X.; Peng, B.; Xie, G.; Chen, Y. Efficient on-chip hotspot removal combined solution of thermoelectric cooler and mini-channel heat sink. Appl. Therm. Eng. 2016, 100, 170178.
- [6] Rao, Z.; Wang, S. A review of power battery thermal energy management. Renew. Sustain. Energy Rev. 2011, 15, 4554-4571.
- [7] Wu, K.; Zhang, Y.; Zeng, Y.Q.; Yang, J. Study on the safety performance of lithium-ion batteries. Adv. Chem. 2011, 23, 401-409.
- [8] Chen, S.; Zhang, G.; Zhu, J.; Feng, X.; Wei, X.; Ouyang, M.; Dai, H. Multi-objective optimization design and experimental investigation for a parallel liquid cooling-based Lithium-ion battery module under fast charging. Appl. Therm. Eng. 2022, 211, 118503.
- [9] Zhuang, D.; Yang, Y.; Ding, G.; Du, X.; Hu, Z. Optimization of Microchannel Heat Sink with Rhombus Fractal-like Units for Electronic Chip Cooling. Int. J. Refrig. 2020, 116, 108118.
- [10] Kose, H.A.; Yildizeli, A.; Cadirci, S. Parametric study and optimization of microchannel heat sinks with various shapes. Appl. Therm. Eng. 2022, 211, 118368.
- [11] Guo, R.; Li, L. Heat dissipation analysis and optimization of lithium-ion batteries with a novel parallel-spiral serpentine channel liquid cooling plate. Int. J. Heat Mass Transf. 2022, 189, 122706.
- [12] Tan, H.; Du, P.; Zong, K.; Meng, G.; Gao, X.; Li, Y. Investigation on the temperature distribution in the two-phase spider netted microchannel network heat sink with non-uniform heat flux. Int. J. Therm. Sci. 2021, 169, 107079.
- [13] Bendsoe, M.P.; Kikuchi, N. Generating optimal topologies in structural design using a homogenization method. Comput. Methods Appl. Mech. Eng. 1988, 71, 197-224.
- [14] Xie, Y.M.; Steven, G.P. Evolutionary Structural Optimization. Springer, London, 1997.
- [15] Bendsoe, M.P.; Sigmund, O. Topology Optimization Theory, Methods and Applications. Springer, Berlin/Heidelberg, 2004.
- [16] Wang, M.Y.; Wang, X.; Guo, D. A level set method for structural topology optimization. Comput. Methods Appl. Mech. Eng. 2003, 192, 227-246.
- [17] Guo, X.; Zhang, W.; Zhang, J.; Yuan, J. Explicit structural topology optimization based on moving morphable components (MMC) with curved skeletons. Comput. Methods Appl. Mech. Eng. 2016, 310, 711748.
- [18] Borrvall, T.; Petersson, J. Topology optimization of fluids in Stokes flow. Int. J. Numer. Methods Fluids 2003, 41, 77-107.
- [19] Gersborg-Hansen, A.; Sigmund, O.; Haber, R.B. Topology optimization of channel flow problems. Struct. Multidiscip. Optim. 2005, 30, 181-192.
- [20] Pietropaoli, M.; Ahlfeld, R.; Montomoli, F.; DErcole, M. Design for Additive Manufacturing: Internal Channel Optimization. J. Eng. Gas Turbines Power 2017, 139, 102101.
- [21] Han, X.-H.; Liu, H.-L.; Xie, G.; Sang, L.; Zhou, J. Topology optimization for spider web heat sinks for electronic cooling. Appl. Therm. Eng. 2021, 195, 117154.
- [22] Chen, F.; Wang, J.; Yang, X. Topology optimization design and numerical analysis on cold plates for lithium-ion battery thermal management. Int. J. Heat Mass Transf. 2022, 183, 122087.
- [23] Hu, D.; Zhang, Z.; Li, Q. Numerical study on flow and heat transfer characteristics of microchannels designed using topology optimization methods. Sci. China Technol. Sci. 2020, 63, 105115.
- [24] Yoon, G.H. Topology optimization for turbulent flow with Spalart-Allmaras model. Comput. Methods Appl. Mech. Eng. 2016, 303, 288-311.
- [25] Yoon, G.H. et al. Topology optimization of turbulent rotating flows using Spalart-Allmaras model. Comput. Methods Appl. Mech. Eng. 2021, 373, 113551.
- [26] Yoon, G.H. Topology optimization method with finite elements based on k-turbulence model. Comput. Methods Appl. Mech. Eng. 2020, 361, 112784.
- [27] Dilgen, C.B.; Dilgen, S.B.; Fuhrman, D.R.; Sigmund, O.; Lazarov, B.S. Topology optimization of turbulent flows. Comput. Methods Appl. Mech. Eng. 2018, 331, 363393.
- [28] Zhao, X.; Zhou, M.; Sigmund, O.; Andreasen, C.S. Poor-man’s approach to topology optimization of cooling channels based on a Darcy flow model. Int. J. Heat Mass Transf. 2018, 116, 110811.
- [29] Li, A. Topology Optimization of Fluid-Solid Conjugate Heat Transfer Structure. Ph.D. Thesis, Dalian University of Technology, 2020.
- [30] Li, B.; Xie, C.; Yin, X.; Lu, R.; Ma, Y.; Liu, H.; Hong, J. Multidisciplinary optimization of liquid cooled heat sinks with compound jet/channel structures arranged in a multipass configuration. Appl. Therm. Eng. 2021, 195, 117159.
- [31] Dilgen, S.B.; Dilgen, C.B.; Fuhrman, D.R.; Sigmund, O.; Lazarov, B.S. Density based topology optimization of turbulent flow heat transfer systems. Struct. Multidiscip. Optim. 2018, 57, 1905-1918.
- [32] Zou, A.; Chuan, R.; Qian, F.; Zhang, W.; Wang, Q.; Zhao, C. Topology optimization for a water-cooled heat sink in micro-electronics based on Pareto frontier. Appl. Therm. Eng. 2022, 207, 118128.
- [33] Lazarov, B.S.; Sigmund, O. Filters in topology optimization based on Helmholtz-type differential equations. Int. J. Numer. Methods Eng. 2011, 86, 765-781.
- [34] Li, H.; Ding, X.H.; Meng, F.Z.; Jing, D.L.; Xiong, M. Optimal design and thermal modelling for liquid-cooled heat sink based on multi-objective topology optimization: An experimental and numerical study. Int. J. Heat Mass Transf. 2019, 144, 118638.

### 10.3 Full references from *Thermal Management for Electrification in Aircraft Engines: Optimization of Coolant System*

- [1] A. Bills, S. Sripad, W. Fredericks, M. Singh, V. Viswanathan, ACS Energy Letters, 2020.
- [2] J. Hayes, K. George, P. Killeen, B. McPherson, K. Olejniczak, T. McNutt, IEEE WiPDA workshop, 2016.
- [3] R. Falck, J. Chin, S. Schnulo, J. Burt, J. Gray, AIAA/ISSMO MDO Conf., 2017.
- [4] S. Jafari; T. Nikolaidis, Applied Sciences, 2018.
- [5] B. Brelje; J. Martins, Progress in Aerospace Sciences, 2019.
- [6] M. Turkyilmazoglu, Int. J. Refrigeration, 2014.
- [7] L. Lu; X. Han; J. Li; J. Hua; M. Ouyang, Journal of Power Sources, 2013.
- [8] S. Cokun; M. Atay, Appl. Thermal Engineering, 2008.
- [9] S. Liu; Y. Huang; J. Wang, Int. J. Thermal Sciences, 2018.
- [10] I. Kaur; P. Singh, Int. J. Heat and Mass Transfer, 2021.
- [11] P. He; C. Mader; J. Martins; K. Maki, Int. J. Heat and Mass Transfer, 2019.
- [12] T. Deng; G. Zhang; Y. Ran, Int. J. Heat and Mass Transfer, 2018.
- [13] C. Mangrulkar et al., Renewable Sustainable Energy Reviews, 2019.
- [14] L. Hghj; D. Nrhave; J. Alexandersen; O. Sigmund; C. Andreasen, Int. J. Heat and Mass Transfer, 2020.
- [15] M. Pietropaoli; F. Montomoli; A. Gaymann, Struct. Multidiscip. Optim., 2019.
- [16] S. Kambampati; J. Gray; H. Kim, Int. J. Heat and Mass Transfer, 2021.
- [17] X. Mo; H. Zhi; Y. Xiao; H. Hua; L. He, Int. J. Heat and Mass Transfer, 2021.
- [18] M. Tomlin; J. Meyer, CAE conf., 2011.
- [19] L. Berrocal et al., Progress in Additive Manufacturing, 2019.
- [20] R. Murphy; C. Imediegwu; R. Hewson; M. Santer, Struct. Multidiscip. Optim., 2021.
- [21] T. Dbouk, Appl. Thermal Engineering, 2017.
- [22] K. Chen et al., International Journal of Energy Research, 2020.
- [23] X. Xu; G. Tong; R. Li, Applied Thermal Engineering, 2020.
- [24] B. Blakey-Milner et al., Materials & Design, 2021.
- [25] D. Jafari; W. Wits, Renewable Sustainable Energy Reviews, 2018.
- [26] J. Liu et al., Structural and Multidisciplinary Optimization, 2018.
- [27] C. Babu et al., Materials Today: Proceedings, 2021.
- [28] C. Sinn; G. Pesch; J. Thming and L. Kiewidt, Int. J. Heat and Mass Transfer, 2019.
- [29] A. Papukchiev; D. Grishchenko; P. Kudinov, Nuclear Engineering and Design, 2020.
- [30] S. Patankar; D. Spalding, Calculation procedure for turbulent heat transfer, 1983.
- [31] ANSYS, ANSYS Fluent User's Guide, 2019.
- [32] A. Demargne et al., Practical and reliable mesh generation for complex real-world geometries, 2014.
- [33] H. Brinkman, Flow, Turbulence and Combustion, 1949.

### 10.4 LIYIJUN thesis references for broader verification

- LIYIJUN thesis has a large reference list (173+ entries) that covers manifold design, flow maldistribution, microchannel optimization, thermal management hardware, and topology-oriented optimization. Use it for depth checks.
- Relevant entries for turbulent/thermal TO development:
  - [3] S. Zeng; B. Kanargi; P.S. Lee, Int. J. Heat Mass Transf., 2018.
  - [61] H.E. Ahmed et al., Int. J. Heat Mass Transf., 2018.
  - [86] N. Gilmore et al., uniform flow via topology optimization and flow visualization, 2021.
  - [141] S.B. Dilgen; C.B. Dilgen; D.R. Fuhrman; O. Sigmund; B.S. Lazarov, Struct. Multidiscip. Optim., 2018.
  - [152] S. Zeng; P.S. Lee, Int. J. Heat Mass Transf., 2019.
  - [155] X. Mo; H. Zhi; Y. Xiao; H. Hua; L. He, Int. J. Heat Mass Transf., 2021.
  - [156] X. Han; H. Liu; G. Xie; L. Sang; J. Zhou, Appl. Therm. Eng., 2021.
  - [167] Y.A. Manaserh et al., Int. J. Heat Mass Transf., 2022.
  - [168] F. Dugast; Y. Favennec; C. Josset; Y. Fan; L. Luo, J. Comput. Phys., 2018.
