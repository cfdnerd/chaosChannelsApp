# PlanTurbulenceOpt: Turbulent Topology Optimization Framework Implementation Plan

## 1. Executive Summary & Goal
The target is to develop a topology optimization (TO) framework capable of optimizing internal microchannel heat sinks under turbulent flow regimes (1-3 m/s inlet velocities) for electronic cooling purposes. 

After reviewing the existing `laminarOptimizer` framework, it serves as a robust foundation utilizing a continuous adjoint formulation combined with MMA (Method of Moving Asymptotes). This plan outlines how we augment `laminarOptimizer` to handle turbulent flow physics based on insights derived from recent academic publications.

## 2. Literature Review & Technical Know-how Synthesis
Based on the repository's journal articles covering both laminar and turbulent TO strategies (Sun et al. 2023, Li Yijun's thesis, and Kaikow et al./ASME GT2022 on aircraft engine thermal management):

### LAMINAR BASELINE:
- The base code (`laminarOptimizer`) successfully uses a continuous adjoint solver solving the Navier-Stokes and convection-diffusion heat equations.
- Porosity is modeled using Brinkman penalization ($$\alpha$$ term) driven by the design variable $$x$$. 
- Gradient-based optimization employs Helmholtz filtering and a projection scheme to push physical properties towards binary 0/1 states.

### TURBULENT EXTENSION STRATEGY:
- **Governing Equations**: The RANS `k-\epsilon` or `k-\omega` SST models are predominantly used in literature for this regime. `k-\epsilon` operates reasonably well under moderate internal flow turbulence and couples straightforwardly with the Brinkman penalty term.
- **Primal Coupling**: As the design variable defines solid regions, the RANS turbulence transport equations must be solved with an appropriate penalty in solid zones to suppress unphysical turbulence generation and ensure numerical stability. The velocity penalties are already in place, but turbulent viscosity and dissipation must be similarly quenched/limited in solid components.
- **Adjoint Formulation (Frozen Turbulence)**: Deriving exact adjoints for highly non-linear RANS equations is notoriously unstable and computationally expensive. Sun et al. and others validate that a **frozen-turbulence assumption** is perfectly acceptable for the adjoint equations. In this approach, the variations of turbulent viscosity ($$\nu_t$$) with respect to the design variable are ignored during sensitivity calculations. The adjoint transport simply utilizes the *effective viscosity* ($$\nu_{eff} = \nu + \nu_t$$) and *effective thermal conductivity* ($$k_{eff} + k_t$$).

## 3. Starting Point & Foundation
The `laminarOptimizer` OpenFOAM solver provides all the mechanics needed for the gradient calculation cycle: primal flow/heat update, sensitivity generation, filtering, objective penalization, and MMA updating.
**Decision:** Yes, `laminarOptimizer` acts as an excellent starting block. We copied this to **`MTOTurbulence`** to preserve the laminar working state and begin modifying for turbulence.

## 4. Required Modifications & Adjustments made in `MTOTurbulence`

### 4.1 Structural & Naming Changes:
- **Directory duplicated**: `laminarOptimizer` -> `MTOTurbulence`.
- **Executable renamed**: Refactored `MTO_ThermalFluid.C` to `MTO_Turbulence.C`.
- **Make/files updated**: Set target compilation output `EXE` to `$(FOAM_USER_APPBIN)/MTO_Turbulence`.

### 4.2 Primal Flow Fixes:
- The base `Primal_U.H` actually contains `+ turbulence->divDevReff(U)`. This natively works for both laminar and turbulent models in OpenFOAM depending on the simulation configuration (`constant/turbulenceProperties`). Thus, the primal code is largely already set up for RANS. 
- OpenFOAM handles solving the `k` and `epsilon` scalar transport equations implicitly during the SIMPLE loop if `turbulence->correct()` is called. We just need to make sure the standard RANS objects are handled effectively.

### 4.3 Adjoint Flow Fixes (Frozen Turbulence Approach):
To implement the frozen turbulence adjoint assumption evaluated from the literature, the following crucial modifications have been integrated directly into the `MTOTurbulence` source code:
- **Adjoint Velocity (`AdjointFlow_Ua.H`)**: 
  - Substituted the laminar diffusion term `- fvm::laplacian(nu,Ua)` with the turbulent effective diffusion term `- fvm::laplacian(turbulence->nuEff(),Ua)`.
- **Adjoint Heat (`AdjointHeat_Ub.H`)**:
  - Substituted `- fvm::laplacian(nu,Ub)` with `- fvm::laplacian(turbulence->nuEff(),Ub)`.
  - *Future Work Note*: Requires the injection of the turbulent Prandtl number ($$Pr_t$$) into the effective thermal conductivity to fully match the adjoint thermal transport in high-Reynolds conditions.
- By incorporating `turbulence->nuEff()`, the adjoint momentum transport now properly accounts for turbulent momentum mixing dynamically generated from the primary flow solutions without trying to differentiate the full `k-epsilon` variables.

### 4.4 Turbulence Variable Penalization & Boundary Mechanics (Crucial Additions):
* **Supressing Turbulence in Solids (`k` & `epsilon` Penalization)**: As the design material dictates solid ($$x=0$$) or fluid ($$x=1$$), it is critical to heavily suppress turbulence variables ($$k$$ and $$\epsilon$$) inside the "solid" domains using the same Brinkman penalizations ($$\alpha$$) used for velocity. Unphysical turbulence generation inside these regions will otherwise pollute the field and break the frozen-turbulence assumption stability. This requires coupling the porosity/permeability directly into the OpenFOAM turbulence model source equations or handling it via `fvOptions`.
* **Standard Wall Functions**: For internal heat sink microchannel topology optimization, attempting to fully resolve the turbulent boundary layer ($$y^+ \approx 1$$) everywhere forces the required mesh resolution beyond what is feasible for 3D TO. Therefore, the setup **must** strictly enforce High-Reynolds standard wall functions (`kqRWallFunction`, `epsilonWallFunction`, `nutkWallFunction`) locally across the walls and optimize focusing on macro-design shifts rather than micro-boundary physics.

## 5. Next Steps for Simulation Configuration
Before running a target test case, the internal case files (e.g., in a future `/app` folder to be copied over) need adjustments:
1. **fvSchemes** & **fvSolution**: Must be appended to support robust bounded upwind discretization and linear solvers for `k` and `epsilon`.
2. **0/ Fields**: Initial boundaries for `k`, `epsilon`, and `nut` must be created with realistically scaled Dirichlet values at the inlet for 1-3 m/s channel flow.
3. **turbulenceProperties**: Changed from `laminar` to `RAS` (specifically `kEpsilon`).

## 6. Stability & Robustness Plan
To prevent divergence when introducing the high non-linearity of turbulence generation:
- *Continuation Strategy*: Optimize on a purely laminar or extremely low-velocity case first, gradually increasing inlet boundary velocity towards fully turbulent operational modes.
- *Relaxation Factors*: Substantially lower field under-relaxation factors initially for `p` (0.2), `U` (0.5), `k` (0.4), and `epsilon` (0.4) during the primary SIMPLE momentum predictor-corrector cycles.
- *Brinkman Ceiling*: Constrain `alphaMax` conservatively to prevent the matrix conditioning from exploding the RANS turbulent term source bounds.
