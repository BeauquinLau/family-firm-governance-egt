# When to Meet the Equifinality? Opening the Black Box of Family Firms’ Governance Heterogeneity through the Perspective of Evolutionary Game Theory

This repository contains the source code for the numerical experiments presented in the associated manuscript. For a better reproduction experience, a detailed illustration is provided below.

## Repository Structure

The repository is organized into two main subdirectories within the `Simulation` folder: one for the MATLAB source code and one for the generated figures.

```text
Simulation/
│
├── Figures_eps/
│   ├── F1_1.eps, F1_2.eps
│   ├── F2_1.eps, F2_2.eps
│   ├── F3.eps
│   ├── F4_1.eps, F4_2.eps, F4_3.eps
│   ├── F5_1.eps, F5_2.eps
│   ├── F6_1.eps, F6_2.eps
│   └── F7_1.eps, F7_2.eps
│
└── MATLAB_code/
    ├── S1_1.m
    ├── S1_2.m
    ├── S2_1.m
    ├── S2_2.m
    ├── Dynamic_trajectories.m
    ├── A.m
    ├── B.m
    ├── C.m
    ├── A_p.m
    ├── A_q.m
    ├── B_p.m
    ├── B_q.m
    ├── C_p.m
    └── C_q.m
```

## System Requirements

*   **MATLAB:** The code has been tested on MATLAB R2025a and is expected to be compatible with newer versions.
*   **Toolboxes:** No special toolboxes are required for execution. The code relies only on the base MATLAB environment.

## Instructions for Reproduction

This repository provides executable scripts to generate all the figures presented in the manuscript.

**General Instruction:** To run any script, please first set MATLAB's current working directory to the `Simulation/MATLAB_code/` folder. All generated figures will be saved as `.eps` files in the `../Figures_eps/` directory relative to the script's location.

### 1. System Dynamics and Stability Analysis (Figures 1 & 2)

These scripts illustrate the system's phase portrait and the stability of equilibrium points under two different sets of model parameters.

*   **`S1_1.m` & `S1_2.m`**: Generate the phase portrait and the stable/unstable manifolds for the baseline parameter set.
    *   **Outputs**: `F1_1.eps`, `F1_2.eps`
*   **`S2_1.m` & `S2_2.m`**: Generate the phase portrait and the stable/unstable manifolds for the second parameter set.
    *   **Outputs**: `F2_1.eps`, `F2_2.eps`

📌 **To run (example for Figure 1):**
```matlab
% In MATLAB, with the current directory set to 'Simulation/MATLAB_code/'
run('S1_1.m');
run('S1_2.m');
```

### 2. Sensitivity Analysis of Parameter `n` (Figure 3)

This script shows how the stable and unstable manifolds around the saddle point shift as the key parameter `n` varies.

*   **`Dynamic_trajectories.m`**: Loops through a range of `n` values and plots their corresponding manifolds.
    *   **Output**: `F3.eps`

📌 **To run:**
```matlab
% In MATLAB, with the current directory set to 'Simulation/MATLAB_code/'
run('Dynamic_trajectories.m');
```

### 3. Evolutionary Trajectories from Key Initial Points (Figures 4-7)

These scripts demonstrate the system's evolutionary paths starting from three key initial points (A, B, and C), showing both the phase portraits and the time evolution of the state variables `p` and `q`.

#### 3.1. Phase Portraits (Figure 4)

*   **`A.m`**: Generates the phase portrait starting from point A.
    *   **Output**: `F4_1.eps`
*   **`B.m`**: Generates the phase portrait starting from point B.
    *   **Output**: `F4_2.eps`
*   **`C.m`**: Generates the phase portrait starting from point C.
    *   **Output**: `F4_3.eps`

#### 3.2. Time Evolution Plots (Figures 5-7)

*   **`A_p.m`** and **`A_q.m`**: Plot the time evolution of `p` and `q` starting from point A.
    *   **Outputs**: `F5_1.eps`, `F5_2.eps`
*   **`B_p.m`** and **`B_q.m`**: Plot the time evolution of `p` and `q` starting from point B.
    *   **Outputs**: `F6_1.eps`, `F6_2.eps`
*   **`C_p.m`** and **`C_q.m`**: Plot the time evolution of `p` and `q` starting from point C.
    *   **Outputs**: `F7_1.eps`, `F7_2.eps`

📌 **To run (example for initial point A):**
```matlab
% In MATLAB, with the current directory set to 'Simulation/MATLAB_code/'
run('A.m');
run('A_p.m');
run('A_q.m');
```
