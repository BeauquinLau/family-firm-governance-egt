# The Gravity of Control: A Dynamic Theory of Socioemotional Wealth and Governance Traps in Family Firms

This repository contains the MATLAB source code for the numerical experiments presented in the associated manuscript. This document provides detailed instructions to ensure the clear and accurate reproduction of all figures.

## Repository Structure

The repository is organized into two main subdirectories: `Figures_PDF/` contains all pre-generated figures from the manuscript for easy reference, and `MATLAB_Code/` contains the source code required to reproduce them. The code is structured to mirror the analyses in the paper.

```
Family_Firm_Governance_Dynamics/
│
├── Figures_PDF/
│   ├── S1_phase_portrait.pdf
│   ├── S2_phase_portrait.pdf
│   ├── Sensitivity_delta_SEW_thm1.pdf
│   ├── Sensitivity_delta_SEW_thm2.pdf
│   ├── Sensitivity_f_PM_thm1.pdf
│   ├── Sensitivity_f_PM_thm2.pdf
│   ├── Sensitivity_n_thm1.pdf
│   ├── Sensitivity_n_thm2.pdf
│   ├── Trajectory_delta_SEW_thm1_S1_0.1_0.1.pdf
│   ├── Trajectory_delta_SEW_thm1_S2_0.6_0.5.pdf
│   ├── Trajectory_delta_SEW_thm1_T1_0.2_0.4.pdf
│   ├── Trajectory_delta_SEW_thm2_S3_0.2_0.9.pdf
│   ├── Trajectory_delta_SEW_thm2_S4_0.6_0.4.pdf
│   ├── Trajectory_delta_SEW_thm2_T2_0.5_0.7.pdf
│   ├── Trajectory_f_PM_thm1_S5_0.05_0.05.pdf
│   ├── Trajectory_f_PM_thm1_S6_0.6_0.7.pdf
│   ├── Trajectory_f_PM_thm1_T3_0.2_0.3.pdf
│   ├── Trajectory_f_PM_thm2_S7_0.2_0.8.pdf
│   ├── Trajectory_f_PM_thm2_S8_0.7_0.3.pdf
│   ├── Trajectory_f_PM_thm2_T4_0.5_0.6.pdf
│   ├── Trajectory_n_thm1_R1_0.7_0.11.pdf
│   ├── Trajectory_n_thm1_S9_0.3_0.1.pdf
│   ├── Trajectory_n_thm1_S10_0.5_0.4.pdf
│   ├── Trajectory_n_thm1_T5_0.2_0.35.pdf
│   ├── Trajectory_n_thm1_T6_0.8_0.06.pdf
│   ├── Trajectory_n_thm2_R2_0.9237_0.8840.pdf
│   ├── Trajectory_n_thm2_S11_0.4_0.7.pdf
│   ├── Trajectory_n_thm2_S12_0.6_0.4.pdf
│   ├── Trajectory_n_thm2_T7_0.3_0.4.pdf
│   └── Trajectory_n_thm2_T8_0.975740_0.942642.pdf
│
└── MATLAB_Code/
    │
    ├── Governance_Landscapes/
    │   ├── S1_phase_portrait.m
    │   └── S2_phase_portrait.m
    │
    └── Sensitivity_Analysis/
        │
        ├── Analysis_delta_SEW/
        │   ├── Sensitivity_delta_SEW_thm1.m
        │   ├── Sensitivity_delta_SEW_thm2.m
        │   └── Trajectories/
        │       ├── Trajectory_delta_SEW_thm1_S1_0.1_0.1.m
        │       ├── Trajectory_delta_SEW_thm1_S2_0.6_0.5.m
        │       ├── Trajectory_delta_SEW_thm1_T1_0.2_0.4.m
        │       ├── Trajectory_delta_SEW_thm2_S3_0.2_0.9.m
        │       ├── Trajectory_delta_SEW_thm2_S4_0.6_0.4.m
        │       └── Trajectory_delta_SEW_thm2_T2_0.5_0.7.m
        │
        ├── Analysis_f_PM/
        │   ├── Sensitivity_f_PM_thm1.m
        │   ├── Sensitivity_f_PM_thm2.m
        │   └── Trajectories/
        │       ├── Trajectory_f_PM_thm1_S5_0.05_0.05.m
        │       ├── Trajectory_f_PM_thm1_S6_0.6_0.7.m
        │       ├── Trajectory_f_PM_thm1_T3_0.2_0.3.m
        │       ├── Trajectory_f_PM_thm2_S7_0.2_0.8.m
        │       ├── Trajectory_f_PM_thm2_S8_0.7_0.3.m
        │       └── Trajectory_f_PM_thm2_T4_0.5_0.6.m
        │
        └── Analysis_n/
            ├── Sensitivity_n_thm1.m
            ├── Sensitivity_n_thm2.m
            └── Trajectories/
                ├── Trajectory_n_thm1_R1_0.7_0.11.m
                ├── Trajectory_n_thm1_S9_0.3_0.1.m
                ├── Trajectory_n_thm1_S10_0.5_0.4.m
                ├── Trajectory_n_thm1_T5_0.2_0.35.m
                ├── Trajectory_n_thm1_T6_0.8_0.06.m
                ├── Trajectory_n_thm2_R2_0.9237_0.8840.m
                ├── Trajectory_n_thm2_S11_0.4_0.7.m
                ├── Trajectory_n_thm2_S12_0.6_0.4.m
                ├── Trajectory_n_thm2_T7_0.3_0.4.m
                └── Trajectory_n_thm2_T8_0.975740_0.942642.m
```

## System Requirements

-   **MATLAB**: The code has been tested on MATLAB R2025b and is expected to be compatible with recent versions.
-   **Toolboxes**: No special toolboxes are required. The code relies only on the base MATLAB environment.

## Instructions for Reproduction

This repository provides executable scripts to generate all figures from the manuscript's numerical simulation and sensitivity analysis sections.

> **General Instruction:** To run any script, please first set MATLAB's current working directory to the script's location. All generated figures will be saved as `.pdf` files in the same directory.

### 1. Governance Landscapes: Visualizing the Main Theorems

These scripts generate the foundational phase portraits that illustrate the two forms of bistability.

-   **Location**: `MATLAB_Code/Governance_Landscapes/`
-   **Script**: `S1_phase_portrait.m`
    -   **Output**: `S1_phase_portrait.pdf` (Corresponds to Figure `\ref{fig:simulations}(a)`).
-   **Script**: `S2_phase_portrait.m`
    -   **Output**: `S2_phase_portrait.pdf` (Corresponds to Figure `\ref{fig:simulations}(b)`).

📌 **To run:**

```matlab
% In MATLAB, navigate to the 'MATLAB_Code/Governance_Landscapes/' directory
run('S1_phase_portrait.m');
run('S2_phase_portrait.m');
```

### 2. The Architecture of Path Dependence: Sensitivity Analysis

These scripts demonstrate how the system's separatrix and evolutionary paths shift as key parameters are varied. Each sub-analysis is organized into its own folder.

#### 2.1. The Legacy Anchor: SEW Premium as a Bifurcation Driver

-   **Location**: `MATLAB_Code/Sensitivity_Analysis/Analysis_delta_SEW/`
-   **Global View - Separatrix Locus (Figure \ref{fig:sensitivity_SEW_global})**:
    -   **Scripts**: `Sensitivity_delta_SEW_thm1.m` and `Sensitivity_delta_SEW_thm2.m`.
    -   **Outputs**: `Sensitivity_delta_SEW_thm1.pdf` (Fig. \ref{fig:sensitivity_SEW_global}(a)) and `Sensitivity_delta_SEW_thm2.pdf` (Fig. \ref{fig:sensitivity_SEW_global}(b)).
-   **Micro-Cases - Trajectory Analysis (Figure \ref{fig:trajectory_cases_SEW})**:
    -   **Location**: `Trajectories/` subfolder.
    -   **Scripts**: `Trajectory_delta_SEW_... .m` (6 files).
    -   **Outputs**: `Trajectory_delta_SEW_... .pdf` files corresponding to panels (a) through (f) of Figure \ref{fig:trajectory_cases_SEW}.

#### 2.2. Trust as a Catalyst: Stewardship Utility as a Path-Shaping Force

-   **Location**: `MATLAB_Code/Sensitivity_Analysis/Analysis_f_PM/`
-   **Global View - Separatrix Locus (Figure \ref{fig:sensitivity_fPM_global})**:
    -   **Scripts**: `Sensitivity_f_PM_thm1.m` and `Sensitivity_f_PM_thm2.m`.
    -   **Outputs**: `Sensitivity_f_PM_thm1.pdf` (Fig. \ref{fig:sensitivity_fPM_global}(a)) and `Sensitivity_f_PM_thm2.pdf` (Fig. \ref{fig:sensitivity_fPM_global}(b)).
-   **Micro-Cases - Trajectory Analysis (Figure \ref{fig:trajectory_cases_fPM})**:
    -   **Location**: `Trajectories/` subfolder.
    -   **Scripts**: `Trajectory_f_PM_... .m` (6 files).
    -   **Outputs**: `Trajectory_f_PM_... .pdf` files corresponding to panels (a) through (f) of Figure \ref{fig:trajectory_cases_fPM}.

#### 2.3. The Disciplining Hand: External Monitoring as a Landscape Shifter

-   **Location**: `MATLAB_Code/Sensitivity_Analysis/Analysis_n/`
-   **Global View - Separatrix Locus (Figure \ref{fig:sensitivity_n_global})**:
    -   **Scripts**: `Sensitivity_n_thm1.m` and `Sensitivity_n_thm2.m`.
    -   **Outputs**: `Sensitivity_n_thm1.pdf` (Fig. \ref{fig:sensitivity_n_global}(a)) and `Sensitivity_n_thm2.pdf` (Fig. \ref{fig:sensitivity_n_global}(b)).
-   **Micro-Cases - Trajectory Analysis (Figures \ref{fig:trajectory_cases_n_thm1} & \ref{fig:trajectory_cases_n_thm2})**:
    -   **Location**: `Trajectories/` subfolder.
    -   **Scripts**: `Trajectory_n_... .m` (10 files).
    -   **Outputs**: `Trajectory_n_... .pdf` files corresponding to the panels in Figures \ref{fig:trajectory_cases_n_thm1} and \ref{fig:trajectory_cases_n_thm2}.

📌 **To run (example for the full analysis of parameter $\Delta_{\mathrm{SEW}}$):**

```matlab
% In MATLAB, navigate to 'MATLAB_Code/Sensitivity_Analysis/Analysis_delta_SEW/'
run('Sensitivity_delta_SEW_thm1.m');
run('Sensitivity_delta_SEW_thm2.m');

% Then navigate to the subfolder to run the micro-cases
cd('Trajectories');
run('Trajectory_delta_SEW_thm1_S1_0.1_0.1.m');
run('Trajectory_delta_SEW_thm1_S2_0.6_0.5.m');
% ... and so on for all trajectory scripts in this folder.
```
