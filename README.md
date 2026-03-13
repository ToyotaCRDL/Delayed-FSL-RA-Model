## 1. Overview
This repository contains:

- Finite Element (FE) simulations using Abaqus
- Neural response simulations using a Leaky Integrate-and-Fire (LIF) model implemented in MATLAB
- Processed results used to generate the figures in the manuscript


## 2. Folder Structure

Project_Root/
|
|----- FE_model/
| |----- stimulus1_probe/
| |----- stimulus2_braille/
| |----- stimulus3_dot_texture/
| |----- stimulus4_single_dot/
|
|----- LIF_model/
| |----- stimulus1_probe/
| |----- stimulus2_braille/
| |----- stimulus3_dot_texture/
| |----- stimulus4_single_dot/

### Folder Description

FE_model/
Abaqus finite element model files (.inp).

LIF_model/
MATLAB scripts (.m) implementing the neural response (LIF) model and corresponding simulation output for each simulation condition.

stimulusXX_XX/
Each folder corresponds to one simulation condition described in the manuscript.


## 3. Software Requirements
Finite Element Simulation
Abaqus Version: Abaqus 2016 and Abaqus 2022

Neural Simulation
MATLAB Version: R2024b
Required toolboxes: Global Optimization Toolbox


## 4. Overall Workflow
The simulation pipeline consists of the following steps:
1. Run Abaqus FE simulation (.inp)
2. Extract strain enelgy density data from FE output (extracted data: abaqusXX.rpt and fXX-X.xlsx, stored in LIF_model/stimulusXX_XX/fingerXX/)
3. Run the MATLAB-based LIF simulation using the extracted FE data as input (.m)
4. Neural response outputs


## 5. Detailed Execution Procedure

### 5.1 Stimulus1
1. Multi-objective Optimization
- Navigate to:
   LIF_model/stimulus1_probe/
- Run:
   `multiObj_Pareto.m`
- Output (variables):
   `popSol`: Identified parameter sets (saved in param.mat)
   `popCost`:  Pareto-optimal solutions. Exported as probe_pareto_sol.csv (stored in result1/). Used to generate Fig. 2 in the manuscript.

2. Probe stimulation
- Navigate to:
   LIF_model/stimulus1_probe/
- Run:
   `mainProbe.m` 
- Output (variables):
   `spikeTime`: Calculated spike timing. Exported as probe_spiketrain_param_XX_ps_XX-XX.csv (stored in result2/). Here, param_XX denotes the parameter set number, and ps_XX-XX denotes the probe indentation speed. Used to generate Fig. 3.

### 5.2 Stimulus2
- Navigate to:
   LIF_model/stimulus2_braille/
- Run:
   `mainBrailleRA.m` (for RA responses)
   `mainBrailleSA1.m` (for SA1 responses)
- Output (variables):
   `allSpikeTime`: Calculated spike timing. Exported as braille_sep_param_XX.csv (stored in result1/). Here, param_XX denotes the parameter set number for RA responses, whereas braille_sep_sa1.csv corresponds to SA1 responses. Used to generate Fig. 6.


### 5.3 Stimulus3
- Navigate to:
   LIF_model/stimulus3_braille/
- Run:
   `mainDotTexRA.m` (for RA responses)
   `mainDotTexSA1.m` (for SA1 responses)
- Output (variables):
   `allSpikeTime`: Calculated spike timing. Exported as dot_texture_sep_param_XX_ds_XX-XX.csv (stored in result1/). Here, param_XX denotes the parameter set number for RA responses, whereas dot_texture_sep_sa1_ds_XX-XX.csv corresponds to SA1 responses. ds_XX-XX denotes the dot spacing. Used to generate Fig. 7.
   `vec_allFiringRateHistogramArray`: Calculated firing rate. Exported as dot_texture_firerate_sim.csv (stored in result2/). Used to generate Fig. 8(A-1). Fig. 9(A) was generated using dot_texture_firerate_diff2.csv. The results in XX_diff.csv are sorted in ascending order of the parameter set number, whereas the results in XX_diff2.csv are sorted in ascending order of ΔFST.
   `vec_tempoVariArray`: Calculated temporal variation. Exported as dot_texture_tempvar_sim.csv (stored in result2/). Used to generate Fig. 8(A-2). Fig. 9(B) was generated using dot_texture_tempvar_diff2.csv.
   `setFiringRateRMSE`: Calculated RMSE for firing rate. Exported as dot_texture_firerate_rmse.csv (stored in result2/). The results in XX_rmse.csv are sorted in ascending order of the parameter set number, whereas the results in XX_rmse2.csv are sorted in ascending order of ΔFST. Fig. 8(B) was generated using XX_rmse2.csv.
   `setTempoVariRMSE`: Calculated RMSE for temporal variation. Exported as dot_texture_tempvar_rmse.csv (stored in result2/). Fig. 8(B) was generated using XX_rmse2.csv.


### 5.4 Stimulus4
- Navigate to:
   LIF_model/stimulus4_single_dot/
- Run:
   `mainSingleDotRA.m` (for RA responses)
   `mainSingleDotSA1.m` (for SA1 responses)
- Output (variables):
   `allSpikeTime`: Calculated spike timing. Exported as single_dot_sep_param_XX_dw_XX-XX.csv (stored in result1/). Here, param_XX denotes the parameter set number for RA responses, whereas single_dot_sep_sa1_dw_XX-XX.csv corresponds to SA1 responses. dw_XX-XX denotes the dot width. Used to generate Fig. 10.
   `setRespWidth`: Calculated response width. Exported as single_dot_respwidth.csv (stored in result2/). The results in XX_respwidth.csv are sorted in ascending order of the parameter set number, whereas the results in XX_respwidth2.csv are sorted in ascending order of ΔFST. In addition, the response widths in XX_respwidth2.csv represent the differences between the simulation results and the experimental data. Fig. 11 was generated using XX_respwidth2.csv.
