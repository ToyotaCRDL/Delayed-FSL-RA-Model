## 1. Overview
This repository contains:

- Finite Element (FE) simulations using Abaqus
- Neural response simulations using a Leaky Integrate-and-Fire (LIF) model implemented in MATLAB
- Processed results used to generate the figures in the manuscript


## 2. Folder Structure

```text
Project_Root/
├── FE_model/
│  ├── stimulus1_probe/
│  ├── stimulus2_braille/
│  ├── stimulus3_dot_texture/
│  └── stimulus4_single_dot/
│
└── LIF_model/
    ├── stimulus1_probe/
    ├── stimulus2_braille/
    ├── stimulus3_dot_texture/
    └── stimulus4_single_dot/
```


### Folder Description

FE_model/
Abaqus finite element model files (`.inp`). Each `.inp` file includes nodal coordinates and element connectivity data for FE geometry models.

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
1\. Multi-objective Optimization
- Navigate to:<br>
&nbsp;&nbsp;&nbsp;   LIF_model/stimulus1_probe/
- Run:<br>
&nbsp;&nbsp;&nbsp;   `multiObj_Pareto.m`
- Output (variables):<br>
&nbsp;&nbsp;&nbsp;   `popSol`: Identified parameter sets (saved in param.mat)<br>
&nbsp;&nbsp;&nbsp;   `popCost`:  Pareto-optimal solutions. Exported as probe_pareto_sol.csv (stored in result1/). Used to generate Fig. 2 in the manuscript.

2\. Probe stimulation
- Navigate to:<br>
&nbsp;&nbsp;&nbsp;   LIF_model/stimulus1_probe/
- The SED outputs from the FE model are stored in `f61-1.xlsx`.
- Run:<br>
&nbsp;&nbsp;&nbsp;   `mainProbe.m` 
- Output (variables):<br>
&nbsp;&nbsp;&nbsp;   `spikeTime`: Calculated spike timing. Exported as probe_spiketrain_param_XX_ps_XX-XX.csv (stored in result2/). Here, param_XX denotes the parameter set number, and ps_XX-XX denotes the probe indentation speed. Used to generate Fig. 3.

### 5.2 Stimulus2
- Navigate to:<br>
&nbsp;&nbsp;&nbsp;   LIF_model/stimulus2_braille/
- The SED outputs from the FE model are stored in `dataSED.mat`.
- Run:<br>
&nbsp;&nbsp;&nbsp;   `mainBrailleRA.m` (for RA responses)<br>
&nbsp;&nbsp;&nbsp;   `mainBrailleSA1.m` (for SA1 responses)
- Output (variables):<br>
&nbsp;&nbsp;&nbsp;   `allSpikeTime`: Calculated spike timing. Exported as braille_sep_param_XX.csv (stored in result1/). Here, param_XX denotes the parameter set number for RA responses, whereas braille_sep_sa1.csv corresponds to SA1 responses. Used to generate Fig. 6.


### 5.3 Stimulus3
- Navigate to:<br>
&nbsp;&nbsp;&nbsp;   LIF_model/stimulus3_dot_texture/
- The SED outputs from the FE model are stored in `dataSED.mat`.
- Run:<br>
&nbsp;&nbsp;&nbsp;   `mainDotTexRA.m` (for RA responses)<br>
&nbsp;&nbsp;&nbsp;   `mainDotTexSA1.m` (for SA1 responses)
- Output (variables):<br>
&nbsp;&nbsp;&nbsp;   `allSpikeTime`: Calculated spike timing. Exported as dot_texture_sep_param_XX_ds_XX-XX.csv (stored in result1/). Here, param_XX denotes the parameter set number for RA responses, whereas dot_texture_sep_sa1_ds_XX-XX.csv corresponds to SA1 responses. ds_XX-XX denotes the dot spacing. Used to generate Fig. 7.<br>
&nbsp;&nbsp;&nbsp;   `setFiringRateRMSE`: Calculated RMSE for firing rate. Exported as dot_texture_firerate_rmse.csv (stored in result2/). The results in XX_rmse.csv are sorted in ascending order of the parameter set number, whereas the results in XX_rmse2.csv are sorted in ascending order of &Delta;FST. Fig. 8 was generated using XX_rmse2.csv.<br>
&nbsp;&nbsp;&nbsp;   `setTempoVariRMSE`: Calculated RMSE for temporal variation. Exported as dot_texture_tempvar_rmse.csv (stored in result2/). Fig. 8 was generated using XX_rmse2.csv.<br>
&nbsp;&nbsp;&nbsp;   `vec_allFiringRateHistogramArray`: Calculated firing rate. Fig. 9(A) was generated using dot_texture_firerate_diff2.csv. The results in XX_diff.csv are sorted in ascending order of the parameter set number, whereas the results in XX_diff2.csv are sorted in ascending order of &Delta;FST.<br>
&nbsp;&nbsp;&nbsp;   `vec_tempoVariArray`: Calculated temporal variation. Exported as dot_texture_tempvar_sim.csv (stored in result2/). Fig. 9(B) was generated using dot_texture_tempvar_diff2.csv.
   

### 5.4 Stimulus4
- Navigate to:<br>
&nbsp;&nbsp;&nbsp;   LIF_model/stimulus4_single_dot/
- The SED outputs from the FE model are stored in `dataSED.mat`.
- Run:<br>
&nbsp;&nbsp;&nbsp;   `mainSingleDotRA.m` (for RA responses)<br>
&nbsp;&nbsp;&nbsp;   `mainSingleDotSA1.m` (for SA1 responses)
- Output (variables):<br>
&nbsp;&nbsp;&nbsp;   `allSpikeTime`: Calculated spike timing. Exported as single_dot_sep_param_XX_dw_XX-XX.csv (stored in result1/). Here, param_XX denotes the parameter set number for RA responses, whereas single_dot_sep_sa1_dw_XX-XX.csv corresponds to SA1 responses. dw_XX-XX denotes the dot width. Used to generate Fig. 10.<br>
&nbsp;&nbsp;&nbsp;   `setRespWidth`: Calculated response width. Exported as single_dot_respwidth.csv (stored in result2/). The results in XX_respwidth.csv are sorted in ascending order of the parameter set number, whereas the results in XX_respwidth2.csv are sorted in ascending order of &Delta;FST. In addition, the response widths in XX_respwidth2.csv represent the differences between the simulation results and the experimental data. Fig. 11 was generated using XX_respwidth2.csv.


### 5.5 Effect of unresponsive period (Stimulus3)
- Navigate to:<br>
&nbsp;&nbsp;&nbsp;   LIF_model/stimulus3_dot_texture_delay/
- The SED outputs from the FE model are stored in `dataSED.mat`.
- Run:<br>
&nbsp;&nbsp;&nbsp;   `mainDotTexRA.m` (for RA responses)
- Output (variables):<br>
&nbsp;&nbsp;&nbsp;   `allSpikeTime`: Calculated spike timing. Exported as dot_texture_sep_param_13_ds_XX-XX_delay_XX.csv (stored in result1/). Used to generate Fig. XX.<br>
&nbsp;&nbsp;&nbsp;   `setFiringRateRMSE`: Calculated RMSE for firing rate. Exported as dot_texture_firerate_rmse.csv (stored in result2/). Used to generate Fig. XX.<br>
&nbsp;&nbsp;&nbsp;   `setTempoVariRMSE`: Calculated RMSE for temporal variation. Exported as dot_texture_tempvar_rmse.csv (stored in result2/). Used to generate Fig. XX.


### 5.6 Effect of unresponsive period (Stimulus4)
- Navigate to:<br>
&nbsp;&nbsp;&nbsp;   LIF_model/stimulus4_single_dot_delay/
- The SED outputs from the FE model are stored in `dataSED.mat`.
- Run:<br>
&nbsp;&nbsp;&nbsp;   `mainSingleDotRA.m` (for RA responses)
- Output (variables):<br>
&nbsp;&nbsp;&nbsp;   `allSpikeTime`: Calculated spike timing. Exported as single_dot_sep_param_13_dw_XX-XX.csv (stored in result1/). Used to generate Fig. XX.<br>
&nbsp;&nbsp;&nbsp;   `setRespWidth`: Calculated response width. Exported as single_dot_respwidth.csv (stored in result2/). The response widths in XX_respwidth2.csv represent the differences between the simulation results and the experimental data. Used to generate Fig. XX..


### 5.7 Effect of neuronal noise (Stimulus3)
- Navigate to:<br>
&nbsp;&nbsp;&nbsp;   LIF_model/stimulus3_dot_texture_noise/
- The SED outputs from the FE model are stored in `dataSED.mat`.
- Run:<br>
&nbsp;&nbsp;&nbsp;   `mainDotTexRA.m` (for RA responses)
- Output (variables):<br>
&nbsp;&nbsp;&nbsp;   `allSpikeTime`: Calculated spike timing. Exported as dot_texture_sep_param_XX_ds_XX-XX_noise_XX.csv (stored in result1/). Used to generate Fig. XX.<br>
&nbsp;&nbsp;&nbsp;   `setFiringRateRMSE`: Calculated RMSE for firing rate. Exported as dot_texture_firerate_rmse.csv (stored in result2/). The results in XX_rmse.csv are sorted in ascending order of the parameter set number, whereas the results in XX_rmse2.csv are sorted in ascending order of &Delta;FST. Fig. XX was generated using XX_rmse2.csv.<br>
&nbsp;&nbsp;&nbsp;   `setTempoVariRMSE`: Calculated RMSE for temporal variation. Exported as dot_texture_tempvar_rmse.csv (stored in result2/). Fig. XX was generated using XX_rmse2.csv.


### 5.8 Effect of neuronal noise (Stimulus4)
- Navigate to:<br>
&nbsp;&nbsp;&nbsp;   LIF_model/stimulus4_single_dot_noise/
- The SED outputs from the FE model are stored in `dataSED.mat`.
- Run:<br>
&nbsp;&nbsp;&nbsp;   `mainSingleDotRA.m` (for RA responses)
- Output (variables):<br>
&nbsp;&nbsp;&nbsp;   `allSpikeTime`: Calculated spike timing. Exported as single_dot_sep_param_XX_dw_XX-XX_noise_XX.csv (stored in result1/). Used to generate Fig. XX.<br>
&nbsp;&nbsp;&nbsp;   `setRespWidth`: Calculated response width. Exported as single_dot_respwidth.csv (stored in result2/). The results in XX_respwidth.csv are sorted in ascending order of the parameter set number, whereas the results in XX_respwidth2.csv are sorted in ascending order of &Delta;FST. In addition, the response widths in XX_respwidth2.csv represent the differences between the simulation results and the experimental data. Fig. XX was generated using XX_respwidth2.csv.


### 5.9 Effect of response width definition (Stimulus4)
- Navigate to:<br>
&nbsp;&nbsp;&nbsp;   LIF_model/stimulus4_single_dot_modi_RW/
- The SED outputs from the FE model are stored in `dataSED.mat`.
- Run:<br>
&nbsp;&nbsp;&nbsp;   `mainSingleDotRA.m` (for RA responses)
- Output (variables):<br>
&nbsp;&nbsp;&nbsp;   `allSpikeTime`: Calculated spike timing. Exported as single_dot_sep_param_XX_dw_XX-XX_noise_XX.csv (stored in result1/). Used to generate Fig. XX.<br>
&nbsp;&nbsp;&nbsp;   `spikeWidth`: Calculated response width. Exported as single_dot_respwidth.csv (stored in result2/). The results in XX_respwidth.csv are sorted in ascending order of the parameter set number, whereas the results in XX_respwidth2.csv are sorted in ascending order of &Delta;FST.

