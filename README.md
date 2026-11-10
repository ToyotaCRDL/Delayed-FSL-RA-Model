## 1. Overview
This repository contains:

- Finite Element (FE) simulations using Abaqus
- Neural response simulations using a Leaky Integrate-and-Fire (LIF) model implemented in MATLAB and Python
- Processed results used to generate the figures in the manuscript


### Versions
- Tag `v1.0`: code and results of the original submission.
- Branch `main`: revised version. Changes from v1.0:
  1. SED threshold: the afferent current is I = β + K<sub>d</sub>|dU/dt| when |dU/dt| > θ and 0 otherwise, with θ = 0.13 J m<sup>-3</sup> ms<sup>-1</sup> for all stimuli. In v1.0 the threshold was applied to the SED increment per time step (1.1 J m<sup>-3</sup>), which corresponded to 2.2 J m<sup>-3</sup> ms<sup>-1</sup> for the probe and 0.44 J m<sup>-3</sup> ms<sup>-1</sup> for scanning.
  2. Time step: 0.0125 ms for all stimuli; the SED from the FE model (output every 2.5 ms) is spline-resampled in time. In v1.0 the time step was 0.5 ms for the probe and 2.5 ms for scanning.
  3. Refractory period: fixed at 2.1 ms. In v1.0 it was optimised within 0–1 ms.
  4. Optimisation: NSGA-II (pymoo) in Python instead of MATLAB's gamultiobj, with the same objectives, bounds, population (150) and number of generations (100). All 150 non-dominated sets of the final population are used (v1.0: 53 sets).
  5. FR and TV (dot pattern): the rate in each 10-ms window is computed by fractional-interval binning as in Connor et al. (1990), i.e. windows lying inside an interspike interval receive the fraction of that interval. In v1.0 these windows were assigned 0.
  6. All RA results in `LIF_model/*/result1/` and `result2/` were replaced. The control analysis (Sections 5.5 and 5.6) uses parameter set #87 (v1.0: #13). The SA1 results are unchanged.

  The results of the revised version were computed with the Python implementation (`python/`, see `python/README.md`). The MATLAB scripts in `LIF_model/` were updated to the same settings; `LIF_model/check_with_python.m` compares them with the Python results.


## 2. Folder Structure

```text
Project_Root/
├── FE_model/
│  ├── stimulus1_probe/
│  ├── stimulus2_braille/
│  ├── stimulus3_dot_texture/
│  └── stimulus4_single_dot/
│
├── LIF_model/
│  ├── stimulus1_probe/
│  ├── stimulus2_braille/
│  ├── stimulus3_dot_texture/
│  └── stimulus4_single_dot/
│
└── python/
    └── results/
```


### Folder Description

FE_model/
Abaqus finite element model files (`.inp`). Each `.inp` file includes nodal coordinates and element connectivity data for FE geometry models.

LIF_model/
MATLAB scripts (.m) implementing the neural response (LIF) model and corresponding simulation output for each simulation condition.

stimulusXX_XX/
Each folder corresponds to one simulation condition described in the manuscript.

python/
Python scripts (.py) implementing the same LIF model, the multi-objective optimisation and all analyses of the revised version. `python/results/` contains the parameter sets and per-set summaries (see `python/README.md`).


## 3. Software Requirements
Finite Element Simulation
Abaqus Version: Abaqus 2016 and Abaqus 2022

Neural Simulation
Python 3.14 with numpy, scipy, pandas, openpyxl and pymoo (versions in `python/requirements.txt`)
MATLAB Version: R2024b
Required toolboxes: Global Optimization Toolbox (optimisation), Signal Processing Toolbox (Sections 5.5 and 5.6)


## 4. Overall Workflow
The simulation pipeline consists of the following steps:
1. Run Abaqus FE simulation (.inp)
2. Extract strain enelgy density data from FE output (extracted data: abaqusXX.rpt and fXX-X.xlsx, stored in LIF_model/stimulusXX_XX/fingerXX/)
3. Run the LIF simulation using the extracted FE data as input: Python (`python/run_all.py` computes all results below; see `python/README.md`) or MATLAB (.m, Section 5)
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
- In the revised version, the 150 parameter sets in param.mat (and `python/results/unitA/params.csv`) were obtained with `python/optimize.py` (NSGA-II). `multiObj_Pareto.m` optimises the same objectives with gamultiobj, which gives a different front.

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
&nbsp;&nbsp;&nbsp;   `allSpikeTime`: Calculated spike timing. Exported as dot_texture_sep_param_87_ds_XX-XX_delay_XX.csv (stored in result1/). Used to generate Supplementary Fig. S7.<br>
&nbsp;&nbsp;&nbsp;   `setFiringRateRMSE`: Calculated RMSE for firing rate. Exported as dot_texture_firerate_rmse.csv (stored in result2/). Used to generate Fig. 12(A).<br>
&nbsp;&nbsp;&nbsp;   `setTempoVariRMSE`: Calculated RMSE for temporal variation. Exported as dot_texture_tempvar_rmse.csv (stored in result2/). Used to generate Fig. 12(A).


### 5.6 Effect of unresponsive period (Stimulus4)
- Navigate to:<br>
&nbsp;&nbsp;&nbsp;   LIF_model/stimulus4_single_dot_delay/
- The SED outputs from the FE model are stored in `dataSED.mat`.
- Run:<br>
&nbsp;&nbsp;&nbsp;   `mainSingleDotRA.m` (for RA responses)
- Output (variables):<br>
&nbsp;&nbsp;&nbsp;   `allSpikeTime`: Calculated spike timing. Exported as single_dot_sep_param_87_dw_XX-XX_delay_XX.csv (stored in result1/). Used to generate Supplementary Fig. S8.<br>
&nbsp;&nbsp;&nbsp;   `setRespWidth`: Calculated response width. Exported as single_dot_respwidth.csv (stored in result2/). The response widths in XX_respwidth2.csv represent the differences between the simulation results and the experimental data. Used to generate Fig. 12(B).


### 5.7 Effect of neuronal noise (Stimulus3)
- Navigate to:<br>
&nbsp;&nbsp;&nbsp;   LIF_model/stimulus3_dot_texture_noise/
- The SED outputs from the FE model are stored in `dataSED.mat`.
- Run:<br>
&nbsp;&nbsp;&nbsp;   `mainDotTexRA.m` (for RA responses)
- Output (variables):<br>
&nbsp;&nbsp;&nbsp;   `allSpikeTime`: Calculated spike timing. Exported as dot_texture_sep_param_XX_ds_XX-XX_noise_XX.csv (stored in result1/). Used to generate Fig. 13.<br>
&nbsp;&nbsp;&nbsp;   `setFiringRateRMSE`: Calculated RMSE for firing rate. Exported as dot_texture_firerate_rmse.csv (stored in result2/). The results in XX_rmse.csv are sorted in ascending order of the parameter set number, whereas the results in XX_rmse2.csv are sorted in ascending order of &Delta;FST. Supplementary Fig. S11(A) was generated using XX_rmse2.csv.<br>
&nbsp;&nbsp;&nbsp;   `setTempoVariRMSE`: Calculated RMSE for temporal variation. Exported as dot_texture_tempvar_rmse.csv (stored in result2/). Supplementary Fig. S11(B) was generated using XX_rmse2.csv.


### 5.8 Effect of neuronal noise (Stimulus4)
- Navigate to:<br>
&nbsp;&nbsp;&nbsp;   LIF_model/stimulus4_single_dot_noise/
- The SED outputs from the FE model are stored in `dataSED.mat`.
- Run:<br>
&nbsp;&nbsp;&nbsp;   `mainSingleDotRA.m` (for RA responses)
- Output (variables):<br>
&nbsp;&nbsp;&nbsp;   `allSpikeTime`: Calculated spike timing. Exported as single_dot_sep_param_XX_dw_XX-XX_noise_XX.csv (stored in result1/). Used to generate Supplementary Fig. S9.<br>
&nbsp;&nbsp;&nbsp;   `setRespWidth`: Calculated response width. Exported as single_dot_respwidth.csv (stored in result2/). The results in XX_respwidth.csv are sorted in ascending order of the parameter set number, whereas the results in XX_respwidth2.csv are sorted in ascending order of &Delta;FST. In addition, the response widths in XX_respwidth2.csv represent the differences between the simulation results and the experimental data. Supplementary Fig. S10 was generated using XX_respwidth2.csv.


### 5.9 Effect of response width definition (Stimulus4)
- Navigate to:<br>
&nbsp;&nbsp;&nbsp;   LIF_model/stimulus4_single_dot_modi_RW/
- The SED outputs from the FE model are stored in `dataSED.mat`.
- Run:<br>
&nbsp;&nbsp;&nbsp;   `mainSingleDotRA.m` (for RA responses)
- Output (variables):<br>
&nbsp;&nbsp;&nbsp;   `allSpikeTime`: Calculated spike timing. Exported as single_dot_sep_param_XX_dw_XX-XX_noise_XX.csv (stored in result1/). Used to generate Supplementary Fig. S6.<br>
&nbsp;&nbsp;&nbsp;   `spikeWidth`: Calculated response width. Exported as single_dot_respwidth.csv (stored in result2/). The results in XX_respwidth.csv are sorted in ascending order of the parameter set number, whereas the results in XX_respwidth2.csv are sorted in ascending order of &Delta;FST.

