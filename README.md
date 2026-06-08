# Online MPLS + Fading SCI — Cooperative Localization Simulation Code

MATLAB simulation code accompanying the paper:

> **Dynamic Ranging-Error Compensation and Consistent Cooperative Localization in TDMA UWB UAV Swarms**

It implements and reproduces all simulation results of a two-module decentralized
cooperative-localization architecture for high-dynamic time-division UWB UAV swarms:

- an **online MPLS ranging front-end** — a sliding-window, Legendre-basis, IRLS-Tukey
  extension of the batch polynomial joint ranging/synchronization model — that outputs,
  per link, a bias-corrected range and a calibrated time-varying variance; and
- a **fading SCI cooperative back-end** — a Split Covariance Intersection filter with
  continuous-time forgetting on the dependent covariance and a determinant-criterion
  mixing-weight search — that consumes that variance through a per-link variance interface.

A five-method ablation (SCI+MPLS / EKF+MPLS / EKF+Inflate+MPLS / CI+MPLS / SDS-TWR+SCI)
isolates front-end ranging quality (drives RMSE) from back-end correlation handling
(drives estimator consistency, ANEES).

---

## 1. Requirements

- **MATLAB R2025a** (any R2019b+ should run; `tiledlayout` is used for the figures).
- **Statistics and Machine Learning Toolbox** — `signrank`, percentile/bootstrap helpers.
- **Parallel Computing Toolbox** — *optional*; the batch scripts use `parfor` and fall back
  to serial execution if no pool is available.

No external dependencies. Output folders are created automatically on first run.

---

## 2. Scenarios

Geometry is defined in `set_parameters.m`, selected by `param.scenario`:

| `scenario` | Experiment | Nodes | Description |
|---|---|---|---|
| `1` | **E1** | 2  | Two nodes on counter-rotating circular orbits; pure ranging validation. |
| `2` | **E2 / E3** | 16 | Heterogeneous cluster + patrol swarm. `param.anchors_mobile = false` → **E2** (2 fixed anchors + 14 mobile); `= true` → **E3** (all 16 mobile, anchor-free). |

Key defaults: polynomial order `L = 3`; MPLS window `= 6` polling rounds; fading time
constant `tau = 2.0 s`; odometry drift `kappa_at = 0.04`, `kappa_ct = 0.075` sqrt(m);
one SCI update per polling round; LoS propagation (no NLOS).

---

## 3. File map

The **runnable scripts live at the top level**; the functions they call live in the
**`lib/`** subfolder. Each top-level script adds `lib/` to the MATLAB path automatically,
so there is no need to `addpath` by hand — just run the scripts from this directory.

**`lib/` — function library** (called by the drivers; not run directly)

| File | Role |
|---|---|
| `lib/set_parameters.m` | Scenario geometry, physics/clock/MAC/odometry defaults. |
| `lib/get_Xtrue.m` | Ground-truth trajectory generator. |
| `lib/mac_protocol_new.m` | TDMA MAC + two-way-ranging timestamp generation. |
| `lib/MPLS.m` | Single-link polynomial least squares (Legendre basis, IRLS-Tukey). |
| `lib/MMPLS.m` | Two-pass network polynomial LS (clock fusion, then distance). |
| `lib/MMPLS_analysis_function.m` | Front-end entry: sliding window + Raw/ZOH/SDS-TWR baselines; saves per-link calibrated ranges/variances. |
| `lib/SCI_Main_Using_MPLS_function_new.m` | Anchored back-end (E2): SCI / EKF / EKF+Inflate / CI. |
| `lib/SCI_Main_AnchorFree_function.m` | Anchor-free back-end (E3): same methods + dead-reckoning baseline. |
| `lib/print_param_table.m` | Parameter-summary printer (called by the batches). |

**Experiment drivers** (entry points — run these)

| File | Produces |
|---|---|
| `Run_Batch_E1_Ranging_Test.m` | E1 ranging accuracy + variance calibration. |
| `Run_Batch_E2_Anchored_Test.m` | E2 anchored cooperative localization (5-way). |
| `Run_Batch_E3_AnchorFree_Test.m` | E3 anchor-free cooperative localization (5-way + DR). |
| `Run_Batch_E2_PacketLoss_Test.m` | Packet-loss robustness sweep (i.i.d. 0–20%). |
| `Validate_Tau_Sweep.m` | Fading-memory T_round-invariance check (continuous-time vs fixed per-step). |
| `Run_Ablation_Lambda.m` | Fading time-constant `tau` sensitivity sweep. |

**Plotting / statistics** (regenerate paper figures and tests from a saved `.mat`)

`Plot_E1_FromMat.m`, `Plot_E2_FromMat.m`, `Plot_E3_FromMat.m`,
`Plot_E2_PacketLoss_Median.m`, `Stat_Wilcoxon_FromMat.m` (pairwise Wilcoxon
signed-rank, Holm–Bonferroni, rank-biserial), and `Measure_Timing.m` (MPLS / SCI
wall-clock timing). These are self-contained — they read a saved `.mat` (or, for
`Measure_Timing`, nothing) and do not depend on `lib/`.

---

## 4. How to run

Run all scripts **from this directory** so the relative output paths resolve.

### Option A — regenerate figures from bundled data (fast, no simulation)

Pre-computed `N_sim = 100` results are bundled in `matfile/` (the `*_sample.mat` files).
To reproduce the paper figures/tables **without** re-running the heavy Monte Carlo, run the
plot / statistics scripts directly:

```matlab
Plot_E1_FromMat                % E1 distance RMSE + variance calibration
Plot_E2_FromMat                % E2 RMSE + ANEES (2x2)
Plot_E3_FromMat                % E3 RMSE + ANEES + DR (2x2)
Plot_E2_PacketLoss_Median      % packet-loss sensitivity
Stat_Wilcoxon_FromMat          % significance TABLES (console + LaTeX, no figure); reads E2 + E3
```

`Stat_Wilcoxon_FromMat` is the only analysis script that produces *tables* rather than a
figure: it runs the pairwise Wilcoxon signed-rank tests behind the paper's significance
claims and prints the corresponding LaTeX table snippets (it needs both the E2 and E3 results).

Each plot/stat script auto-selects the **latest** matching `.mat` in `matfile/` — the
bundled sample if you have not run a batch yourself.

### Option B — full reproduction (re-run the Monte Carlo)

Each of the four E1/E2/E3/packet-loss batches is self-contained: it runs the full Monte
Carlo **and regenerates its figures at the end** by calling the matching `Plot_*` script
internally, so **no separate plotting step is needed** for the figures. The one analysis the
batches do *not* auto-run is the cross-experiment significance table: `Stat_Wilcoxon_FromMat`
reads both the E2 and E3 results, so run it after those two batches. `Validate_Tau_Sweep` and
`Run_Ablation_Lambda` are self-contained analysis scripts that produce their own figures.

```matlab
Run_Batch_E1_Ranging_Test      % E1
Run_Batch_E2_Anchored_Test     % E2
Run_Batch_E3_AnchorFree_Test   % E3
Run_Batch_E2_PacketLoss_Test   % packet loss
Validate_Tau_Sweep             % fading-memory robustness (paper Sec. 4.5)
Run_Ablation_Lambda            % tau sensitivity (paper Sec. 4.5)
Stat_Wilcoxon_FromMat          % significance tables (run after the E2 + E3 batches)
```

A full 16-node, `N_sim = 100`, two-phase batch is computationally heavy (use the Parallel
Computing Toolbox). For a quick check, reduce `N_sim` at the top of the batch script (the
sweep arrays can stay).

### Output locations

Every run writes timestamped (`yy-mm-dd-HH-MM`) outputs under this directory:

| Folder | Contents |
|---|---|
| `matfile/` | `results_<exp>_<case>_<timestamp>.mat` — raw results, re-loaded by the plot/stat scripts. |
| `photo/<Exp>/<timestamp>/` | Figures (`.pdf` + `.png`). |
| `csv/<Exp>/<timestamp>/` | Per-condition summary tables. |
| `json/<Exp>/<timestamp>/` | Machine-readable results + run metadata. |

`<Exp>` is `E1_RangingTest`, `E2_AnchoredTest`, `E3_AnchorFreeTest`, or `E2_PacketLoss`.
Only `matfile/*_sample.mat` is tracked in git; all other generated output is ignored.

---

## 5. Reproducibility notes

- **Monte Carlo:** `N_sim = 100` trials per operating point, each seeded by its index; all
  compared methods are reset to a common per-trial seed before fusion, so within a trial
  they share identical trajectory/clock/noise realizations (this makes the Wilcoxon tests
  properly paired).
- **Sweeps:** two-phase — Phase 1 sweeps velocity `v ∈ {5,…,40}` m/s at `T_round = 64` ms;
  Phase 2 sweeps `T_round ∈ {16,…,160}` ms at `v = 25` m/s; cross-point `(25 m/s, 64 ms)`.
  The `v = 5` point is simulated but excluded from the reported envelope (initialization
  transient dominates).
- **Metrics:** evaluated over `[0.4 s, T_sim]`, one sample per TDMA slot; ANEES reported as
  the cross-trial median (the per-trial distribution is right-skewed), RMSE as the mean.
- **LoS only.** The paper assumes line-of-sight; NLOS is out of scope and disabled.
- The back-end functions also contain flag-gated legacy/ablation branches (an alternative
  noise split, mixing-weight criterion, and forgetting form). These are not dead code — the
  ablation scripts exercise them as reference baselines (`Run_Ablation_Lambda` runs an
  `A_old` legacy core; `Validate_Tau_Sweep` compares per-step against continuous forgetting).
  The production configuration follows Pierre's structural noise split (Eq. 21) and
  determinant mixing criterion with continuous-time `tau` forgetting.

---

## 6. License and citation

Released under the MIT License (see `LICENSE`). If you use this code, please cite the
accompanying paper and the archived repository; machine-readable metadata is in `CITATION.cff`.
