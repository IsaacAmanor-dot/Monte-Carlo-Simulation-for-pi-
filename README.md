# Monte Carlo Estimation of π  
AMS 595 / DCS 525 — Project 1  

## Overview
This project implements three approaches to estimating π using Monte Carlo simulation in MATLAB.  
The algorithm relies on the fact that the ratio of points falling inside a unit quarter-circle to those in a unit square approximates π/4. Scaling this ratio yields an estimate for π.  

The work demonstrates:
1. **Fixed-sample simulation (for-loop)** — precision and cost analysis.  
2. **Adaptive simulation (while-loop)** — automatic stopping once a target number of significant figures is achieved.  
3. **Interactive function with visualization** — live display of sampling and auto-stop at user-specified precision.  

---

## Tasks

### Task 1: For-loop with fixed N
- Generates random points for increasing $N$.  
- Computes $\hat{\pi}$, absolute error, and runtime.  
- Produces plots:  
  - π vs N  
  - Absolute error vs N (log-log)  
  - Runtime vs N (log-log)  
  - Absolute error vs runtime (precision–cost tradeoff)  

### Task 2: While-loop with stability-based auto-stop
- Uses streaming batches of random points.  
- Halts when successive estimates of π agree to **s significant figures** for a fixed number of checks.  
- Does **not** use the true value of π for stopping.  
- Outputs results to CSV with columns:  
  `s_sigfigs, N_used, pi_hat, time_s, iters_batches, checks`.  

### Task 3: Function with live visualization
- User inputs target precision `s`.  
- While-loop sampling with stability-based auto-stop.  
- Interactive plot:
  - Blue = points inside circle  
  - Red = points outside  
  - Overlay shows current N, running estimate, and stabilized π value  
- Returns estimated π to the requested precision.  

---



**Example Results**

| Estimated π vs N | Execution Time vs N |
|------------------|----------------------|
| ![π vs N](figures/fig_task1_min_pi_vs_N.pdf) | ![Execution Time vs N](figures/fig_task1_min_time_vs_N.pdf) |


Command window:
