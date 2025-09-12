# Monte Carlo Estimation of π  
AMS 595 / DCS 525 — Project 1  

## Overview
This repository explores how Monte Carlo methods can be used to approximate π through random sampling.  
By generating points inside a unit square and counting how many fall within the quarter circle, we obtain a probabilistic estimator for π.  

The project demonstrates three complementary approaches:
1. **Fixed-sample simulation (for-loop):** analyze accuracy, error, and runtime tradeoffs.  
2. **Adaptive simulation (while-loop):** stop automatically once a desired precision is achieved.  
3. **Interactive simulation (function with visualization):** watch π emerge in real time with auto-stop at user-specified precision.  

---

## Tasks

### 🔹 Task 1: Fixed-Sample Monte Carlo (For-Loop)
- Runs simulations for increasing values of $N$.  
- Computes $\hat{\pi}$, absolute error, and execution time.  
- Produces key performance plots:  
  - π vs N  
  - Absolute error vs N (log-log)  
  - Execution time vs N (log-log)  
  - Error vs runtime (precision–cost tradeoff)  

### 🔹 Task 2: Adaptive Monte Carlo (While-Loop with Auto-Stop)
- Streams points in **batches** to remain memory-efficient.  
- Stops automatically when $\hat{\pi}$ stabilizes to a user-defined number of **significant figures $s$**.  
- Does **not** rely on the true value of π for stopping.  
- Saves results to CSV with fields:  
  `$s$, $N_{\text{used}}$, $\hat{\pi}$, runtime (s), batches, checks`.  

### 🔹 Task 3: Interactive Monte Carlo (Function with Visualization)
- User specifies target precision $s$.  
- While-loop sampling with stability-based auto-stop.  
- Live visualization:  
  - **Blue points** = inside the circle  
  - **Red points** = outside the circle  
  - On-screen overlay: $N$, running $\hat{\pi}$, and stabilized value  
- Returns the final estimate of π at the requested precision.  

---

✨ This project highlights both the **power** and the **limitations** of Monte Carlo methods: they are simple and flexible, yet require large sample sizes to achieve high precision. Through adaptive stopping and visualization, the simulations balance numerical rigor with interpretability.  
