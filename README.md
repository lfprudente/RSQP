# **Sequential Quadratic Programming on Riemannian Manifolds**

This repository contains MATLAB codes related to the numerical experiments from the article:

>  O. P Ferreira, L. F. Prudente, Jen-Chih Yao, Xiaopeng Zhao, Sequential quadratic programming framework for nonlinear optimization on Riemannian manifolds, 2025.

**Date:** October 2025

---

## 📂 **Contents**

### Algorithms
- **`RSQP.m`** — RSQP with PHR (LStype=1) or ℓ₁ (LStype=2) merit function.  
- **`auglag.m`** — Riemannian Augmented Lagrangian method using Manopt’s solvers: trust-regions (InSolver=1) or Riemannian BFGS (InSolver=2).

### RSQP Core Routines
- **`build_subproblem.m`** — Builds the reduced quadratic subproblem in local tangent coordinates.  
- **`make_SPD.m`** — Safeguards the reduced Hessian by ensuring positive definiteness.


### Problem Routines
- **`evalf.m`** — Objective function.  
- **`evalg.m`** — Euclidean gradient of the objective function.  
- **`evalh.m`** — Euclidean Hessian of the objective function.  
- **`evalcc.m`** — Constraints.  
- **`evalnc.m`** — Gradients of the constraints.  
- **`evalhc.m`** — Hessians of the constraints.  

### Main Script
- **`main.m`** — Controls execution of the numerical experiments. The user selects the problem (nnPCA, packing, or classification) and the algorithm (RSQP-PHR, RSQP-ℓ₁, AL-TR, or AL-BFGS).

---
## ✅ Requirements

- MATLAB R2024b (24.2.0.2712019) or later.  
- Optimization Toolbox (for `quadprog`).  
- [**Manopt 8.0**](https://www.manopt.org/) — MATLAB toolbox for optimization on manifolds (GPLv3, Nicolas Boumal).
---

## 🚀 Quick Start

1. Open MATLAB in the main repository folder.  
2. Install dependencies by running:  
   ```matlab
   install_dependencies
3. Run the main script: 
    ```matlab
    main
---

## 📦 **Third-Party Codes**

This repository includes third-party free software:

- **Manopt 8.0**  
    - **Description:** MATLAB toolbox for optimization on manifolds.  
    - **Author:** Nicolas Boumal  
    - **Website:** [https://www.manopt.org/](https://www.manopt.org/)  
    - **License:** GNU General Public License (GPL) version 3  

⚠️ **Note:** Please review the specific licenses of the third-party codes before modifying or redistributing them. License information can be found in code comments or separate text files in their corresponding directories.

---

## 📄 **License**

This project is licensed under the **GNU General Public License (GPL) version 3**. For more details, see the **LICENSE** file included in this repository.

---
