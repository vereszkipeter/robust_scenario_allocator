# GEMINI.md - Project Context & Architecture Guide (v2.0)

**Project Name:** Robust Scenario Allocator (RSA)
**Domain:** Strategic Asset Allocation (SAA), Quantitative Finance, Macro-Econometrics
**Tech Stack:** R (Pure), `targets`, `bvars`, `rmgarch`, `CVXR`, `PortfolioAnalytics`, `HierPortfolios`

---

## 1. Executive Summary & Paradigm Shift

This project implements a **Generative Strategy Ensemble** system for a 5-year ($T=60$ months) SAA horizon.

**Core Philosophy (v2.0):**
We have moved away from "Predicting Returns" (which fails on small $N=220$ data) to **"Allocating to Strategies"**.
Instead of asking the Super Learner to predict the next month's return using Regression, we ask it to find the optimal convex combination of robust **Base Strategies** (e.g., Min-CDaR, HRP, Risk Parity) evaluated over **Generative Scenarios** that contain macro-structural breaks.

**The Equation:**
$$ \text{Final Portfolio} = \text{Meta-Optimization}(\text{Base Strategies}(\text{History})) \times \text{Evaluated on}(\text{RS-BVAR Scenarios}) $$

---

## 2. Quantitative Architecture

The pipeline consists of 5 logical layers, moving from Physics (Data) to Decisions.

### Layer 1: The Generative Engine (The Physics)
This layer creates the "Virtual Future" ($X_{scen}$) which includes regimes not necessarily seen in the immediate history.
*   **Macro-Economy:** **Regime-Switching BVAR (RS-BVAR)** (`bvars` package). Captures $P$-measure dynamics of Growth, Inflation, and Rates.
*   **Dependence:** **Dynamic t-Copula (DCC-GARCH)** (`rmgarch` package). Separates marginal volatility from tail dependence.
*   **Output:** A tensor of 10,000 physical scenario paths for 60 months for all assets.

### Layer 2: The Base Learners (The Strategies)
These are *not* regression models. They are autonomous allocation algorithms that output weights $w_k$.
*   **Implementation:** Using `PortfolioAnalytics` and `HierPortfolios`.
*   **Strategy 1: Min-CVaR (LP):** Minimizes Expected Shortfall. Convex, robust to fat tails.
*   **Strategy 2: Min-CDaR (LP):** Minimizes Conditional Drawdown at Risk. Focuses on path-dependent loss.
*   **Strategy 3: HRP/HERC:** Hierarchical Risk Parity. Ignores expected returns, focuses on clustering structure. Robust to covariance noise.
*   **Strategy 4: ERC:** Equal Risk Contribution. The "Diversification Anchor".
*   **Strategy 5: Max-Sharpe (QP):** Traditional MV, but strictly with **Ledoit-Wolf shrinkage** for stability.

### Layer 3: Anchoring & Views (The Future Adjustment)
We refine the *probability* of the generated scenarios using **Meucci’s Entropy Pooling**.
*   **Method:** **Vorobets Sequential Entropy Pooling (SeqEP)**.
*   **Views:**
    *   **Term Premium:** Explicitly modeled as $TP = Yield_{model} - E^P[ShortRate]$. View: $TP \in [-50bps, +150bps]$.
    *   **Equilibrium:** Implied returns from ERC are used as a reference point for the "Normal" regime.
*   **Solver:** `CVXR` (ECOS solver for Exponential Cone).
*   **Result:** A probability vector $p_{post}$ for the 10,000 scenarios.

### Layer 4: The Meta-Learner (The Allocator)
This is the "Super Learner" step.
*   **Input:** The cumulative P&L of the *Base Strategies* if they were held through the *Generative Scenarios*.
*   **Optimization:** A convex optimization problem (LP or QP) over the *strategies*.
    *   Objective: Maximize Utility or Minimize Meta-CVaR.
    *   Constraints: Convex combination ($\sum \theta_k = 1, \theta_k \ge 0$).
    *   Regularization: Entropy penalty on weights ($\sum \theta \ln \theta$) to prevent concentration in one strategy.
*   **Robustness:** **Distributionally Robust Optimization (DRO)**. We optimize for the worst-case distribution within a Wasserstein ball around the Entropy-Pooled probabilities.

---

## 3. Technical Implementation Guidelines (R)

### Constraints
*   **Pure R:** No Python.
*   **Pipeline:** Use `targets`.
*   **Reproducibility:** `renv`.

### Key Functionalities
1.  **`R/02_models.R`**: Must handle `bsvars` for MS-BVAR and `rmgarch` for DCC.
2.  **`R/03_strategies.R`**: New file! Needs to wrap `PortfolioAnalytics` functions into standard "Base Learners" that accept history and return weights.
3.  **`R/04_scenarios.R`**: Generates the 3D tensor.
4.  **`R/05_views.R`**: Implements Vorobets SeqEP via `CVXR`.
5.  **`R/06_meta.R`**: Performs the DRO optimization over the strategies.

### Mathematical "Gotchas" (Context for Code Assist)
1.  **P vs Q Measure:** The BVAR generates P-measure (Real World) rates. The "Risk Free" asset in the portfolio is a rolling reinvestment strategy at these rates, *not* a static yield.
2.  **Zero-Mean Trap:** HRP/ERC ignores returns. When combining them in the Meta-Learner, we must evaluate them on scenarios that *have* drifts (from the BVAR), otherwise, the Meta-Learner can't distinguish profitable from unprofitable strategies.
3.  **Look-Ahead Bias:** Base strategies must calculate weights using *only* data available at $t=0$ (or expanding window history). They cannot "see" the BVAR scenarios. Only the Meta-Learner "sees" the scenarios to judge the strategies.

---

## 4. Roadmap for the Agent

1.  **Data:** Load Macro + Asset data (Lag macro by 1 month!).
2.  **Physics:** Fit RS-BVAR and DCC-Copula. Generate `scenarios.rds`.
3.  **Strategies:** Calculate current weights $w_{k}$ for all 5 Base Strategies based on historical data.
4.  **Backtest-on-Future:** Multiply $w_{k}$ by the `scenarios` to get 5 sets of 10,000 future equity curves.
5.  **Meta-Optimize:** Run DRO/Entropy Pooling on these 5 equity curves to find the optimal $\theta$ mix.
6.  **Report:** Plot the Fan Chart of the final RSA portfolio.