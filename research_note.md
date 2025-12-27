---
title: "Robust Scenario Allocator (RSA): A Generative Strategy Ensemble Framework for Strategic Asset Allocation"
author: "Vereszki Péter"
date: "2025-12-26"
output:
  html_document:
    toc: true
    number_sections: true
  pdf_document:
    toc: true
---

# 1. Executive Summary

This research note outlines the theoretical and architectural framework for the **Robust Scenario Allocator (RSA)** project. The project addresses the challenge of **Strategic Asset Allocation (SAA)** over a 5-year horizon using "Small Data" ($N \approx 220$ monthly observations).

**Paradigm Shift:** We reject the traditional "Predictive Super Learner" approach (predicting returns via regression) due to the low signal-to-noise ratio in monthly financial data. Instead, we propose a **Generative Strategy Ensemble**.
1.  **Physics:** A **Regime-Switching BVAR** and **Dynamic t-Copula** engine generates a "Virtual Future" ($P$-measure scenarios) capturing structural breaks and tail dependence.
2.  **Strategies:** A library of robust, convex base strategies (Min-CVaR, HRP, Min-CDaR) is evaluated on these scenarios.
3.  **Meta-Learning:** A **Distributionally Robust Optimization (DRO)** layer finds the optimal convex combination of these strategies, anchored by **Entropy Pooling** views.

This framework guarantees a **"Convex-Hull Oracle"** property: the resulting allocation is globally optimal within the convex hull of the base strategies under the modeled uncertainty.

---

# 2. Context and Motivation

## 2.1 The "Small Data" Problem in SAA
Strategic Asset Allocation requires stable estimates of long-term parameters. However, with only ~220 monthly observations (approx. 18 years), standard Mean-Variance optimization is notoriously unstable ("Error Maximization"). Epistemic uncertainty (parameter uncertainty) often dominates aleatoric uncertainty (market noise).

## 2.2 The Fallacy of Direct Return Prediction
Applying Super Learner (stacking) to directly predict $r_{t+1}$ using ML regressors (Random Forest, Ridge) fails in SAA because:
*   The $R^2$ of monthly returns is near zero.
*   Minimizing MSE (prediction error) does not translate to maximizing Utility or minimizing Drawdown (portfolio objectives).
*   **Correction:** We move from **Prediction** to **Allocation**. We do not predict *returns*; we predict *which allocation logic* (Strategy) performs best across a range of plausible future regimes.

---

# 3. The Generative Engine (The "Physics")

To evaluate strategies robustly, we cannot rely on historical bootstrapping (which assumes stationarity). We must generate a "Virtual Future" that respects macro-dynamics and tail risks.

## 3.1 Macro-Dynamics: Regime-Switching BVAR
We model the joint dynamics of Macro variables (Inflation, Growth, Rates) and Asset proxies using a **Bayesian Structural VAR with Markov Switching Heteroskedasticity**.

$$ Y_t = \mu_{S_t} + \sum_{p=1}^P A_{p, S_t} Y_{t-p} + \epsilon_t, \quad \epsilon_t \sim N(0, \Sigma_{S_t}) $$

*   **Regimes ($S_t$):** Captures distinct states (e.g., "Great Moderation" vs. "Stagflation").
*   **P-Measure:** The BVAR is trained on history, thus producing real-world probability paths, crucial for SAA.
*   **Implementation:** `bsvars` (R package) for Structural BVARs.

## 3.2 Dependence Structure: Dynamic t-Copula
Linear correlation fails to capture **Tail Dependence** (the tendency of assets to crash together). We separate marginals from dependence via Sklar's Theorem:
$$ F(x_1, ..., x_d) = C(F_1(x_1), ..., F_d(x_d)) $$

We employ a **DCC-GARCH** model with a multivariate **Student-t Copula**:
*   **Marginals:** AR(1)-GJR-GARCH(1,1) with Skewed-t distribution (captures volatility clustering and asymmetry).
*   **Copula:** Time-varying correlation matrix $R_t$ and degrees of freedom $\nu$ (captures joint tail risk).
*   **Implementation:** `rmgarch` (R package) by Alexios Ghalanos.

---

# 4. Base Learners: The Convex Strategy Library

The "Base Learners" are defined as autonomous allocation algorithms $S_k: \mathcal{I}_t \to w_t$. We select strategies that are **convex** (globally solvable via LP/QP) or **deterministic** to ensure stability.

## 4.1 Min-CVaR (Minimum Expected Shortfall)
**Logic:** Minimizes the average loss in the worst $\alpha\%$ of cases.
**Formulation (LP):**
$$ \min_{w, \zeta, u} \zeta + \frac{1}{(1-\alpha)T} \sum_{t=1}^T u_t $$
$$ \text{s.t. } u_t \ge -w^T r_t - \zeta, \quad u_t \ge 0, \quad w \in \mathcal{W} $$
**Implementation:** `PortfolioAnalytics` (using `ROI` solvers).

## 4.2 Min-CDaR (Conditional Drawdown at Risk)
**Logic:** Minimizes the tail of the drawdown distribution. More relevant for SAA than volatility.
**Formulation (LP):** Similar to CVaR, but operating on the cumulative drawdown function $D_t(w)$.
**Implementation:** `PerformanceAnalytics` / Custom LP via `CVXR`.

## 4.3 HRP (Hierarchical Risk Parity)
**Logic:** Applies graph theory (clustering) to handle covariance instability. It does not require matrix inversion.
1.  **Tree Clustering:** Distance metric $d_{ij} = \sqrt{2(1-\rho_{ij})}$.
2.  **Quasi-Diagonalization:** Reordering the covariance matrix.
3.  **Recursive Bisection:** Allocating based on cluster variance.
**Implementation:** `HierPortfolios` or `ClusterPortfolios`.

## 4.4 ERC (Equal Risk Contribution)
**Logic:** Finds $w$ such that $w_i (\Sigma w)_i = w_j (\Sigma w)_j$. Anchors the portfolio in maximum diversification.
**Implementation:** `riskParityPortfolio`.

## 4.5 Robust Mean-Variance (Max Sharpe)
**Logic:** The classic approach, but strictly using **Ledoit-Wolf Shrinkage** for $\Sigma$ to mitigate estimation error.
**Formulation (QP):** $\max w^T \mu - \frac{\gamma}{2} w^T \Sigma_{LW} w$.
**Implementation:** `PortfolioAnalytics` (demo_max_Sharpe).

---

# 5. Views and Anchoring: Entropy Pooling

To incorporate expert views (e.g., Term Premium logic) and anchor the generative scenarios, we use **Entropy Pooling (EP)**.

## 5.1 The Principle
We seek a new probability vector $\tilde{p}$ that satisfies linear constraints ($A\tilde{p} \le b$) while remaining as close as possible to the prior $p$ (from the BVAR/Copula).
**Objective (Min Relative Entropy / KL Divergence):**
$$ \tilde{p} = \arg\min_{x} \sum_{j=1}^J x_j \ln \left( \frac{x_j}{p_j} \right) $$

## 5.2 Vorobets Sequential Entropy Pooling (SeqEP)
Standard EP on the full horizon breaks time-consistency (martingale property). **SeqEP** decomposes the problem via the Chain Rule for Relative Entropy, applying views sequentially at each time step $t$.

**Key View applied in RSA:**
*   **Term Premium:** $E^{\tilde{p}}[Yield_{10y}] - E^{\tilde{p}}[ShortRate] \in [\delta_{min}, \delta_{max}]$.
*   **Implementation:** `CVXR` using the `ECOS` solver (Exponential Cone support).

---

# 6. The Meta-Learner: Distributionally Robust Optimization

The "Super Learner" is formulated as a portfolio optimization problem over the equity curves of the Base Strategies.

## 6.1 The Meta-Problem
Let $R_{k, s}$ be the cumulative return of Strategy $k$ in Scenario $s$ (weighted by $\tilde{p}_s$). We seek meta-weights $\theta$.

Instead of optimizing for the single distribution $\tilde{p}$, we use **DRO** to optimize for the worst-case distribution $Q$ within a **Wasserstein Ambiguity Set** $\mathcal{U}_\epsilon(\tilde{p})$.

$$ \max_{\theta} \min_{Q \in \mathcal{U}_\epsilon(\tilde{p})} E_Q \left[ U \left( \sum_{k=1}^K \theta_k R_{k} \right) \right] $$

## 6.2 Convex-Hull Oracle Property
Since the Base Strategies are convex (LP/QP) and the Meta-Learner is a convex optimization (DRO/SOCP), the final solution is guaranteed to be a **global optimum** within the feasible region defined by the base strategies. There are no local minima traps common in non-convex Neural Networks.

## 6.3 Implementation
*   **Solver:** `CVXR` (SOCP mode).
*   **Regularization:** We apply an entropy penalty on $\theta$ to prevent the Meta-Learner from putting 100% weight on a single strategy ("Winner-takes-all" instability).

---

# 7. Validation Framework

Since we optimize on generated scenarios, we need rigorous OOS validation on historical data.

1.  **Walk-Forward Validation:** Expanding window (Initial: 60 months).
2.  **Metrics:**
    *   **CDaR (95%):** Primary risk metric.
    *   **Probabilistic Sharpe Ratio (PSR):** Adjusts Sharpe for non-normality and track record length.
    *   **Deflated Sharpe Ratio (DSR):** Adjusts for "Multiple Testing Bias" (selection bias among strategies).
3.  **Sanity Checks:**
    *   **Turnover:** Is the strategy trading too much?
    *   **Leverage:** Is it implicitly leveraging via the strategies?

---

# 8. Selected Bibliography & R Packages

*   **Super Learner:** Polley, E. C., & van der Laan, M. J. (2010). *Super Learner In Prediction*.
*   **Entropy Pooling:** Meucci, A. (2008). *Fully Flexible Views: Theory and Practice*.
*   **Sequential EP:** Vorobets, T. (2021). *Sequential Entropy Pooling*.
*   **DCC-GARCH:** Ghalanos, A. (2019). *rmgarch: Multivariate GARCH models*.
*   **HRP:** López de Prado, M. (2016). *Building Diversified Portfolios that Outperform Out-of-Sample*.
*   **BVARs:** Woźniak, T. (2022). *bsvars: Bayesian Structural Vector Autoregressions*.
*   **Optimization:** Fu, A., Narasimhan, B., & Boyd, S. (2020). *CVXR: An R Package for Disciplined Convex Optimization*.