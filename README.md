# Robust Scenario Allocator (RSA)

[![License: CC BY-NC-SA 4.0](https://img.shields.io/badge/License-CC%20BY--NC--SA%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by-nc-sa/4.0/)

Ez a projekt egy **Stratégiai Eszközallokációs (SAA)** keretrendszer, amely a modern kvantitatív pénzügyek és a gépi tanulás robusztus módszereit ötvözi R környezetben. A rendszer nem hozamokat próbál jósolni (ami kis adaton instabil), hanem **robosztus stratégiákat versenyeztet** generatív makrogazdasági szcenáriókon.

## 🚀 Projekt Filozófia

A **Robust Scenario Allocator** szakít a hagyományos "átlag-szórás" modellekkel és a naiv gépi tanulással.
*   **Generatív, nem Prediktív:** A jövőt nem egyetlen pontbecsléssel írjuk le, hanem egy Regime-Switching BVAR és Copula által generált, 10,000 elemű szcenárió-térrel.
*   **Stratégia Ensemble:** A "Base Learnerek" nem regressziós modellek, hanem teljes értékű allokációs algoritmusok (pl. HRP, Min-CDaR).
*   **Matematikai Robusztusság:** A döntési rétegben Distributionally Robust Optimization (DRO) és Entropy Pooling gondoskodik arról, hogy a modell kezelje a bizonytalanságot.

## 🏗 Architektúra és Módszertan

A rendszer 5 logikai rétegből épül fel:

1.  **Generatív Motor (Physics):**
    *   **RS-BVAR:** Makrogazdasági pályák (GDP, Infláció, Kamatok) generálása a P-mérték alatt.
    *   **Dynamic t-Copula:** A farok-kockázatok és változó korrelációk modellezése.
2.  **Base Strategies (The Players):**
    *   **Min-CVaR & Min-CDaR:** Konvex optimalizálás (LP) a tail risk és drawdown minimalizálására.
    *   **HRP & ERC:** Hierarchikus és Kockázati Paritás alapú diverzifikáció.
    *   **Shrunk Mean-Variance:** Klasszikus megközelítés Ledoit-Wolf zsugorítással.
3.  **Views & Anchoring:**
    *   **Entropy Pooling (SeqEP):** A szcenáriók valószínűségeinek finomhangolása (pl. Term Premium nézetek beépítése).
    *   **Black-Litterman Anchor:** Egyensúlyi hozamok használata referenciaként.
4.  **Meta-Learner (The Judge):**
    *   Egy **konvex optimalizáló**, amely meghatározza a stratégiák optimális súlyozását a generált szcenáriók alapján.
    *   Célfüggvény: Robusztus hasznosság maximalizálás (DRO).

## 📂 Könyvtárszerkezet

A projekt a `targets` pipeline kezelőt használja a reprodukálhatóság érdekében.

*   `R/`: A modellezési logika függvényei.
    *   `02_models.R`: BVAR és GARCH modellek.
    *   `03_strategies.R`: A Base Stratégiák implementációja (PortfolioAnalytics).
    *   `05_views.R`: Entropy Pooling (CVXR).
    *   `06_meta.R`: A Meta-optimalizáció.
*   `data/`: Bemeneti adatok (FRED makro + Yahoo ETF).
*   `_targets.R`: A pipeline vezérlőfájlja.
*   `config.yml`: Paraméterek (Horizont: 60 hó, Szimulációk: 10k).

## 🛠 Futtatás

A projekt tisztán R alapú (`renv` környezettel).

Indításhoz futtasd a `run.R` scriptet, vagy konzolban:

```r
targets::tar_make()