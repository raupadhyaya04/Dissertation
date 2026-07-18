# ECU33092 Final Project – Nonlinear Dynamics in Cryptocurrency Pricing

## Research Question

What macro indicators drive cryptocurrency prices, and does this relationship differ **non-linearly** across market regimes?

The analysis focuses on: **Gold Futures**, **USD Index (DXY)**, **S&P 500 Returns**, and **VIX** as macro covariates against the returns of 10 major cryptocurrencies (BTC, ETH, BNB, XRP, SOL, ADA, DOT, LINK, AVAX, DOGE).

---

## Project Structure

```
Final Project/
├── main.do                      ← Master script: runs the full pipeline
├── 1_clean_data.do              ← Data cleaning & merging (Step 1)
├── 2_nonlinear_analysis.do      ← All nonlinear models & outputs (Step 2)
├── Data/
│   ├── Crypto Data/csvs/        ← Raw hourly OHLCV data (Binance, 10 coins)
│   └── Macro Data/csvs/         ← Yahoo Finance + FRED macro data
└── Cleaned Data/                ← All outputs written here
    ├── crypto_macro_daily.dta   ← Master merged panel (coin × day)
    ├── Graphs/                  ← All .png charts
    └── Tables/                  ← All .csv results tables
```

---

## Data Sources

| Dataset                       | Source        | Frequency | Variables                                          |
| ----------------------------- | ------------- | --------- | -------------------------------------------------- |
| `macro_yahoo_data_filled.csv` | Yahoo Finance | Daily     | USD Index (DXY), Gold Futures (GC=F), S&P 500, VIX |
| `macro_fred_data_filled.csv`  | FRED          | Monthly   | Fed Funds Rate, CPI, M2, 10-yr Treasury, Crude Oil |
| `{COIN}USDT.csv`              | Binance API   | Hourly    | OHLCV (10 coins, 2017–2025)                        |

---

## Methodology & Modules

### Module A – Descriptive & Preliminary

- **Summary statistics** (mean, SD, skewness, excess kurtosis)
- **Jarque-Bera test** for normality of returns
- **ARCH-LM test** (Engle 1982) for conditional heteroscedasticity
- **BDS test** (Brock, Dechert & Scheinkman 1987/1996) for nonlinear serial dependence in AR(1) residuals

### Module B – Regime Detection ⭐

- **Markov-Switching AR(1) – 2 states** (Hamilton 1989): Bull vs Bear regimes
  - Regime means, variances, transition probabilities, expected durations
  - Smoothed regime probability plots
- **Markov-Switching AR(1) – 3 states**: Bull / Neutral / Bear
- **Markov-Switching Regression** with macro regressors (Gold, USD, S&P500, VIX) — regime-specific sensitivities
- **TAR (Threshold Autoregression)**: VIX as threshold variable, tests whether crypto return autocorrelation differs in high/low stress
- **LSTAR linearity test** (Teräsvirta 3rd-order Taylor expansion): formally tests whether a logistic smooth transition is needed

### Module C – Volatility Dynamics

- **GARCH(1,1)** with Student-t errors — persistence of crypto volatility
- **EGARCH(1,1)** — tests for leverage effects (asymmetric response to negative shocks)
- **Rolling 60-day OLS betas** vs each macro variable — time-varying macro sensitivities

### Module D – Correlation Regime Analysis

- **Rolling 90-day Pearson correlations** between crypto returns and macro returns
- **Regime-conditioned correlations**: compare ρ in Low-VIX (< 25) vs High-VIX (≥ 25) regimes

### Module E – Panel Nonlinear Regressions

- **Panel Fixed-Effects** baseline: pooled regression across all 10 coins
- **Regime × Macro interactions**: tests whether macro betas change in high-VIX regime
- **Coin-level MS-Regression**: regime-specific beta estimates extracted for each coin

### Module F – Visualisations

- Normalised log-price trajectories for all 10 coins
- Macro variable time-series charts
- Cross-coin return correlation matrix

---

## Required Stata Packages

Run once before executing the pipeline:

```stata
ssc install mswitch    // Markov-Switching (Kim & Nelson / Hamilton)
ssc install estout     // esttab / estout for regression tables
ssc install bdstest    // BDS nonlinearity test
ssc install outreg2    // Alternative table output (optional)
```

---

## How to Run

1. Open **Stata** (v15+ recommended for `mswitch`)
2. Run:
   ```stata
   do "/Users/raupadhyaya04/Documents/GitHub/ECU33092/Final Project/main.do"
   ```
3. All outputs are written to `Cleaned Data/`:
   - **`crypto_macro_daily.dta`** — the master analysis dataset
   - **`Graphs/`** — regime probability charts, rolling beta plots, correlation charts
   - **`Tables/`** — CSV results from all tests and models
   - **`analysis_log.smcl`** — full Stata log

---

## Key Nonlinear Dynamics Concepts Used

| Method              | What it identifies                                               |
| ------------------- | ---------------------------------------------------------------- |
| Markov-Switching    | Unobserved discrete regimes (structural breaks in mean/variance) |
| TAR / SETAR         | Asymmetric dynamics above/below a threshold                      |
| LSTAR               | Smooth transition between regimes (gradual regime change)        |
| EGARCH              | Asymmetric volatility (leverage effects)                         |
| BDS Test            | General nonlinear dependence beyond linear ARMA                  |
| Rolling Beta        | Time-varying (non-constant) macro sensitivities                  |
| Regime Correlations | Whether macro–crypto links change in crisis periods              |
