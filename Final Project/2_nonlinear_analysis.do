/*=============================================================================
  2_nonlinear_analysis.do
  -----------------------------------------------------------------------------
  PURPOSE : Nonlinear dynamics analysis of cryptocurrency returns against
            macro indicators (Gold Futures, USD Index, S&P 500, VIX) using:

    MODULE A – Descriptive & Preliminary
      A1. Summary statistics and distributional analysis (normality, fat tails)
      A2. BDS test for nonlinear dependence (Brock, Dechert & Scheinkman)
      A3. Autocorrelation / ARCH-LM tests (preconditions for nonlinearity)

    MODULE B – Regime Detection
      B1. Markov-Switching (MS) Model:
            – MS-AR(1) in returns  (two & three regimes)
            – MS regression: crypto return ~ macro regressors
          Prints: transition probabilities, expected durations, regime probs
      B2. Threshold / TAR model (threshold autoregression)
            – Threshold variable: lagged VIX (market stress regime)
            – Enders-Granger threshold test

    MODULE C – Volatility Dynamics
      C1. ARCH / GARCH(1,1) in returns  (benchmark linear volatility)
      C2. EGARCH(1,1)  – leverage / asymmetry in crypto volatility
      C3. Rolling 60-day betas vs macro indicators (time-varying sensitivity)

    MODULE D – Correlation Regime Analysis
      D1. DCC-GARCH concept approximation using rolling Pearson correlation
          between crypto returns and each macro variable (30 / 90 / 180-day)
      D2. Regime-conditioned correlations: compare correlations in
          "high-VIX" (VIX > 25) vs "low-VIX" regime

    MODULE E – Panel Nonlinear Regressions
      E1. Fixed-effects panel regression baseline
      E2. Panel Markov-switching: run coin-by-coin and collect regime params
      E3. Interaction terms: how macro effect differs by crypto regime state

    OUTPUTS:
      – Logs saved to  "Cleaned Data/analysis_log.smcl"
      – Graphs saved to "Cleaned Data/Graphs/"
      – Results tables saved to "Cleaned Data/Tables/"

  REQUIRES: "Cleaned Data/crypto_macro_daily.dta"  (built by 1_clean_data.do)
            Stata packages: mswitch, bdstest, arch, xtset (all base or ssc)
=============================================================================*/

clear all
set more off
set scheme s2color    // clean graphs

* ── Root & output paths ─────────────────────────────────────────────────────
local root   "/Users/raupadhyaya04/Documents/GitHub/ECU33092/Final Project"
local cdata  "`root'/Cleaned Data"
local output "`root'/Output"
local graphs "`output'/Graphs"
local tables "`output'/Tables"

* Create output folders if they do not exist
capture mkdir "`output'"
capture mkdir "`graphs'"
capture mkdir "`tables'"

* Start log
log using "`output'/analysis_log.smcl", replace smcl

* ── Load master panel ────────────────────────────────────────────────────────
use "`cdata'/crypto_macro_daily.dta", clear

* Re-declare panel structure
xtset coin_id date, delta(1)

* ── Create Sector Categorizations ────────────────────────────────────────────
gen sector = ""
replace sector = "Majors" if inlist(coin, "BTCUSDT", "ETHUSDT")
replace sector = "Alt_L1" if inlist(coin, "BNBUSDT", "SOLUSDT", "ADAUSDT", "AVAXUSDT", "DOTUSDT")
replace sector = "Infrastructure" if inlist(coin, "XRPUSDT", "LINKUSDT")
replace sector = "Meme" if coin == "DOGEUSDT"

* ── Convenience locals ───────────────────────────────────────────────────────
// Full coin list loaded for final categorical analysis
local coins    "BTCUSDT ETHUSDT BNBUSDT XRPUSDT SOLUSDT ADAUSDT DOTUSDT LINKUSDT AVAXUSDT DOGEUSDT"
// local coins    "BTCUSDT ETHUSDT" // (Uncomment for TEST MODE to save run-time)
local macrovars "usd_index gold_futures sp500 vix"
local macrorts  "r_usd r_gold r_sp500 r_vix"

********************************************************************************
* MODULE A – DESCRIPTIVE & PRELIMINARY DIAGNOSTICS
********************************************************************************
display as result _newline(2) "════════════════════════════════════════"
display as result " MODULE A: Descriptive & Preliminary"
display as result "════════════════════════════════════════"

*── A1. Summary statistics ──────────────────────────────────────────────────
estpost summarize ln_return abs_return usd_index gold_futures sp500 vix ///
    r_usd r_gold r_sp500 r_vix, detail
esttab using "`tables'/A1_summary_stats.csv", ///
    cells("mean sd min p25 p50 p75 max skewness kurtosis") ///
    label replace title("Summary Statistics – All Coins Pooled")

* Per-coin summary
foreach coin of local coins {
    display as result "  → Coin: `coin'"
    sum ln_return if coin == "`coin'", detail
}

*── A2. Distributional tests (per coin) ─────────────────────────────────────
* Jarque-Bera normality test & skewness/kurtosis test
tempname jb_results
file open `jb_results' using "`tables'/A2_normality_tests.csv", write replace
file write `jb_results' "coin,obs,skewness,kurtosis,jb_chi2,jb_p" _n

foreach coin of local coins {
    qui sum ln_return if coin == "`coin'", detail
    local sk  = r(skewness)
    local ku  = r(kurtosis)
    local n   = r(N)
    local jb  = `n'/6 * (`sk'^2 + (`ku'-3)^2/4)
    local pv  = chi2tail(2, `jb')
    file write `jb_results' "`coin',`n'," %6.4f (`sk') "," %6.4f (`ku') "," %8.2f (`jb') "," %6.4f (`pv') _n
    display "  `coin': Skew=`sk', Kurt=`ku', JB=`jb' (p=`pv')"
}
file close `jb_results'

* Distribution histograms – all coins on one page
graph bar (mean) abs_return, over(coin, label(angle(45))) ///
    ytitle("Mean |Return|") title("Mean Absolute Daily Returns by Coin") ///
    note("Proxy for realised volatility")
graph export "`graphs'/A2_mean_abs_returns.png", replace width(1200)

*── A3. ARCH-LM and autocorrelation tests ──────────────────────────────────
* Per coin: test for ARCH effects (precondition for nonlinear vol models)
tempname arch_results
file open `arch_results' using "`tables'/A3_arch_lm_tests.csv", write replace
file write `arch_results' "coin,lags,chi2,p_value,decision" _n

foreach coin of local coins {
    preserve
    keep if coin == "`coin'"
    qui tsset date
    * Fit a simple AR(1) to get residuals first
    qui regress ln_return l.ln_return
    predict double uhat_`coin' if e(sample), residuals
    * ARCH-LM test on residuals (5 lags)
    qui estat archlm, lags(5)
    matrix arch_mat = r(arch)
    matrix p_mat = r(p)
    local chi2_arch = arch_mat[1,1]
    local p_arch    = p_mat[1,1]
    local decision  = cond(`p_arch' < 0.05, "ARCH effects present", "No ARCH")
    file write `arch_results' "`coin',5," %8.4f (`chi2_arch') "," %6.4f (`p_arch') ",`decision'" _n
    display "  ARCH-LM `coin': chi2=`chi2_arch' p=`p_arch'  (`decision')"
    drop uhat_`coin'
    restore
}
file close `arch_results'

*── A4. BDS Test for Nonlinear Dependence ───────────────────────────────────
* BDS test: H0 = iid.  Rejection → nonlinear dependence in the series.
* Run on residuals from AR(1) to isolate nonlinear component.
* NOTE: requires 'bdstest' or the built-in 'wntestb' (different test).
* We use 'wntestb' (white noise test) as a base fallback, and bdstest if present.
display as result _newline "  A4. BDS tests (nonlinear dependence)"
tempname bds_results
file open `bds_results' using "`tables'/A4_bds_tests.csv", write replace
file write `bds_results' "coin,test,statistic,p_value,note" _n

* Check if bdstest is available
capture which bdstest
local bds_available = (_rc == 0)

foreach coin of local coins {
    preserve
    keep if coin == "`coin'"
    sort date
    tsset date

    qui regress ln_return l.ln_return
    predict double res_bds if e(sample), residuals

    if `bds_available' {
        * bdstest syntax: bdstest varname, reps(#)
        * embedding dimensions tested: 2 through 5 (default)
        capture bdstest res_bds, reps(200)
        if _rc == 0 {
            * r(bds) is a matrix: rows = dimension, cols = (z-stat, std err, p-val)
            matrix BDS = r(bds)
            local nr = rowsof(BDS)
            forvalues i = 1/`nr' {
                local z_val = BDS[`i', 1]
                local p_val = BDS[`i', 3]
                file write `bds_results' "`coin',BDS_m`i'," %8.4f (`z_val') "," %6.4f (`p_val') ",z-stat" _n
            }
            display "  BDS `coin': m=2 z=`=BDS[1,1]' p=`=BDS[1,3]'"
        }
        else {
            file write `bds_results' "`coin',BDS,NA,NA,bdstest error" _n
        }
    }
    else {
        * Fallback 1: Portmanteau white-noise test on squared residuals
        * (tests for ARCH-type nonlinear serial dependence)
        qui gen res_sq = res_bds^2
        capture wntestq res_sq, lags(10)
        if _rc == 0 {
            local chi2_wn = r(stat)
            local p_wn    = r(p)
            file write `bds_results' "`coin',WNtest_sq_resid," %8.4f (`chi2_wn') "," %6.4f (`p_wn') ",chi2(10) on sq.resid" _n
            display "  WN-test (sq resid) `coin': chi2=`chi2_wn' p=`p_wn'"
        }
        else {
            file write `bds_results' "`coin',BDS,NA,NA,bdstest not installed - run: net install bdstest from(http://fmwww.bc.edu/repec/bocode/b)" _n
            display as error "  Install bdstest: net install bdstest, from(http://fmwww.bc.edu/repec/bocode/b)"
        }
        capture drop res_sq
    }

    drop res_bds
    restore
}
file close `bds_results'

*── A4. Panel Unit Root Tests (Stationarity Proof) ──────────────────────────
display as result "  A4. Fisher-type Panel Unit Root Tests (Stationarity)"

tempname unitroot_file
file open `unitroot_file' using "`tables'/A4_unit_roots.csv", write replace
file write `unitroot_file' "Variable,Statistic_P,p_value,Stationary" _n

foreach var in ln_return r_gold r_usd r_sp500 r_vix {
    capture quietly xtunitroot fisher `var', dfuller lags(1)
    if _rc == 0 {
        local p_stat = r(P)
        local p_val  = r(p_P)
        local stat_txt = cond(`p_val' < 0.05, "Yes", "No")
        file write `unitroot_file' "`var'," %8.4f (`p_stat') "," %8.4f (`p_val') ",`stat_txt'" _n
    }
}
file close `unitroot_file'

********************************************************************************
* MODULE B – REGIME DETECTION
********************************************************************************
display as result _newline(2) "════════════════════════════════════════"
display as result " MODULE B: Markov-Switching & Threshold Models"
display as result "════════════════════════════════════════"

*── B1. Markov-Switching Models ─────────────────────────────────────────────
* Hamilton (1989) Markov-switching dynamic regression:
*   (i)  MS-DR 2-state with regime-switching variance (bear / bull)
*   (ii) MS-DR 3-state with regime-switching variance (bear / neutral / bull)
*   (iii) MS-DR regression with macro regressors – two regimes
*
* NOTE: mswitch ar with arswitch fails on crypto daily data due to numerical
*       derivative issues at AR(1) initialisation. We use mswitch dr instead
*       (regime-switching in the mean and variance; same economic content).
*
* Parameter names in Stata 18 mswitch dr output:
*   _b[State1:_cons]        = regime 1 mean
*   _b[State2:_cons]        = regime 2 mean
*   exp(_b[lnsigma1:_cons]) = regime 1 std dev  (lnsigma stored as log)
*   exp(_b[lnsigma2:_cons]) = regime 2 std dev
*   _b[p11:_cons]  stores logit(p12), so p11 = 1 - invlogit(_b[p11:_cons])
*   _b[p21:_cons]  stores logit(p22), so p22 = invlogit(_b[p21:_cons])
*   e(aic), e(ll) — stored; e(bic) not stored, compute manually

tempname ms_results
file open `ms_results' using "`tables'/B1_markov_switching_results.csv", write replace
file write `ms_results' "coin,model,states,mu_s1,mu_s2,mu_s3,sigma_s1,sigma_s2,p11,p22,AIC,BIC" _n

foreach coin of local coins {

    preserve
    keep if coin == "`coin'"
    * Drop first obs: no prior price on day 1 → missing ln_return, which
    * causes the EM algorithm's numerical derivative to fail
    drop if missing(ln_return)
    sort date
    tsset date

    display as result _newline "  ── Markov-Switching: `coin' ──"

    * ── (i) MS-DR 2-state, regime-switching mean + variance ───────────────
    capture noisily mswitch dr ln_return, states(2) varswitch nolog iter(20)

    if _rc == 0 {
        estat transition               // transition probability matrix
        estat duration                 // expected regime durations

        * Extract parameters (Stata 18 mswitch dr naming)
        local mu1  = _b[State1:_cons]
        local mu2  = _b[State2:_cons]
        local s1   = exp(_b[lnsigma1:_cons])
        local s2   = exp(_b[lnsigma2:_cons])
        * mswitch stores logit(p12) in p11:_cons → p11 = 1 - invlogit(raw)
        * mswitch stores logit(p22) in p21:_cons → p22 = invlogit(raw)
        local p11  = 1 - invlogit(_b[p11:_cons])
        local p22  = invlogit(_b[p21:_cons])
        local p21  = 1 - `p22'

        local aic  = e(aic)
        local ll   = e(ll)
        local n    = e(N)
        local bic  = (-2*`ll' + ln(`n')*6) / `n'   // 6 params; per-obs scale

        file write `ms_results' "`coin',MS-DR,2," %8.6f (`mu1') "," %8.6f (`mu2') ///
            ",NA," %8.6f (`s1') "," %8.6f (`s2') "," %6.4f (`p11') "," %6.4f (`p22') ///
            "," %8.2f (`aic') "," %8.2f (`bic') _n

        * Predicted smoothed regime probabilities
        predict double pr_s1_2st double pr_s2_2st if e(sample), pr

        * Plot regime probabilities over time
        twoway (area pr_s1_2st date, color(red%40) yaxis(1)) ///
               (line ln_return date, yaxis(2) lcolor(navy) lwidth(thin)), ///
            ytitle("Prob(Regime 1 – Bear)", axis(1)) ///
            ytitle("Log Return", axis(2)) ///
            title("MS 2-State Regime Probability: `coin'") ///
            xtitle("Date") legend(label(1 "P(Bear Regime)") label(2 "Log Return")) ///
            note("Red shading = Bear regime probability")
        graph export "`graphs'/B1_MS2state_`coin'.png", replace width(1400)

        drop pr_s1_2st pr_s2_2st
    }
    else {
        file write `ms_results' "`coin',MS-DR,2,NA,NA,NA,NA,NA,NA,NA,NA,NA" _n
        display as error "  mswitch failed for `coin' 2-state"
    }

    * ── (ii) MS-DR 3-state ─────────────────────────────────────────────
    capture noisily mswitch dr ln_return, states(3) varswitch nolog iter(20)

    if _rc == 0 {
        estat transition
        estat duration

        local mu1  = _b[State1:_cons]
        local mu2  = _b[State2:_cons]
        local mu3  = _b[State3:_cons]
        local s1   = exp(_b[lnsigma1:_cons])
        local s2   = exp(_b[lnsigma2:_cons])
        local p11  = 1 - invlogit(_b[p11:_cons])
        local p22  = invlogit(_b[p21:_cons])
        local ll   = e(ll)
        local n    = e(N)
        local aic  = e(aic)
        local bic  = (-2*`ll' + ln(`n')*12) / `n'  // 12 params for 3-state

        file write `ms_results' "`coin',MS-DR,3," %8.6f (`mu1') "," %8.6f (`mu2') ///
            "," %8.6f (`mu3') "," %8.6f (`s1') "," %8.6f (`s2') ///
            "," %6.4f (`p11') "," %6.4f (`p22') ///
            "," %8.2f (`aic') "," %8.2f (`bic') _n

        * Regime probability plot (smoothed probabilities)
        predict double pr_bear double pr_neut double pr_bull if e(sample), pr

        twoway (area pr_bear date, color(red%50))  ///
               (area pr_neut date, color(gray%40)) ///
               (area pr_bull date, color(green%40)), ///
            ytitle("Smoothed Regime Probability") xtitle("Date") ///
            title("MS 3-State Regime Probabilities: `coin'") ///
            legend(label(1 "Bear") label(2 "Neutral") label(3 "Bull"))
        graph export "`graphs'/B1_MS3state_`coin'.png", replace width(1400)

        drop pr_bear pr_neut pr_bull
    }
    else {
        file write `ms_results' "`coin',MS-DR,3,NA,NA,NA,NA,NA,NA,NA,NA,NA" _n
        display as error "  mswitch 3-state failed for `coin'"
    }

    * ── (iii) MS-DR Regression with Macro Regressors (2-state) ───────────
    * Stata 18 syntax: regressors go directly in varlist, varswitch for
    * regime-switching coefficients on those regressors
    capture noisily mswitch dr ln_return r_usd r_gold r_sp500 r_vix, ///
        states(2) varswitch nolog iter(20)

    if _rc == 0 {
        local aic = e(aic)
        local ll  = e(ll)
        local n   = e(N)
        local bic = (-2*`ll' + ln(`n')*12) / `n'   // 2 states × 6 params each
        file write `ms_results' "`coin',MS-Reg-Macro,2,NA,NA,NA,NA,NA,NA,NA," ///
            %8.2f (`aic') "," %8.2f (`bic') _n

        predict double pr_macro_s1 double dummy_pr_macro_s2 if e(sample), pr
        drop dummy_pr_macro_s2

        twoway (area pr_macro_s1 date, color(maroon%50)) ///
               (line ln_return date, lcolor(navy) lwidth(vthin) yaxis(2)), ///
            title("MS Regression (Macro Regressors) 2-State: `coin'") ///
            ytitle("P(Regime 1)", axis(1)) ytitle("Log Return", axis(2)) ///
            note("Regressors: ΔUSD, ΔGold, ΔS&P500, ΔVIX")
        graph export "`graphs'/B1_MS_macro_`coin'.png", replace width(1400)

        drop pr_macro_s1
    }
    else {
        file write `ms_results' "`coin',MS-Reg-Macro,2,NA,NA,NA,NA,NA,NA,NA,NA,NA" _n
        display as error "  mswitch macro-reg failed for `coin'"
    }

    restore
}
file close `ms_results'

*── B2. Threshold / TAR Models ──────────────────────────────────────────────
* Enders-Granger TAR: test whether the speed of adjustment to equilibrium
* differs across regimes defined by VIX threshold.
* We use a simple Self-Exciting TAR (SETAR) approximation:
*   ln_return = alpha_1*(I=1)*ln_return[t-1] + alpha_2*(I=0)*ln_return[t-1]
*   where I=1 if VIX[t-1] > threshold (default median VIX)

display as result _newline "  B2. Threshold / TAR models"

qui sum vix, detail
local vix_median = r(p50)
local vix_75     = r(p75)

gen byte high_vix = (vix > `vix_median') if !missing(vix)
gen byte crisis_vix = (vix > 25) if !missing(vix)      // common crisis threshold

* TAR interaction terms
gen ar_low_vix  = l.ln_return * (1 - high_vix)   // low-stress regime
gen ar_high_vix = l.ln_return * high_vix          // high-stress regime

tempname tar_results
file open `tar_results' using "`tables'/B2_tar_threshold_results.csv", write replace
file write `tar_results' "coin,alpha_low_vix,se_low,t_low,alpha_high_vix,se_high,t_high,F_linearity,p_linearity" _n

foreach coin of local coins {
    preserve
    keep if coin == "`coin'"
    sort date
    tsset date

    * TAR regression: two-regime AR(1) with VIX as threshold variable
    qui regress ln_return ar_low_vix ar_high_vix, noconstant
    local a_low  = _b[ar_low_vix]
    local se_low = _se[ar_low_vix]
    local t_low  = `a_low' / `se_low'
    local a_high = _b[ar_high_vix]
    local se_hi  = _se[ar_high_vix]
    local t_high = `a_high' / `se_hi'

    * Linearity test: H0 alpha_low = alpha_high (F-test)
    qui test ar_low_vix = ar_high_vix
    local F_lin = r(F)
    local p_lin = r(p)

    file write `tar_results' "`coin'," %8.6f (`a_low') "," %8.6f (`se_low') "," ///
        %6.3f (`t_low') "," %8.6f (`a_high') "," %8.6f (`se_hi') "," ///
        %6.3f (`t_high') "," %8.4f (`F_lin') "," %6.4f (`p_lin') _n

    display "  TAR `coin': alpha_low=`a_low' alpha_high=`a_high' F=`F_lin' p=`p_lin'"

    restore
}
file close `tar_results'

* Smooth Transition TAR (STAR) – logistic transition function
* LSTAR: transition function = 1 / (1 + exp(-gamma*(VIX - c)))
* This is estimated by NLS (nl) or approximated by cubic Taylor expansion

display as result "  B2b. LSTAR linearity test (Teräsvirta 3rd-order Taylor)"

tempname lstar_results
file open `lstar_results' using "`tables'/B2b_lstar_linearity.csv", write replace
file write `lstar_results' "coin,H01_p,H02_p,H03_p,overall_F_p,reject_linearity" _n

foreach coin of local coins {
    preserve
    keep if coin == "`coin'"
    sort date
    tsset date

    * Teräsvirta auxiliary test: regress ln_return on AR(1), VIX-AR(1) interactions
    * (approximation of LSTAR nonlinearity test)
    gen vix_lag1    = l.vix
    gen vix_lag1_sq = vix_lag1^2
    gen vix_lag1_cu = vix_lag1^3
    gen lret_lag1   = l.ln_return

    gen int1 = lret_lag1 * vix_lag1
    gen int2 = lret_lag1 * vix_lag1_sq
    gen int3 = lret_lag1 * vix_lag1_cu

    capture qui regress ln_return lret_lag1 int1 int2 int3
    if _rc == 0 {
        qui test int1
        local p1 = r(p)
        qui test int2
        local p2 = r(p)
        qui test int3
        local p3 = r(p)
        qui test int1 int2 int3
        local pF = r(p)
        local rej = cond(`pF' < 0.05, "YES", "NO")
        file write `lstar_results' "`coin'," %6.4f (`p1') "," %6.4f (`p2') "," %6.4f (`p3') "," %6.4f (`pF') ",`rej'" _n
        display "  LSTAR `coin': overall F p=`pF'  Reject linearity: `rej'"
    }
    restore
}
file close `lstar_results'

********************************************************************************
* MODULE C – VOLATILITY DYNAMICS
********************************************************************************
display as result _newline(2) "════════════════════════════════════════"
display as result " MODULE C: ARCH/GARCH Volatility Models"
display as result "════════════════════════════════════════"

*── C1 & C2. GARCH(1,1) and EGARCH(1,1) per coin ───────────────────────────
tempname garch_results
file open `garch_results' using "`tables'/C_garch_results.csv", write replace
file write `garch_results' "coin,model,omega,alpha,beta,alpha+beta,gamma_leverage,AIC,BIC" _n

foreach coin of local coins {
    preserve
    keep if coin == "`coin'"
    sort date
    tsset date

    display as result "  ── GARCH: `coin' ──"

    * ── GARCH(1,1) ──────────────────────────────────────────────────────
    capture {
        arch ln_return, arch(1) garch(1) distribution(t) nolog

        if _rc == 0 {
            local om = _b[ARCH:_cons]
            local a  = _b[ARCH:l.arch]
            local b  = _b[ARCH:l.garch]
            local ab = `a' + `b'
            local aic = -2*e(ll) + 2*e(k)
            local bic = -2*e(ll) + ln(e(N))*e(k)

            file write `garch_results' "`coin',GARCH11," %10.8f (`om') "," ///
                %8.6f (`a') "," %8.6f (`b') "," %8.6f (`ab') ",NA," ///
                %8.2f (`aic') "," %8.2f (`bic') _n

            * Conditional variance
            predict double cond_var_g_`coin' if e(sample), variance
            rename cond_var_g_`coin' condvar_garch

            twoway (line condvar_garch date, lcolor(navy) lwidth(thin)), ///
                title("GARCH(1,1) Conditional Variance: `coin'") ///
                ytitle("Conditional Variance") xtitle("Date")
            graph export "`graphs'/C1_GARCH_condvar_`coin'.png", replace width(1400)
            drop condvar_garch
        }
    }

    * ── EGARCH(1,1) – leverage effect ───────────────────────────────────
    capture {
        arch ln_return, earch(1) egarch(1) distribution(t) nolog

        if _rc == 0 {
            local om = _b[ARCH:_cons]
            local a  = _b[ARCH:l.earch]
            local b  = _b[ARCH:l.egarch]
            local gm = _b[ARCH:l.earch_a]   // asymmetry/leverage term
            local aic = -2*e(ll) + 2*e(k)
            local bic = -2*e(ll) + ln(e(N))*e(k)

            file write `garch_results' "`coin',EGARCH11," %10.8f (`om') "," ///
                %8.6f (`a') "," %8.6f (`b') ",NA," %8.6f (`gm') "," ///
                %8.2f (`aic') "," %8.2f (`bic') _n

            display "  EGARCH leverage (gamma): `gm'"
        }
    }

    restore
}
file close `garch_results'

*── C3. Rolling 60-day betas: crypto return on macro returns ─────────────────
display as result "  C3. Fast Rolling betas (60-day) using asreg"

* Run rolling regression for all coins at once using asreg
quietly: bys coin_id: asreg ln_return r_gold r_usd r_sp500 r_vix, window(date 60)

* Create time-series plots
foreach coin of local coins {
    capture {
        twoway (line _b_r_gold date if coin == "`coin'", lcolor(gold) lwidth(medthick)) ///
               (line _b_r_usd  date if coin == "`coin'", lcolor(blue) lwidth(medthick)) ///
               (line _b_r_sp500 date if coin == "`coin'", lcolor(green) lwidth(medthick)) ///
               (line _b_r_vix  date if coin == "`coin'", lcolor(red) lwidth(medthick)), ///
            yline(0, lcolor(black) lpattern(dash)) ///
            ytitle("Rolling 60-day Beta") xtitle("Date") ///
            title("Time-Varying Macro Sensitivities: `coin'") ///
            legend(label(1 "β Gold") label(2 "β USD") ///
                   label(3 "β S&P500") label(4 "β VIX")) ///
            note("OLS rolling window = 60 days")
        graph export "`graphs'/C3_rolling_betas_`coin'.png", replace width(1400)
    }
}

* Optionally save rolling betas to a CSV directly
preserve
keep date coin _b_r_gold _b_r_usd _b_r_sp500 _b_r_vix
rename (_b_r_gold _b_r_usd _b_r_sp500 _b_r_vix) (beta_r_gold beta_r_usd beta_r_sp500 beta_r_vix)
drop if missing(beta_r_gold)
export delimited using "`tables'/C3_rolling_betas.csv", replace
restore

capture drop _RMSE _R2 _Adj_R2 _b_r_gold _b_r_usd _b_r_sp500 _b_r_vix _b_cons _Nobs

********************************************************************************
* MODULE D – CORRELATION REGIME ANALYSIS
********************************************************************************
display as result _newline(2) "════════════════════════════════════════"
display as result " MODULE D: Regime-Conditioned Correlations"
display as result "════════════════════════════════════════"

*── D1. Rolling correlations (30 / 90 / 180-day) ────────────────────────────
display as result "  D1. Fast Rolling correlations (90-day) using rangestat"

* Calculate rolling 90-day correlations for all coins and macro variables
foreach mvar in r_gold r_usd r_sp500 r_vix {
    quietly: rangestat (corr) ln_return `mvar', interval(date -89 0) by(coin_id)
    rename corr_x corr_`mvar'
    drop corr_nobs
}

foreach coin of local coins {
    capture {
        twoway (line corr_r_gold  date if coin == "`coin'", lcolor(gold))  ///
               (line corr_r_usd   date if coin == "`coin'", lcolor(blue))  ///
               (line corr_r_sp500 date if coin == "`coin'", lcolor(green)) ///
               (line corr_r_vix   date if coin == "`coin'", lcolor(red)),   ///
            yline(0, lcolor(black) lpattern(dash)) ///
            yline(0.3 0.5, lcolor(gray) lpattern(shortdash)) ///
            ytitle("90-day Rolling Correlation") xtitle("Date") ///
            title("Rolling Correlations vs Macro: `coin'") ///
            legend(label(1 "ρ Gold") label(2 "ρ USD") ///
                   label(3 "ρ S&P500") label(4 "ρ VIX")) ///
            yscale(range(-1 1)) ylabel(-1(0.5)1) ///
            note("90-day rolling Pearson correlation")
        graph export "`graphs'/D1_rolling_corr_`coin'.png", replace width(1400)
    }
}

capture drop corr_r_gold corr_r_usd corr_r_sp500 corr_r_vix

*── D2. Regime-conditioned correlations: High-VIX vs Low-VIX ────────────────
display as result "  D2. Regime-conditioned correlations"

tempname regcorr_file
file open `regcorr_file' using "`tables'/D2_regime_correlations.csv", write replace
file write `regcorr_file' "coin,macro_var,regime,N,rho,p_value" _n

foreach coin of local coins {
    foreach mvar in r_gold r_usd r_sp500 r_vix {
        * Low VIX regime (VIX ≤ 25)
        qui correlate ln_return `mvar' if coin == "`coin'" & crisis_vix == 0
        local rho_low = r(rho)
        local n_low   = r(N)
        local t_low   = `rho_low' * sqrt(`n_low'-2) / sqrt(1 - `rho_low'^2)
        local p_low   = 2 * ttail(`n_low'-2, abs(`t_low'))
        file write `regcorr_file' "`coin',`mvar',low_vix,`n_low'," ///
            %8.6f (`rho_low') "," %6.4f (`p_low') _n

        * High VIX regime (VIX > 25)
        qui correlate ln_return `mvar' if coin == "`coin'" & crisis_vix == 1
        local rho_hi = r(rho)
        local n_hi   = r(N)
        local t_hi   = `rho_hi' * sqrt(`n_hi'-2) / sqrt(1 - `rho_hi'^2)
        local p_hi   = 2 * ttail(`n_hi'-2, abs(`t_hi'))
        file write `regcorr_file' "`coin',`mvar',high_vix,`n_hi'," ///
            %8.6f (`rho_hi') "," %6.4f (`p_hi') _n
    }
}
file close `regcorr_file'

* Summary bar chart: average correlation by regime
import delimited "`tables'/D2_regime_correlations.csv", clear
reshape wide rho p_value n, i(coin macro_var) j(regime) string
gen rho_diff = rholow_vix - rhohigh_vix

graph bar rholow_vix rhohigh_vix, over(macro_var) over(coin, label(angle(45))) ///
    legend(label(1 "Low VIX") label(2 "High VIX")) ///
    title("Correlation with Macro Indicators by VIX Regime") ///
    ytitle("Pearson ρ") yline(0, lcolor(black) lpattern(dash)) ///
    note("High VIX = VIX > 25 (crisis/stress regime)")
graph export "`graphs'/D2_regime_correlations_bar.png", replace width(1600)

use "`cdata'/crypto_macro_daily.dta", clear
xtset coin_id date, delta(1)

* ── Create Sector Categorizations ────────────────────────────────────────────
gen sector = ""
replace sector = "Majors" if inlist(coin, "BTCUSDT", "ETHUSDT")
replace sector = "Alt_L1" if inlist(coin, "BNBUSDT", "SOLUSDT", "ADAUSDT", "AVAXUSDT", "DOTUSDT")
replace sector = "Infrastructure" if inlist(coin, "XRPUSDT", "LINKUSDT")
replace sector = "Meme" if coin == "DOGEUSDT"

* Re-create variables that were lost when reloading the main dataset
qui sum vix, detail
local vix_median = r(p50)
gen byte high_vix = (vix > `vix_median') if !missing(vix)
gen byte crisis_vix = (vix > 25) if !missing(vix)

********************************************************************************
* MODULE E – PANEL NONLINEAR REGRESSIONS
********************************************************************************
display as result _newline(2) "════════════════════════════════════════"
display as result " MODULE E: Panel Nonlinear Regressions"
display as result "════════════════════════════════════════"

*── E1. Panel Fixed-Effects Baseline ─────────────────────────────────────────
display as result "  E1. Panel FE baseline"

* Hausman Test for Random vs Fixed Effects
display as result "  E1a. Hausman Test (FE vs RE)"
quietly xtreg ln_return r_gold r_usd r_sp500 r_vix treasury10y fed_funds l.ln_return, re
estimates store re_model
quietly xtreg ln_return r_gold r_usd r_sp500 r_vix treasury10y fed_funds l.ln_return, fe
estimates store fe_model
capture noisily hausman fe_model re_model

* Standard FE regression: crypto return ~ macro returns + controls
* Added Cross-Sectional Clustered SEs via vce(cluster coin_id)
xtreg ln_return r_gold r_usd r_sp500 r_vix ///
    treasury10y fed_funds l.ln_return, ///
    fe vce(cluster coin_id)
capture outreg2 using "`tables'/E1_panel_fe_baseline.doc", ///
    replace title("Panel FE: Crypto Returns on Macro (Clustered SE)") ///
    ctitle("All Coins")

* Alternative: estout / esttab
eststo m_fe_all: xtreg ln_return r_gold r_usd r_sp500 r_vix ///
    treasury10y fed_funds l.ln_return, fe vce(cluster coin_id)

* Run by sector
levelsof sector, local(sectors)
foreach sec of local sectors {
    capture eststo m_fe_`sec': xtreg ln_return r_gold r_usd r_sp500 r_vix ///
        treasury10y fed_funds l.ln_return if sector == "`sec'", fe vce(cluster coin_id)
}

esttab m_fe_all m_fe_* using "`tables'/E1_panel_fe_baseline.csv", ///
    b(4) se(4) star(* 0.10 ** 0.05 *** 0.01) ///
    title("Panel FE: Crypto Returns on Macro Variables (by Sector)") replace

*── E1b. Panel Fixed-Effects with Lagged Macro (Causality Check) ─────────────
display as result "  E1b. Lagged Macro Returns (Temporal Precedence)"

eststo m_fe_lag: xtreg ln_return l.r_gold l.r_usd l.r_sp500 l.r_vix ///
    l.treasury10y l.fed_funds l.ln_return, fe vce(cluster coin_id)

esttab m_fe_all m_fe_lag using "`tables'/E1b_causality_lags.csv", ///
    b(4) se(4) star(* 0.10 ** 0.05 *** 0.01) ///
    title("Panel FE: Contemporaneous vs Lagged Macro Impact") ///
    mtitles("Contemporaneous" "Lagged (t-1)") replace

*── E2. Nonlinear Panel: Regime × Macro interactions ─────────────────────────
display as result "  E2. Regime × Macro interaction panel"

* Interaction of high_vix regime with macro returns
gen vix_x_rgold  = high_vix * r_gold
gen vix_x_rusd   = high_vix * r_usd
gen vix_x_rsp500 = high_vix * r_sp500
gen vix_x_rvix   = high_vix * r_vix

label variable vix_x_rgold  "High VIX × ΔGold"
label variable vix_x_rusd   "High VIX × ΔUSD"
label variable vix_x_rsp500 "High VIX × ΔS&P500"
label variable vix_x_rvix   "High VIX × ΔVIX"

eststo m_interact: xtreg ln_return ///
    r_gold r_usd r_sp500 r_vix ///
    vix_x_rgold vix_x_rusd vix_x_rsp500 vix_x_rvix ///
    high_vix treasury10y fed_funds l.ln_return, ///
    fe vce(cluster coin_id)

esttab m_fe_all m_interact using "`tables'/E2_panel_interaction.csv", ///
    b(4) se(4) star(* 0.10 ** 0.05 *** 0.01) ///
    title("Panel FE: Nonlinear Regime Interactions") ///
    mtitles("Baseline FE" "Regime × Macro") replace

display as result "  Interaction model: test joint significance of regime effects"
testparm vix_x_rgold vix_x_rusd vix_x_rsp500 vix_x_rvix

*── E3. Coin-by-coin MS regression: collect regime-specific macro betas ──────
display as result "  E3. Coin-level MS-Regression summary"

tempname ms_macro_file
file open `ms_macro_file' using "`tables'/E3_ms_macro_betas.csv", write replace
file write `ms_macro_file' "coin,regime,beta_gold,beta_usd,beta_sp500,beta_vix,mu,sigma,AIC" _n

foreach coin of local coins {
    preserve
    keep if coin == "`coin'"
    drop if missing(ln_return)
    sort date
    tsset date

    * MS-DR with macro regressors and regime-switching coefficients
    capture noisily mswitch dr ln_return r_gold r_usd r_sp500 r_vix, ///
        states(2) varswitch nolog iter(20)

    if _rc == 0 {
        * Stata 18: regime-specific betas named State1:varname, State2:varname
        forvalues s = 1/2 {
            local bg  = _b[State`s':r_gold]
            local bu  = _b[State`s':r_usd]
            local bsp = _b[State`s':r_sp500]
            local bv  = _b[State`s':r_vix]
            local mu  = _b[State`s':_cons]
            local sig = exp(_b[lnsigma`s':_cons])
            local aic = e(aic)

            file write `ms_macro_file' "`coin',`s'," ///
                %8.6f (`bg') "," %8.6f (`bu') "," %8.6f (`bsp') "," ///
                %8.6f (`bv') "," %8.6f (`mu') "," %8.6f (`sig') "," ///
                %8.2f (`aic') _n
        }
    }
    else {
        file write `ms_macro_file' "`coin',NA,NA,NA,NA,NA,NA,NA,NA" _n
        display as error "  MS-macro failed for `coin'"
    }
    restore
}
file close `ms_macro_file'

********************************************************************************
* MODULE F – SUMMARY VISUALISATIONS
********************************************************************************
display as result _newline(2) "════════════════════════════════════════"
display as result " MODULE F: Summary Visualisations"
display as result "════════════════════════════════════════"

*── F1. Crypto log-price trajectories ──────────────────────────────────────
* Normalise all log prices to 0 at start for comparison
bysort coin_id (date): gen ln_price_norm = ln_price - ln_price[1]

twoway ///
    (line ln_price_norm date if coin == "BTCUSDT",  lcolor(orange))  ///
    (line ln_price_norm date if coin == "ETHUSDT",  lcolor(blue))    ///
    (line ln_price_norm date if coin == "BNBUSDT",  lcolor(green))   ///
    (line ln_price_norm date if coin == "SOLUSDT",  lcolor(purple))  ///
    (line ln_price_norm date if coin == "XRPUSDT",  lcolor(red))     ///
    (line ln_price_norm date if coin == "ADAUSDT",  lcolor(teal))    ///
    (line ln_price_norm date if coin == "DOGEUSDT", lcolor(brown))   ///
    (line ln_price_norm date if coin == "DOTUSDT",  lcolor(lime))    ///
    (line ln_price_norm date if coin == "LINKUSDT", lcolor(dknavy))  ///
    (line ln_price_norm date if coin == "AVAXUSDT", lcolor(cranberry)), ///
    title("Normalised Log-Price Trajectories: All Cryptocurrencies") ///
    ytitle("Log-Price (rebased to 0 at first observation)") ///
    xtitle("Date") ///
    legend(label(1 "BTC") label(2 "ETH") label(3 "BNB") label(4 "SOL") ///
           label(5 "XRP") label(6 "ADA") label(7 "DOGE") label(8 "DOT") ///
           label(9 "LINK") label(10 "AVAX") size(small)) ///
    note("Source: Binance hourly OHLCV data, collapsed to daily")
graph export "`graphs'/F1_logprice_trajectories.png", replace width(1600)

*── F2. Macro variables overview ──────────────────────────────────────────
* Use unique dates (one row per date)
preserve
keep date usd_index gold_futures sp500 vix treasury10y fed_funds oil_price
duplicates drop date, force
sort date

twoway (line gold_futures date, lcolor(gold) lwidth(medthick)), ///
    title("Gold Futures Price") ytitle("USD/oz") xtitle("Date")
graph export "`graphs'/F2a_gold_price.png", replace width(1200)

twoway (line usd_index date, lcolor(navy) lwidth(medthick)), ///
    title("USD Index (DXY)") ytitle("Index Level") xtitle("Date")
graph export "`graphs'/F2b_usd_index.png", replace width(1200)

twoway (line sp500 date, lcolor(green) lwidth(medthick)), ///
    title("S&P 500 Index Level") ytitle("Index Level") xtitle("Date")
graph export "`graphs'/F2c_sp500.png", replace width(1200)

twoway (line vix date, lcolor(red) lwidth(medthick)), ///
       yline(25, lcolor(gray) lpattern(dash) lwidth(thin)) ///
    title("VIX (CBOE Volatility Index)") ytitle("VIX Level") xtitle("Date") ///
    legend(off) note("Dashed line = VIX = 25 (crisis threshold)")
graph export "`graphs'/F2d_vix.png", replace width(1200)

restore

*── F3. Cross-coin return correlation heatmap ─────────────────────────────
* Pivot to wide for correlations
preserve
keep coin date ln_return
reshape wide ln_return, i(date) j(coin) string
corr ln_returnBTCUSDT ln_returnETHUSDT ln_returnBNBUSDT ln_returnSOLUSDT ///
     ln_returnXRPUSDT ln_returnADAUSDT ln_returnDOGEUSDT ln_returnDOTUSDT ///
     ln_returnLINKUSDT ln_returnAVAXUSDT
restore

display as result _newline "✓ All modules complete."
display as result "  Graphs → `graphs'/"
display as result "  Tables → `tables'/"

log close
