* run_B2_to_F.do  - Modules B2, C, D, E, F

local cdata "/Users/raupadhyaya04/Documents/GitHub/ECU33092/Final Project/Cleaned Data"
local tables "`cdata'/Tables"
local graphs  "`cdata'/Graphs"
local coins   "BTCUSDT ETHUSDT BNBUSDT XRPUSDT SOLUSDT ADAUSDT DOTUSDT LINKUSDT AVAXUSDT DOGEUSDT"

capture mkdir "`graphs'"
capture mkdir "`tables'"

log using "`cdata'/../run_B2_to_F.log", replace text

display as result "--- MODULES B2-F $(c(current_date)) ---"

* Load data & generate shared variables
use "`cdata'/crypto_macro_daily.dta", clear
xtset coin_id date

qui sum vix, detail
local vix_median = r(p50)
gen byte high_vix   = (vix > `vix_median') if !missing(vix)
gen byte crisis_vix = (vix > 25)           if !missing(vix)
gen ar_low_vix  = l.ln_return * (1 - high_vix)
gen ar_high_vix = l.ln_return * high_vix

* B2. TAR
display as result _n "--- B2. TAR ---"

file open ftar using "`tables'/B2_tar_threshold_results.csv", write replace
file write ftar "coin,alpha_low_vix,se_low,t_low,alpha_high_vix,se_high,t_high,F_linearity,p_linearity" _n

foreach coin of local coins {
    preserve
    keep if coin == "`coin'"
    sort date
    tsset date
    qui regress ln_return ar_low_vix ar_high_vix, noconstant
    local al = _b[ar_low_vix]
    local sl = _se[ar_low_vix]
    local tl = `al'/`sl'
    local ah = _b[ar_high_vix]
    local sh = _se[ar_high_vix]
    local th = `ah'/`sh'
    qui test ar_low_vix = ar_high_vix
    local Fv = r(F)
    local pv = r(p)
    file write ftar "`coin'," %8.6f (`al') "," %8.6f (`sl') "," %6.3f (`tl') "," %8.6f (`ah') "," %8.6f (`sh') "," %6.3f (`th') "," %8.4f (`Fv') "," %6.4f (`pv') _n
    display "  `coin': al=`al'  ah=`ah'  F_p=`pv'"
    restore
}
file close ftar
display as result "  TAR done"

* B2b. LSTAR linearity test
display as result _n "--- B2b. LSTAR ---"

file open flstar using "`tables'/B2b_lstar_linearity.csv", write replace
file write flstar "coin,H01_p,H02_p,H03_p,overall_F_p,reject_linearity" _n

foreach coin of local coins {
    preserve
    keep if coin == "`coin'"
    sort date
    tsset date
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
        file write flstar "`coin'," %6.4f (`p1') "," %6.4f (`p2') "," %6.4f (`p3') "," %6.4f (`pF') ",`rej'" _n
        display "  `coin': F_p=`pF'  Reject: `rej'"
    }
    else {
        file write flstar "`coin',NA,NA,NA,NA,NA" _n
        display as error "  `coin': LSTAR failed"
    }
    restore
}
file close flstar
display as result "  LSTAR done"

* C. GARCH(1,1) and EGARCH(1,1)
display as result _n "--- C. GARCH/EGARCH ---"

file open fgarch using "`tables'/C_garch_results.csv", write replace
file write fgarch "coin,model,omega,alpha,beta,alpha_plus_beta,gamma_leverage,AIC,BIC" _n

foreach coin of local coins {
    preserve
    keep if coin == "`coin'"
    sort date
    tsset date
    display "  GARCH: `coin'"

    capture noisily arch ln_return, arch(1) garch(1) distribution(t) nolog
    if _rc == 0 {
        local om = _b[ARCH:_cons]
        local al = _b[ARCH:l.arch]
        local be = _b[ARCH:l.garch]
        local ab = `al' + `be'
        local aic = e(aic)
        local bic = e(bic)
        file write fgarch "`coin',GARCH11," %10.8f (`om') "," %8.6f (`al') "," %8.6f (`be') "," %8.6f (`ab') ",NA," %8.2f (`aic') "," %8.2f (`bic') _n
        display "  GARCH OK  a+b=`ab'"
    }
    else {
        file write fgarch "`coin',GARCH11,NA,NA,NA,NA,NA,NA,NA" _n
        display as error "  GARCH(1,1) failed: `coin'"
    }

    capture noisily arch ln_return, earch(1) egarch(1) distribution(t) nolog
    if _rc == 0 {
        local om = _b[ARCH:_cons]
        local al = _b[ARCH:l.earch]
        local be = _b[ARCH:l.egarch]
        capture local gam = _b[ARCH:l.earch_a]
        if _rc != 0 local gam = .
        local aic = e(aic)
        local bic = e(bic)
        file write fgarch "`coin',EGARCH11," %10.8f (`om') "," %8.6f (`al') "," %8.6f (`be') ",NA," %8.6f (`gam') "," %8.2f (`aic') "," %8.2f (`bic') _n
        display "  EGARCH OK  gam=`gam'"
    }
    else {
        file write fgarch "`coin',EGARCH11,NA,NA,NA,NA,NA,NA,NA" _n
        display as error "  EGARCH(1,1) failed: `coin'"
    }
    restore
}
file close fgarch
display as result "  GARCH done"

* C3. Rolling 60-day betas vs S&P500
display as result _n "--- C3. Rolling betas ---"
capture mkdir "`cdata'/rolling_betas"

foreach coin of local coins {
    preserve
    keep if coin == "`coin'"
    sort date
    tsset date

    capture rolling beta_sp500=_b[r_sp500], window(60) saving("`cdata'/rolling_betas/beta_`coin'.dta", replace) nodots: regress ln_return r_sp500

    if _rc == 0 {
        use "`cdata'/rolling_betas/beta_`coin'.dta", clear
        rename end date
        format date %td
        twoway (line beta_sp500 date, lcolor(navy) lwidth(thin)), yline(0, lcolor(red) lpattern(dash)) title("60-day Rolling Beta vs S&P500: `coin'") ytitle("Beta S&P500") xtitle("Date")
        graph export "`graphs'/C3_rolling_beta_sp500_`coin'.png", replace width(1200)
        display "  `coin': rolling beta saved"
    }
    else {
        display as error "  `coin': rolling beta failed"
    }
    restore
}
display as result "  C3 done"

* D1. Rolling 90-day correlations
display as result _n "--- D1. Rolling correlations ---"
capture mkdir "`cdata'/corr_data"

foreach coin of local coins {
    preserve
    keep if coin == "`coin'"
    sort date
    tsset date

    foreach mvar in r_gold r_usd r_sp500 r_vix {
        capture rolling corr_`mvar'=r(rho), window(90) saving("`cdata'/corr_data/tmp_`coin'_`mvar'.dta", replace) nodots: correlate ln_return `mvar'
    }

    capture {
        use "`cdata'/corr_data/tmp_`coin'_r_gold.dta", clear
        foreach mvar in r_usd r_sp500 r_vix {
            merge 1:1 end using "`cdata'/corr_data/tmp_`coin'_`mvar'.dta", nogenerate
        }
        rename end date
        format date %td
        twoway (line corr_r_gold date, lcolor(gold)) (line corr_r_usd date, lcolor(blue)) (line corr_r_sp500 date, lcolor(green)) (line corr_r_vix date, lcolor(red)), yline(0, lcolor(black) lpattern(dash)) ytitle("90-day Rolling Corr") xtitle("Date") title("Rolling Corr vs Macro: `coin'") legend(label(1 "Gold") label(2 "USD") label(3 "SP500") label(4 "VIX")) yscale(range(-1 1)) ylabel(-1(0.5)1)
        graph export "`graphs'/D1_rolling_corr_`coin'.png", replace width(1400)
        foreach mvar in r_gold r_usd r_sp500 r_vix {
            capture erase "`cdata'/corr_data/tmp_`coin'_`mvar'.dta"
        }
        display "  `coin': rolling corr saved"
    }
    restore
}
display as result "  D1 done"

* D2. Regime-conditioned correlations
display as result _n "--- D2. Regime correlations ---"

use "`cdata'/crypto_macro_daily.dta", clear
xtset coin_id date
qui sum vix, detail
gen byte crisis_vix = (vix > 25) if !missing(vix)

file open fd2 using "`tables'/D2_regime_correlations.csv", write replace
file write fd2 "coin,macro_var,regime,N,rho,p_value" _n

foreach coin of local coins {
    foreach mvar in r_gold r_usd r_sp500 r_vix {
        qui correlate ln_return `mvar' if coin == "`coin'" & crisis_vix == 0
        local rlo = r(rho)
        local nlo = r(N)
        if `nlo' > 2 {
            local plo = 2 * ttail(`nlo'-2, abs(`rlo'*sqrt(`nlo'-2)/sqrt(1-`rlo'^2)))
        }
        else {
            local plo = .
        }
        file write fd2 "`coin',`mvar',low_vix,`nlo'," %8.6f (`rlo') "," %6.4f (`plo') _n

        qui correlate ln_return `mvar' if coin == "`coin'" & crisis_vix == 1
        local rhi = r(rho)
        local nhi = r(N)
        if `nhi' > 2 {
            local phi = 2 * ttail(`nhi'-2, abs(`rhi'*sqrt(`nhi'-2)/sqrt(1-`rhi'^2)))
        }
        else {
            local phi = .
        }
        file write fd2 "`coin',`mvar',high_vix,`nhi'," %8.6f (`rhi') "," %6.4f (`phi') _n
    }
}
file close fd2
display as result "  D2 done"

* E. Panel Fixed Effects
display as result _n "--- E. Panel FE ---"

use "`cdata'/crypto_macro_daily.dta", clear
xtset coin_id date
qui sum vix, detail
gen byte crisis_vix = (vix > 25) if !missing(vix)

xtreg ln_return r_gold r_usd r_sp500 r_vix, fe vce(robust)
estimates store FE_baseline

xtreg ln_return r_gold r_usd r_sp500 r_vix c.r_vix#i.crisis_vix, fe vce(robust)
estimates store FE_crisis

file open ffe using "`tables'/E_panel_results.csv", write replace
file write ffe "model,beta_gold,beta_usd,beta_sp500,beta_vix,N,R2_within" _n

estimates restore FE_baseline
local bg  = _b[r_gold]
local bu  = _b[r_usd]
local bsp = _b[r_sp500]
local bv  = _b[r_vix]
local eN  = e(N)
local r2  = e(r2_w)
file write ffe "FE_baseline," %8.6f (`bg') "," %8.6f (`bu') "," %8.6f (`bsp') "," %8.6f (`bv') "," %g (`eN') "," %6.4f (`r2') _n

estimates restore FE_crisis
local bg  = _b[r_gold]
local bu  = _b[r_usd]
local bsp = _b[r_sp500]
local bv  = _b[r_vix]
local eN  = e(N)
local r2  = e(r2_w)
file write ffe "FE_crisis_int," %8.6f (`bg') "," %8.6f (`bu') "," %8.6f (`bsp') "," %8.6f (`bv') "," %g (`eN') "," %6.4f (`r2') _n

file close ffe
display as result "  E done"

* F. Summary visualisations
display as result _n "--- F. Visualisations ---"

use "`cdata'/crypto_macro_daily.dta", clear
xtset coin_id date
bysort coin_id (date): gen ln_price_norm = ln_price - ln_price[1]

twoway (line ln_price_norm date if coin=="BTCUSDT", lcolor(orange)) (line ln_price_norm date if coin=="ETHUSDT", lcolor(blue)) (line ln_price_norm date if coin=="BNBUSDT", lcolor(green)) (line ln_price_norm date if coin=="SOLUSDT", lcolor(purple)) (line ln_price_norm date if coin=="XRPUSDT", lcolor(red)) (line ln_price_norm date if coin=="ADAUSDT", lcolor(teal)) (line ln_price_norm date if coin=="DOGEUSDT", lcolor(brown)) (line ln_price_norm date if coin=="DOTUSDT", lcolor(lime)) (line ln_price_norm date if coin=="LINKUSDT", lcolor(dknavy)) (line ln_price_norm date if coin=="AVAXUSDT", lcolor(sienna)), title("Normalised Log-Price Trajectories") ytitle("Cumulative Log Return") xtitle("Date") legend(label(1 "BTC") label(2 "ETH") label(3 "BNB") label(4 "SOL") label(5 "XRP") label(6 "ADA") label(7 "DOGE") label(8 "DOT") label(9 "LINK") label(10 "AVAX") rows(2)) note("Source: Binance daily OHLCV data")
graph export "`graphs'/F1_logprice_trajectories.png", replace width(1600)
display "  F1 saved"

preserve
keep date usd_index gold_futures sp500 vix
duplicates drop date, force
sort date

twoway (line gold_futures date, lcolor(gold) lwidth(medthick)), title("Gold Futures") ytitle("USD/oz") xtitle("Date")
graph export "`graphs'/F2a_gold_price.png", replace width(1200)

twoway (line usd_index date, lcolor(navy) lwidth(medthick)), title("USD Index (DXY)") ytitle("Index") xtitle("Date")
graph export "`graphs'/F2b_usd_index.png", replace width(1200)

twoway (line sp500 date, lcolor(green) lwidth(medthick)), title("S&P 500") ytitle("Index") xtitle("Date")
graph export "`graphs'/F2c_sp500.png", replace width(1200)

twoway (line vix date, lcolor(red) lwidth(medthick)), yline(25, lcolor(gray) lpattern(dash)) title("VIX") ytitle("VIX") xtitle("Date") legend(off) note("Dashed = VIX 25")
graph export "`graphs'/F2d_vix.png", replace width(1200)
restore
display "  F2a-d saved"

display as result _n "--- ALL MODULES B2-F COMPLETE ---"
display as result "  Tables -> `tables'/"
display as result "  Graphs -> `graphs'/"

log close
