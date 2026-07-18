/*=============================================================================
  1_clean_data.do
  -----------------------------------------------------------------------------
  PURPOSE : Build a clean, balanced daily panel dataset that merges:
              (a) 10 cryptocurrency closing prices  (hourly → daily)
              (b) Yahoo Finance macro series        (daily: USD Index, Gold
                  Futures, S&P 500, VIX)
              (c) FRED macro series                 (monthly → daily via
                  forward-fill: Fed Funds, CPI, M2, 10-yr Treasury, Oil)

  OUTPUT  : "Cleaned Data/crypto_macro_daily.dta"
              – long-format panel (coin × date) with log-prices, log-returns,
                and all macro covariates.

  REQUIRES: Data already saved as .dta in their respective sub-folders by
            data_conversion.do  (or the csvs/ folders are used directly here).
=============================================================================*/

clear all
set more off

* ── Root path ──────────────────────────────────────────────────────────────
local root "/Users/raupadhyaya04/Documents/GitHub/ECU33092/Final Project"
cd "`root'"

********************************************************************************
* SECTION 1 – CRYPTO DATA: hourly → daily (close price & volume)
********************************************************************************

* List of coins
local coins BTCUSDT ETHUSDT BNBUSDT XRPUSDT SOLUSDT ADAUSDT DOTUSDT LINKUSDT AVAXUSDT DOGEUSDT

* We will save a daily .dta for each coin, then append all into a panel
local first_coin = 1

foreach coin of local coins {

    import delimited "`root'/Data/Crypto Data/csvs/`coin'.csv", ///
        clear varnames(1) stringcols(1)

    * ── Date parsing ─────────────────────────────────────────────────────
    * timestamp format: "2017-08-17 04:00:00"
    gen double date_num = clock(timestamp, "YMDhms")
    format date_num %tc
    gen date = dofc(date_num)          // convert to daily date
    format date %td

    * ── Keep only needed columns (keep date_num for intra-day sorting) ───
    keep date date_num close volume
    destring close volume, replace force

    * ── Collapse hourly → daily (OHLC close = last price of the day) ──────
    * Use last observation within each day as the daily closing price,
    * and sum volume for the day.
    sort date date_num     // sort within day by time
    bysort date (date_num): gen _last  = (_n == _N)
    bysort date (date_num): gen _first = (_n == 1)
    gen  daily_close  = close  if _last  == 1
    gen  daily_volume = volume
    gen  daily_open   = close  if _first == 1

    collapse (lastnm) daily_close (sum) daily_volume ///
             (firstnm) daily_open, by(date)

    rename daily_close  price_close
    rename daily_volume volume
    rename daily_open   price_open

    * ── Returns ──────────────────────────────────────────────────────────
    sort date
    gen ln_price  = ln(price_close)
    gen ln_return = ln_price - ln_price[_n-1]

    * ── Volatility proxy: |return| and squared return ────────────────────
    gen abs_return = abs(ln_return)
    gen sq_return  = ln_return^2

    * ── Coin identifier ──────────────────────────────────────────────────
    gen coin = "`coin'"

    * Save temporary file
    save "`root'/Cleaned Data/tmp_`coin'.dta", replace

    display as result "✓ `coin' processed"
}

* ── Append all coins into a long panel ──────────────────────────────────────
use "`root'/Cleaned Data/tmp_BTCUSDT.dta", clear
foreach coin in ETHUSDT BNBUSDT XRPUSDT SOLUSDT ADAUSDT DOTUSDT LINKUSDT AVAXUSDT DOGEUSDT {
    append using "`root'/Cleaned Data/tmp_`coin'.dta"
}

* ── Panel ID ──────────────────────────────────────────────────────────────
encode coin, gen(coin_id)
label variable coin_id "Coin identifier (numeric)"
order coin coin_id date

sort coin_id date

* Save intermediate crypto panel
save "`root'/Cleaned Data/crypto_daily_panel.dta", replace
display as result "✓ Crypto daily panel saved"

* ── Clean up temp files ───────────────────────────────────────────────────
foreach coin of local coins {
    erase "`root'/Cleaned Data/tmp_`coin'.dta"
}

********************************************************************************
* SECTION 2 – YAHOO MACRO DATA: daily (USD Index, Gold Futures, S&P500, VIX)
********************************************************************************
* CSV header: Date,DX-Y.NYB,GC=F,^GSPC,^VIX
* Special characters in header names cause Stata to mangle them unpredictably.
* Use varnames(nonames) and assign names by column position — always reliable.

import delimited "`root'/Data/Macro Data/csvs/macro_yahoo_data_filled.csv", ///
    clear colrange(1:5) varnames(nonames)

rename (v1 v2 v3 v4 v5) (date_str usd_index gold_futures sp500 vix)
drop if date_str == "Date"    // remove header row

* Parse date
gen date = date(date_str, "YMD")
format date %td
drop date_str

destring usd_index gold_futures sp500 vix, replace force

* ── Macro derived variables ───────────────────────────────────────────────
sort date

* Log-levels
gen ln_gold    = ln(gold_futures)
gen ln_usd     = ln(usd_index)
gen ln_sp500   = ln(sp500)
gen ln_vix     = ln(vix)

* Daily log-returns
gen r_gold  = ln(gold_futures) - ln(gold_futures[_n-1])
gen r_usd   = ln(usd_index)    - ln(usd_index[_n-1])
gen r_sp500 = ln(sp500)        - ln(sp500[_n-1])
gen r_vix   = ln(vix)          - ln(vix[_n-1])

* VIX level as fear gauge (keep raw)
label variable vix          "CBOE VIX (implied vol)"
label variable usd_index    "USD Index (DXY)"
label variable gold_futures "Gold Futures price (USD/oz)"
label variable sp500        "S&P 500 Index level"
label variable r_gold       "Daily log-return: Gold"
label variable r_usd        "Daily log-return: USD Index"
label variable r_sp500      "Daily log-return: S&P 500"
label variable r_vix        "Daily log-return: VIX"

save "`root'/Cleaned Data/yahoo_macro_daily.dta", replace
display as result "✓ Yahoo macro data saved"

********************************************************************************
* SECTION 3 – FRED MACRO DATA: monthly → daily (forward-fill)
********************************************************************************
* CSV header row: ,Fed Funds Rate,Consumer Price Index,M2 Money Supply,...
* Dates are stored as the first day of each month: "2016-01-01" (YMD format)
* Use varnames(nonames) to avoid header mangling, then drop header row.

import delimited "`root'/Data/Macro Data/csvs/macro_fred_data_filled.csv", ///
    clear colrange(1:6) varnames(nonames)

rename (v1 v2 v3 v4 v5 v6) (date_str fed_funds cpi m2 treasury10y oil_price)

* Drop header row (v1 will contain the literal string from the CSV header)
drop if missing(date_str) | date_str == "" | length(date_str) < 6
* Header row has non-date text in date_str column
drop if regexm(date_str, "[A-Za-z]")   // drop any row where date_str has letters

destring fed_funds cpi m2 treasury10y oil_price, replace force

* Parse full date "2016-01-01" using YMD
gen date = date(date_str, "YMD")
format date %td
drop date_str

* Verify: should be ~120 monthly obs
display "FRED obs after cleaning: " _N

* ── Expand monthly → daily using forward-fill ─────────────────────────────
* Save monthly anchor dates
sum date
local dmin = r(min)
local dmax = r(max)
local dmax_ext = `dmax' + 31   // cover to end of last month

preserve
    keep date fed_funds cpi m2 treasury10y oil_price
    rename date date_monthly_start
    save "`root'/Cleaned Data/tmp_fred_monthly.dta", replace
restore

* Build daily spine
clear
local ndays = `dmax_ext' - `dmin' + 1
set obs `ndays'
gen date = `dmin' + _n - 1
format date %td

* Merge with monthly data on the first-of-month observations
gen date_monthly_start = date   // will only match at month boundaries

merge m:1 date_monthly_start using "`root'/Cleaned Data/tmp_fred_monthly.dta", ///
    nogenerate keep(master match)

* Forward-fill (locf) all FRED variables
sort date
foreach v in fed_funds cpi m2 treasury10y oil_price {
    replace `v' = `v'[_n-1] if missing(`v') & _n > 1
}

* Drop spine column
drop date_monthly_start

* Label
label variable fed_funds    "Federal Funds Rate (%)"
label variable cpi           "Consumer Price Index"
label variable m2            "M2 Money Supply (USD bn)"
label variable treasury10y   "10-Year Treasury Yield (%)"
label variable oil_price     "Crude Oil Price (USD/bbl)"

save "`root'/Cleaned Data/fred_macro_daily.dta", replace
erase "`root'/Cleaned Data/tmp_fred_monthly.dta"
display as result "✓ FRED macro data (daily interpolated) saved"

********************************************************************************
* SECTION 4 – MERGE: Crypto Panel + Yahoo + FRED
********************************************************************************

use "`root'/Cleaned Data/crypto_daily_panel.dta", clear

* ── Merge Yahoo (daily) macro ─────────────────────────────────────────────
merge m:1 date using "`root'/Cleaned Data/yahoo_macro_daily.dta", ///
    keep(master match) nogenerate

* ── Merge FRED (daily interpolated) macro ────────────────────────────────
merge m:1 date using "`root'/Cleaned Data/fred_macro_daily.dta", ///
    keep(master match) nogenerate

* ── Additional derived variables ─────────────────────────────────────────

* Real return: crypto return minus inflation proxy (monthly CPI change)
* (CPI change is slow-moving; use as deflator concept)
bysort coin_id (date): gen d_ln_cpi = ln(cpi) - ln(cpi[_n-1])
gen real_ln_return = ln_return - d_ln_cpi

* Relative volatility: crypto abs-return normalised by VIX
gen rel_vol = abs_return / (vix / 100)

* Fear/Greed proxy: VIX * USD Index interaction term
gen fear_index = vix * usd_index

* Gold-to-VIX ratio (safe-haven demand indicator)
gen gold_vix_ratio = gold_futures / vix

* ── Drop pre-2017 observations (most coins only start 2017-2018) ─────────
drop if date < td(01jan2017)

* ── Drop weekends for macro data completeness (optional: keep all) ────────
* We keep all days as crypto trades 24/7; macro vars forward-filled over weekends

* ── Handle remaining missings ────────────────────────────────────────────
* Drop observations where core crypto price is missing
drop if missing(price_close)

* ── Variable labels ──────────────────────────────────────────────────────
label variable date          "Calendar date"
label variable coin          "Cryptocurrency ticker"
label variable coin_id       "Coin identifier (numeric)"
label variable price_close   "Daily closing price (USDT)"
label variable price_open    "Daily opening price (USDT)"
label variable volume        "Daily trading volume (base coin)"
label variable ln_price      "Log daily closing price"
label variable ln_return     "Log daily return"
label variable abs_return    "Absolute log daily return (realized vol proxy)"
label variable sq_return     "Squared log daily return"
label variable real_ln_return "Real log return (CPI-adjusted)"
label variable rel_vol       "Crypto abs-return / (VIX/100)"
label variable fear_index    "VIX × USD Index interaction"
label variable gold_vix_ratio "Gold Futures / VIX ratio"

* ── Set panel structure ───────────────────────────────────────────────────
sort coin_id date
xtset coin_id date, delta(1)

* ── Summary of final panel ───────────────────────────────────────────────
display as result _newline "=== FINAL PANEL SUMMARY ==="
xtdescribe
summarize price_close ln_return usd_index gold_futures sp500 vix ///
          fed_funds treasury10y oil_price cpi m2

* ── Save ─────────────────────────────────────────────────────────────────
save "`root'/Cleaned Data/crypto_macro_daily.dta", replace

* Count variables via ds
quietly ds
local nvars : word count `r(varlist)'
display as result _newline "✓ Master panel saved: Cleaned Data/crypto_macro_daily.dta"
display as result "  Observations: " _N
display as result "  Variables:    " `nvars'

* ── Export a summary CSV for quick inspection ────────────────────────────
export delimited "`root'/Cleaned Data/crypto_macro_daily_preview.csv", replace
