use "Cleaned Data/crypto_macro_daily.dta", clear
keep if coin == "BTCUSDT"
 tsset date
qui regress ln_return l.ln_return
predict double res_bds if e(sample), residuals
bdstest res_bds, reps(200)

