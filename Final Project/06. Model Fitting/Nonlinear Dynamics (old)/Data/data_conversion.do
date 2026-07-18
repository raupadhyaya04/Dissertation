clear all
set more off
cd "/Users/raupadhyaya04/Documents/GitHub/ECU33092/Final Project"


cd "/Users/raupadhyaya04/Documents/GitHub/ECU33092/Final Project/Data/Macro Data" // csv files got moved to a specific csv folder afterwards
import delimited "macro_fred_data_filled.csv", clear
save "macro_fred_data_filled.dta", replace

import delimited "macro_yahoo_data_filled.csv", clear
save "macro_yahoo_data_filled.dta", replace

cd "/Users/raupadhyaya04/Documents/GitHub/ECU33092/Final Project/Data/Crypto Data" // csv files got moved to a specific csv folder afterwards
import delimited "ADAUSDT.csv", clear
save "ADAUSDT.dta", replace

import delimited "AVAXUSDT.csv", clear
save "AVAXUSDT.dta", replace

import delimited "BNBUSDT.csv", clear
save "BNBUSDT.dta", replace

import delimited "DOGEUSDT.csv", clear
save "DOGEUSDT.dta", replace

import delimited "DOTUSDT.csv", clear
save "DOTUSDT.dta", replace

import delimited "ETHUSDT.csv", clear
save "ETHUSDT.dta", replace

import delimited "LINKUSDT.csv", clear
save "LINKUSDT.dta", replace

import delimited "SOLUSDT.csv", clear
save "SOLUSDT.dta", replace

import delimited "XRPUSDT.csv", clear
save "XRPUSDT.dta", replace

import delimited "BTCUSDT.csv", clear
save "BTCUSDT.dta", replace
