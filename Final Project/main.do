/*=============================================================================
  main.do  –  ECU33092 Final Project
  -----------------------------------------------------------------------------
  TOPIC   : Nonlinear Dynamics in Cryptocurrency Pricing vs Macro Indicators
  AUTHOR  : [Your Name]
  DATE    : March 2026

  PIPELINE:
    1_clean_data.do          → builds Cleaned Data/crypto_macro_daily.dta
    2_nonlinear_analysis.do  → runs all nonlinear dynamics models & produces
                               graphs (Cleaned Data/Graphs/) and
                               tables (Cleaned Data/Tables/)

  REQUIRED SSC PACKAGES (install once before running):
    ssc install mswitch     // Markov-switching models
    ssc install bdstest     // BDS nonlinearity test
    ssc install esttab      // Regression output tables  (part of estout)
    ssc install estout      // estout / esttab
    ssc install outreg2     // Alternative table output (optional)
    net install st0239.pkg  // rolling (included in Stata 12+, usually built-in)

  HOW TO RUN:
    Open Stata, set working directory to the Final Project folder, then:
    do main.do
=============================================================================*/

clear all
set more off

* ── Global root path (edit if moving the project) ────────────────────────────
global root "/Users/raupadhyaya04/Documents/GitHub/ECU33092/Final Project"
cd "$root"

* ── Create output directories ───────────────────────────────────────────────
capture mkdir "$root/Cleaned Data"
capture mkdir "$root/Cleaned Data/Graphs"
capture mkdir "$root/Cleaned Data/Tables"

* ── Install required packages (only runs if not already installed) ────────────
* mswitch: built into Stata 14+ (no install needed); for older versions use SSC
capture which mswitch
if _rc != 0 {
    capture ssc install mswitch, replace
    if _rc != 0 {
        display as error "mswitch not found. Requires Stata 14+ (built-in) or:"
        display as error "  net from http://www.stata.com/users/bbdg"
        display as error "Continuing without mswitch – Module B will be skipped."
    }
}

* estout / esttab
capture which esttab
if _rc != 0 ssc install estout, replace

* bdstest: hosted on SSC as part of 'sts' suite OR via net install
* NOT available as standalone 'bdstest' on SSC — use the correct package name
capture which bdstest
if _rc != 0 {
    * Try the correct SSC package name: 'sts' contains bdstest in some Stata versions
    capture ssc install sts, replace
    if _rc != 0 {
        * Fallback: install directly from Baum's BC archive
        capture net install bdstest, ///
            from("http://fmwww.bc.edu/repec/bocode/b") replace
        if _rc != 0 {
            display as result "bdstest not installed – BDS test (Module A4) will be skipped."
            display as result "To install manually: net install bdstest, from(http://fmwww.bc.edu/repec/bocode/b)"
        }
    }
}

* outreg2 (optional – used in Module E as a fallback)
capture which outreg2
if _rc != 0 capture ssc install outreg2, replace

* ── STEP 1: Data Cleaning & Merging ─────────────────────────────────────────
display as result _newline(2) "============================================"
display as result " STEP 1: Data Cleaning"
display as result "============================================"
do "$root/1_clean_data.do"

* ── STEP 2: Nonlinear Dynamics Analysis ─────────────────────────────────────
display as result _newline(2) "============================================"
display as result " STEP 2: Nonlinear Dynamics Analysis"
display as result "============================================"
do "$root/2_nonlinear_analysis.do"

display as result _newline(2) "============================================"
display as result " ✓ PIPELINE COMPLETE"
display as result " Outputs in: $root/Cleaned Data/"
display as result "============================================"
