*** Moderniziation code for Beck, Levine, and Levkov (2010)
*** Author: Dui Juniper Davidsson
*** Implements: Goodman-Bacon decomposition, Callaway-Sant'Anna, Borusyak-Jaravel-Spiess, and Rambachan-Roth Honest DiD


#delimit;
clear all;
set more off;

capture log close;
log using modernization.log, replace;

cd "/Users/computerboi/Downloads/Levkovdata";

use "macro_workfile.dta", clear;



** SECTION 0: DATA PREPARATION;


*** Panel setup;
tsset statefip wrkyr;

*** Construct outcome variables (same as original Table II);
replace p10 = 1 if p10==0;
generate logistic_gini = log(gini/(1-gini));
generate log_gini      = log(gini);
generate log_theil     = log(theil);
generate log_9010      = log(p90)-log(p10);
generate log_7525      = log(p75)-log(p25);

*** Create first-treatment-year variable for modern estimators;

generate first_treat = branch_reform;

*** For CS and BJS: always-treated states must be dropped from estimation;

generate first_treat_modern = first_treat;
replace  first_treat_modern = . if first_treat < 1976;

*** Controls;
local Xs gsp_pc_growth prop_blacks prop_dropouts prop_female_headed unemploymentrate;


*** SECTION 1: Replication of Original Results (Table II, Panel A);


tabulate wrkyr, gen(wrkyr_dumm);

xtreg log_gini _intra wrkyr_dumm*, fe i(statefip) robust cluster(statefip);
estimates store twfe_loggini;
xtreg logistic_gini _intra wrkyr_dumm*, fe i(statefip) robust cluster(statefip);
estimates store twfe_logistic;

xtreg log_theil _intra wrkyr_dumm*, fe i(statefip) robust cluster(statefip);
estimates store twfe_logtheil;

xtreg log_9010 _intra wrkyr_dumm*, fe i(statefip) robust cluster(statefip);
estimates store twfe_log9010;

xtreg log_7525 _intra wrkyr_dumm*, fe i(statefip) robust cluster(statefip);
estimates store twfe_log7525;



*** SECTION 2: Goodman-Bacon Decomposition;

*** First check for balance;
xtset statefip wrkyr;

preserve;
    generate branch_reform_binned = branch_reform;
    replace branch_reform_binned = 1978 if branch_reform >= 1976 & branch_reform <= 1980;
    replace branch_reform_binned = 1983 if branch_reform >= 1981 & branch_reform <= 1985;
    replace branch_reform_binned = 1988 if branch_reform >= 1986 & branch_reform <= 1990;
    replace branch_reform_binned = 1993 if branch_reform >= 1991 & branch_reform <= 1995;
    replace branch_reform_binned = 1998 if branch_reform >= 1996 & branch_reform <= 2000;
    
    generate _intra_binned = (wrkyr >= branch_reform_binned);
    xtset statefip wrkyr;
    
*** Run bacondecomp with minimal gropt — just strip the default junk ***;
    bacondecomp log_gini _intra_binned, ddetail 
        gropt(
            legend(
                order(
                    1 "Later vs. Earlier Treated (forbidden)" 
                    2 "Earlier vs. Later Treated (clean)"
                    3 "Treated vs. Always Treated (forbidden)"
                ) 
                size(small) 
                rows(3) 
                position(6) 
                region(lcolor(white))
            )
            xtitle("Weight", size(medium))
            ytitle("2x2 DD Estimate", size(medium))
            graphregion(fcolor(white) lcolor(white) margin(medium))
            plotregion(margin(medium))
            note("")
            caption("")
            subtitle("")
        );
    
*** Stripping out the leftover "Overall DD Estimate";
    graph display, xsize(7) ysize(5);
    
    graph save bacon_loggini_full, replace;
    graph export bacon_loggini_full.png, replace width(1400) height(1000);
restore;




*** SECTION 3: Callaway-Sant'Anna (CSA) (2021);

csdid log_gini, ivar(statefip) time(wrkyr) gvar(first_treat_modern) notyet;
estimates store cs_loggini;

*** Report the overall ATT;
estat simple;

*** Report the event-study;
estat event, estore(cs_es_loggini);

*** Plot the CSA event study;
event_plot cs_es_loggini, 
    default_look
    graph_opt(
        xtitle("Years since deregulation") 
        ytitle("Effect on log Gini") 
        title("Callaway-Sant'Anna Event Study: Log Gini")
        graphregion(fcolor(white) lcolor(white))
    )
    stub_lag(Tp#) stub_lead(Tm#) together;
graph save cs_eventstudy_loggini, replace;
graph export cs_eventstudy_loggini.png, replace width(1200);


*** Overall ATT for logistic_gini;
csdid logistic_gini, ivar(statefip) time(wrkyr) gvar(first_treat_modern) notyet;
estat simple;
estat event, estore(cs_es_logistic);

*** Overall ATT for log_theil;
csdid log_theil, ivar(statefip) time(wrkyr) gvar(first_treat_modern) notyet;
estat simple;
estat event, estore(cs_es_logtheil);

*** Overall ATT for log_9010;
csdid log_9010, ivar(statefip) time(wrkyr) gvar(first_treat_modern) notyet;
estat simple;
estat event, estore(cs_es_log9010);

*** Overall ATT for log_7525;
csdid log_7525, ivar(statefip) time(wrkyr) gvar(first_treat_modern) notyet;
estat simple;
estat event, estore(cs_es_log7525);

*** CSA for every 5th percentile (modernized Figure 2);

tempfile cs_percentile_results;

generate cs_pct_b  = .;
generate cs_pct_se = .;
generate cs_pct_p  = .;
generate pctile_id = .;

local row = 0;
foreach p in 5 10 15 20 25 30 35 40 45 50 55 60 65 70 75 80 85 90 95 {;
    local row = `row' + 1;
    
*** Create log percentile variable;
    capture drop log_p`p'_temp;
    generate log_p`p'_temp = log(p`p');
    
*** Run CSA;
    quietly csdid log_p`p'_temp, ivar(statefip) time(wrkyr) gvar(first_treat_modern) notyet;
    quietly estat simple;
    
*** Store results;
    matrix b = r(b);
    matrix V = r(V);
    replace cs_pct_b  = b[1,1]          in `row';
    replace cs_pct_se = sqrt(V[1,1])    in `row';
    replace cs_pct_p  = 2*normal(-abs(b[1,1]/sqrt(V[1,1]))) in `row';
    replace pctile_id = `p'             in `row';
    
    drop log_p`p'_temp;
};

*** Plot modernized Figure 2;
preserve;
    keep if pctile_id != .;
    keep pctile_id cs_pct_b cs_pct_se cs_pct_p;
    duplicates drop;
    sort pctile_id;
    
    generate sig5 = (cs_pct_p < 0.05);
    
    twoway (bar cs_pct_b pctile_id if sig5==1, sort fcolor(navy) lcolor(navy) barwidth(3))
           (bar cs_pct_b pctile_id if sig5==0, sort fcolor(navy) lcolor(navy) barwidth(3) fintensity(30)),
           ytitle("Effect on log income (CSA ATT)") ytitle(, size(small))
           ylabel(, labsize(small) angle(horizontal) nogrid)
           xtitle("Percentile of income distribution") xtitle(, size(small))
           xlabel(5(5)95, labsize(small))
           legend(order(1 "Significant at 5%" 2 "Not significant") size(small))
           title("Callaway-Sant'Anna: Effect by Percentile (Modernized Figure 2)")
           graphregion(fcolor(white) lcolor(white));
    graph save cs_figure2, replace;
    graph export cs_figure2.png, replace width(1200);
restore;

drop cs_pct_b cs_pct_se cs_pct_p pctile_id;

*** SECTION 4: Borusyak-Jaravel-Spiess (BJS) (2024);

preserve;
    drop if first_treat < 1976;
    
    *** Overall ATTs from BJS;
    did_imputation log_gini statefip wrkyr first_treat, cluster(statefip) autosample;
    estimates store bjs_att_loggini;
    did_imputation logistic_gini statefip wrkyr first_treat, cluster(statefip) autosample;
    estimates store bjs_att_logistic;
    did_imputation log_theil statefip wrkyr first_treat, cluster(statefip) autosample;
    estimates store bjs_att_logtheil;
    did_imputation log_9010 statefip wrkyr first_treat, cluster(statefip) autosample;
    estimates store bjs_att_log9010;
    did_imputation log_7525 statefip wrkyr first_treat, cluster(statefip) autosample;
    estimates store bjs_att_log7525;
    
** Event studies (keep existing code);
    did_imputation log_gini statefip wrkyr first_treat, 
        horizons(0/15) pretrends(10) cluster(statefip) autosample;
    estimates store bjs_loggini;
    
restore;

preserve;
    drop if first_treat < 1976;
    
    *** Event study from BJS;
    did_imputation log_gini statefip wrkyr first_treat, 
        horizons(0/15) pretrends(10) cluster(statefip) autosample;
    estimates store bjs_loggini;
    
    event_plot bjs_loggini,
        default_look
        graph_opt(
            xtitle("Years since deregulation")
            ytitle("Effect on log Gini")
            title("BJS Imputation Event Study: Log Gini")
            graphregion(fcolor(white) lcolor(white))
        )
        stub_lag(tau#) stub_lead(pre#) together;
    graph save bjs_eventstudy_loggini, replace;
    graph export bjs_eventstudy_loggini.png, replace width(1200);
    
    did_imputation logistic_gini statefip wrkyr first_treat, 
        horizons(0/15) pretrends(10) cluster(statefip) autosample;
    estimates store bjs_logistic;
    
    did_imputation log_theil statefip wrkyr first_treat, 
        horizons(0/15) pretrends(10) cluster(statefip) autosample;
    estimates store bjs_logtheil;
    
    did_imputation log_9010 statefip wrkyr first_treat, 
        horizons(0/15) pretrends(10) cluster(statefip) autosample;
    estimates store bjs_log9010;
    
    did_imputation log_7525 statefip wrkyr first_treat, 
        horizons(0/15) pretrends(10) cluster(statefip) autosample;
    estimates store bjs_log7525;
    
restore;

*** SECTION 5: Combined Event Study Plot (TWFE vs CSA vs BJS);

*** Re-run original TWFE event study (from Figure3.do);
generate _tintra = wrkyr - branch_reform;
replace _tintra = -10 if _tintra < -10;
replace _tintra = 15  if _tintra > 15;

forvalues j = 10(-1)1 {;
    generate dm`j' = (_tintra == -`j');
};
forvalues j = 1/15 {;
    generate dp`j' = (_tintra == `j');
};

xtreg log_gini dm10-dm1 dp1-dp15 wrkyr_dumm*, fe i(statefip) robust cluster(statefip);
estimates store twfe_es_loggini;

*** Plot all three estimates ***;
event_plot cs_es_loggini,
    default_look
    graph_opt(
        xtitle("Years since deregulation")
        ytitle("Effect on log Gini")
        title("Callaway-Sant'Anna")
        graphregion(fcolor(white) lcolor(white))
        yline(0, lcolor(gs10) lpattern(dash))
    )
    stub_lag(Tp#) stub_lead(Tm#) together;
graph save panel_cs, replace;

event_plot bjs_loggini,
    default_look
    graph_opt(
        xtitle("Years since deregulation")
        ytitle("Effect on log Gini")
        title("BJS Imputation")
        graphregion(fcolor(white) lcolor(white))
        yline(0, lcolor(gs10) lpattern(dash))
    )
    stub_lag(tau#) stub_lead(pre#) together;
graph save panel_bjs, replace;

event_plot twfe_es_loggini,
    default_look
    graph_opt(
        xtitle("Years since deregulation")
        ytitle("Effect on log Gini")
        title("TWFE (Original)")
        graphregion(fcolor(white) lcolor(white))
        yline(0, lcolor(gs10) lpattern(dash))
    )
    stub_lag(dp#) stub_lead(dm#) together;
graph save panel_twfe, replace;

graph combine panel_twfe.gph panel_cs.gph panel_bjs.gph,
    rows(1) cols(3)
    title("Event Study Comparison: TWFE vs CS vs BJS")
    graphregion(fcolor(white) lcolor(white))
    xsize(12) ysize(4);
graph save combined_panels, replace;
graph export combined_panels.png, replace width(1800);

drop dm1 dm2 dm3 dm4 dm5 dm6 dm7 dm8 dm9 dm10 dp1 dp2 dp3 dp4 dp5 dp6 dp7 dp8 dp9 dp10 dp11 dp12 dp13 dp14 dp15 _tintra;


*** SECTION 6: RAMBACHAN AND ROTH (2023) HONEST DiD;

*** Re-run CS and store the event study;
csdid log_gini, ivar(statefip) time(wrkyr) gvar(first_treat_modern) notyet;
estat event, estore(cs_for_honest);

*** Run Honest DiD with relative magnitudes approach;
*** M = 0: linear extrapolation of pre-trend;
*** M = 1: post-treatment violations no bigger than pre-treatment;
*** M = 2: post-treatment violations up to twice as big;

*** Compute honest CI for the average post-treatment effect;
honestdid, pre(1/10) post(1/15) mvec(0 0.5 1 1.5 2)
    delta(rm) alpha(0.05);

*** SECTION 7: TWFE Robustness — Drop Always-Treated;

*** Re-run Table II dropping the 12 always-treated states;
*** This checks how much always-treated states influence the TWFE estimate;

preserve;
    drop if first_treat < 1976;
    
    tabulate wrkyr, gen(yr_dumm_sub);
    
    xtreg log_gini _intra yr_dumm_sub*, fe i(statefip) robust cluster(statefip);
    estimates store twfe_loggini_drop;
    
    xtreg logistic_gini _intra yr_dumm_sub*, fe i(statefip) robust cluster(statefip);
    estimates store twfe_logistic_drop;
    
    xtreg log_theil _intra yr_dumm_sub*, fe i(statefip) robust cluster(statefip);
    estimates store twfe_logtheil_drop;
    
    xtreg log_9010 _intra yr_dumm_sub*, fe i(statefip) robust cluster(statefip);
    estimates store twfe_log9010_drop;
    
    xtreg log_7525 _intra yr_dumm_sub*, fe i(statefip) robust cluster(statefip);
    estimates store twfe_log7525_drop;
    
    estout twfe_logistic_drop twfe_loggini_drop twfe_logtheil_drop twfe_log9010_drop twfe_log7525_drop,
    keep(_intra)
    cells(b(star fmt(3)) se(par)) stats(r2 N, fmt(2 0))
    starlevel(* 0.10 ** 0.05 *** 0.01)
    title("Table II Robustness: TWFE Dropping Always-Treated States");
    
restore;


*** SECTION 8: CSA Sensitivity to Control Group Composition;


*** Baseline CSA — full sample, all not-yet-treated as controls;
csdid log_gini, ivar(statefip) time(wrkyr) gvar(first_treat_modern) notyet;
estat simple;
estimates store cs_baseline;

*** CSA restricted to pre-1994 window;
tempfile fulldata_restricted;
save `fulldata_restricted', replace;

*** Create a restricted treatment variable: cohorts treated after 1993;
generate first_treat_restricted = first_treat_modern;
replace first_treat_restricted = . if first_treat_modern > 1993;

*** Drop observations after 1993;
drop if wrkyr > 1993;

*** Run CSA on the restricted panel;
csdid log_gini, ivar(statefip) time(wrkyr) gvar(first_treat_restricted) notyet;
estat simple;
estimates store cs_pre1994;

*** Restore full data;
use `fulldata_restricted', clear;

*** CSA restricted to pre-1991 window;
tempfile fulldata_restricted2;
save `fulldata_restricted2', replace;

generate first_treat_restricted2 = first_treat_modern;
replace first_treat_restricted2 = . if first_treat_modern > 1990;

drop if wrkyr > 1990;

csdid log_gini, ivar(statefip) time(wrkyr) gvar(first_treat_restricted2) notyet;
estat simple;
estimates store cs_pre1991;


log close;
