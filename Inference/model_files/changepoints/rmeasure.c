int alpha_IH1_int = round(alpha_IH1);
int alpha_IH2_int = round(alpha_IH2); 
int alpha_IC_int = round(alpha_IC);

// reporting covariates
const double *hosp_report = &hosp_reporting_1;
const double *icu_report = &icu_reporting_1;

// States
const double *new_deaths = &new_deaths_1_1; // latent deaths
double *ObsDeaths = &ObsDeaths_1; // observed incident deaths in each region
double *ObsHosp = &ObsHosp_1; // observed total census hospitalizations in each region
double *ObsICU = &ObsICU_1; // observed total census hospitalizations in each region
const double *IFRtrack = &IFRtrack_1_1; 
const double *IH1 = &IH1_1_1_1;
const double *IH2 = &IH2_1_1_1;
const double *IC = &IC_1_1_1;
const double *HFRtrack = &HFRtrack_1_1;
const double *DeathReportTrack = &DeathReportTrack_1_1;


for (int region=0; region<n_regions; region += 1){
    double total_hospital_reporting= hosp_report[region];
    double icu_reporting= icu_report[region];
    double death_reporting = DeathReportTrack[region];
    
    if (icu_reporting >= reporting_cap){
        icu_reporting = reporting_cap;
    }
    if (total_hospital_reporting >= reporting_cap){
        total_hospital_reporting = reporting_cap;
    }


    // Nonhospitalized deaths
    ObsDeaths[region] = rnbinom_mu(10, new_deaths[region] * death_reporting);

    // Census hospitalization
    double agg_hosp = get_sum(region, alpha_IH1_int, IH1) + get_sum(region, alpha_IH2_int, IH2);
    ObsHosp[region] = rpois(agg_hosp * total_hospital_reporting);

    // Census ICU
    double agg_icu = get_sum(region, alpha_IC_int, IC);
    ObsICU[region] = rpois(agg_icu * icu_reporting);
}