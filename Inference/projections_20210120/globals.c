//Function to make a line
double get_par_value( double p0, double pf, double t0, double tf, double t_now ) {
   double slope = (pf - p0) / (tf - t0);
   double value;
   if ( t_now >= tf){
       value = pf;
   } else if (t_now <= t0){
       value = p0;
   } else{
       value = p0 + slope * (t_now - t0);
   }
  return value;
}
// Function untransform a logit parameter
double unlogit( double logit_value, double pmax, double pmin) {
  return (logit_value * (pmax - pmin)) + pmin;
}

// Function to get the sum of a series of subdivided compartments
int get_sum(int reg, int alpha_int, const double* state){
  double agg = 0;
    int Start = reg * alpha_int;
    for (int k=0; k<alpha_int; k++){
        int offset = Start + k;
        agg += state[offset];
    }
    return(agg);
}

int get_changepoint_index(int reg, const double* changepoints_){
  int start = 0;
  for (int r=1; r<=reg; r++){
      int z = round(changepoints_[r-1]);
      start = round(start + z);
  }
  return start;
}

double calc_time_varying_param( int coffset, 
	int n_timepoints, 
	double t_now, 
	const double* changepoint_vals, 
	const double* vals, 
	const double vmax, 
	const double vmin){
  double final_val;
  // if we're at a point before the first timepoint, just assign to first ICU value
  if (t_now <= changepoint_vals[coffset]){
    final_val = unlogit(vals[coffset], vmax, vmin);
  } else if(t_now >= changepoint_vals[coffset + n_timepoints - 1]){
    final_val = unlogit(vals[coffset + n_timepoints - 1], vmax, vmin);
  } else{
    for (int i=coffset + 1; i <= coffset + n_timepoints - 1; i++){
      if(t_now < changepoint_vals[i]){
        double b1 = unlogit(vals[i-1], vmax, vmin);
        double b2 = unlogit(vals[i], vmax, vmin);
        double t1 = changepoint_vals[i-1];
        double t2 = changepoint_vals[i];
        final_val = get_par_value(b1, b2, t1, t2, t_now);
        break;
      }
    }
  }
  return final_val;
}

double dnbinom_states(const double* ObsState, double agg_state, int reg, double reporting_prob, double dispersion){
	double local_lik = 0;
	if (ISNA(ObsState[reg])) {
      local_lik = 0;
    }
    else{   
      local_lik = dnbinom_mu(ObsState[reg], dispersion, agg_state * reporting_prob, 1);
     }
     return(local_lik);
}

double dpois_states(const double* ObsState, double agg_state, int reg, double reporting_prob){
	double local_lik = 0;
	if (ISNA(ObsState[reg])) {
      local_lik = 0;
    }
    else{   
      local_lik = dpois(ObsState[reg], agg_state * reporting_prob, 1);
    }
    return(local_lik);
}

double dbinom_states(const double* ObsState, double agg_state, int reg, double reporting_prob){
	double local_lik = 0;
	if (ISNA(ObsState[reg])) {
      local_lik = 0;
    }
    else{
		if(agg_state < ObsState[reg]){
	        local_lik = -1e3;
	    } else{
	        local_lik = dbinom(ObsState[reg], agg_state, reporting_prob, 1);
	    }
     }
     return(local_lik);
}

void update_states(int reg, int alpha_int, double* dState, double* State, double DT, double rate, double inflow){
	int Start = reg * alpha_int;
	double Pr = 1-exp(-rate * alpha_int * DT);
	dState[0] = rbinom(State[Start], Pr);
	State[Start] += inflow - dState[0];
	for (int i=1; i<alpha_int; i++){
		int offset = Start + i;
		dState[i] = rbinom(State[offset], Pr);
		State[offset] += dState[i - 1] - dState[i];
	}
}