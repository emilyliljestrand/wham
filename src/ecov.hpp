#include <limits>

// Check dimensions before any latent-covariate or process-parameter access.
template<class Type>
void validate_Ecov_re(vector<int> Ecov_model, matrix<Type> Ecov_process_pars,
  matrix<Type> Ecov_re, vector<int> Ecov_use_re){
  int n_Ecov = Ecov_model.size();
  if(Ecov_use_re.size() != n_Ecov)
    Rf_error("Ecov_use_re must have one entry per Ecov");
  for(int i = 0; i < n_Ecov; i++) if(Ecov_use_re(i) == 1 && Ecov_model(i) > 0){
    int n_pars = Ecov_model(i) == 2 ? 3 : 2;
    if(Ecov_re.cols() <= i || Ecov_process_pars.cols() <= i || Ecov_process_pars.rows() < n_pars)
      Rf_error("Ecov random effects or process parameters have incompatible dimensions");
  }
}

template<class Type>
void validate_Ecov_years(vector<int> Ecov_model, matrix<Type> Ecov_re,
  vector<int> Ecov_use_re, vector<int> years_use){
  for(int y = 0; y < years_use.size(); y++){
    if(years_use(y) < 0 || years_use(y) >= Ecov_re.rows())
      Rf_error("Ecov years_use contains an out-of-bounds year");
    if(y >= 2 && years_use(y) == 0){
      for(int i = 0; i < Ecov_model.size(); i++)
        if(Ecov_model(i) == 1 && Ecov_use_re(i) == 1)
          Rf_error("Ecov random walk requires a preceding year");
    }
  }
}

template<class Type>
matrix<Type> get_nll_Ecov(vector<int> Ecov_model, matrix<Type> Ecov_process_pars, matrix<Type> Ecov_re, vector<int> Ecov_use_re, 
  vector<int> years_use){
  /* 
     get nll contribtions for any Ecovs
       Ecov_model: which time series model to use (RW, AR1)
        Ecov_process_pars: parameters for time series models for each Ecov (columns)
          Ecov_re: random effects representing latent covariates
      Ecov_use_re: whether to include random effects in the likelihood
        years_use: possibly a subset of years to use for evaluating likelihood (and simulating values). normally = 0,....,n_years_model-1
  */
  using namespace density; // necessary to use AR1, SCALE, SEPARABLE
  int n_Ecov = Ecov_model.size();
  int n_y = years_use.size();
  validate_Ecov_re(Ecov_model, Ecov_process_pars, Ecov_re, Ecov_use_re);
  validate_Ecov_years(Ecov_model, Ecov_re, Ecov_use_re, years_use);
  //int n_y = Ecov_re.rows();
  matrix<Type> nll_Ecov(n_y, n_Ecov); // nll contribution each Ecov_re
  nll_Ecov.setZero();

  for(int i = 0; i < n_Ecov; i++){ // loop over Ecovs
    if(n_y == 0 || (Ecov_model(i) == 1 && n_y < 2)) continue;
    // Ecov model option 1: RW
    if((Ecov_model(i) == 1) & (Ecov_use_re(i)==1)){
      Type Ecov_sig = exp(Ecov_process_pars(1,i)); // sd (sig_x in Eq1, pg 1262, Miller et al. 2016)
      Type Ecov1 = Ecov_process_pars(0,i);; // Ecov_x in year 1 (fixed effect)
      nll_Ecov(1,i) -= dnorm(Ecov_re(years_use(1),i), Ecov1, Ecov_sig, 1); // Ecov_re(0,i) set to NA
      for(int y = 2; y < n_y; y++){
        nll_Ecov(y,i) -= dnorm(Ecov_re(years_use(y),i), Ecov_re(years_use(y)-1,i), Ecov_sig, 1);
      }
    }

    // Ecov model option 2: AR1
    if((Ecov_model(i) == 2) & (Ecov_use_re(i)==1)){
      //Type Ecov_mu = Ecov_process_pars(0,i);  // mean
      Type Ecov_phi = geninvlogit(Ecov_process_pars(2,i), Type(-1), Type(1), Type(1)); // autocorrelation
      //Type Ecov_phi = -1 + 2/(1 + exp(-Ecov_process_pars(2,i))); // autocorrelation
      Type Ecov_sig = exp(Ecov_process_pars(1,i)) *pow(1 - pow(Ecov_phi,2),-0.5); // marginal sd
      vector<Type> re_i(n_y);
      for(int y = 0; y < n_y; y++) re_i(y) = Ecov_re(years_use(y),i);
      //vector<Type> re_i = Ecov_re.col(i);
      nll_Ecov(0,i) += SCALE(AR1(Ecov_phi), Ecov_sig)(re_i);
    }
  }
  return(nll_Ecov);
}

template<class Type>
matrix<Type> get_Ecov(vector<int> Ecov_model, matrix<Type> Ecov_process_pars, matrix<Type> Ecov_re, vector<int> Ecov_use_re){
  /* 
     get Ecov time series spanning years of environmental data and population model
       Ecov_model: which time series model to use (RW, AR1)
        Ecov_process_pars: parameters for time series models for each Ecov (columns)
          Ecov_re: random effects representing latent covariates
      Ecov_use_re: whether to include random effects in the likelihood
  */

  int n_Ecov = Ecov_model.size();
  int n_y = Ecov_re.rows();
  validate_Ecov_re(Ecov_model, Ecov_process_pars, Ecov_re, Ecov_use_re);
  if(n_y == 0) for(int i = 0; i < n_Ecov; i++)
    if(Ecov_model(i) == 1 && Ecov_use_re(i) == 1)
      Rf_error("Ecov random walk requires at least one process year");
  matrix<Type> Ecov_x(n_y, n_Ecov); // 'true' estimated Ecov (x_t in Miller et al. 2016 CJFAS)
  Ecov_x.setZero();

  // Ecov_model == 0) no Ecov
  for(int i = 0; i < n_Ecov; i++){ // loop over Ecovs
    // Ecov model option 1: RW
    if((Ecov_model(i) == 1) & (Ecov_use_re(i)==1)){
      Type Ecov1 = Ecov_process_pars(0,i); // Ecov_x in year 1 (fixed effect)
      Ecov_x(0,i) = Ecov1;
      for(int y = 1; y < n_y; y++) Ecov_x(y,i) = Ecov_re(y,i);
    }

    // Ecov model option 2: AR1
    if((Ecov_model(i) == 2) & (Ecov_use_re(i)==1)){
      Type Ecov_mu = Ecov_process_pars(0,i); // mean
      for(int y = 0; y < n_y; y++) Ecov_x(y,i) = Ecov_mu + Ecov_re(y,i);
    }
  } // end loop over Ecovs
  return(Ecov_x);
}

template<class Type>
matrix<Type> simulate_Ecov_re(vector<int> Ecov_model, matrix<Type> Ecov_process_pars, matrix<Type> Ecov_re, vector<int> Ecov_use_re, 
  vector<int> years_use){
  /* 
     simulate random effects for Ecov time series
       Ecov_model: which time series model to use (RW, AR1)
        Ecov_process_pars: parameters for time series models for each Ecov (columns)
          Ecov_re: random effects representing latent covariates
      Ecov_use_re: whether to simulate random effects (and whether to include likelihood contributions)
        years_use: possibly a subset of years to use for evaluating likelihood (and simulating values). normally = 0,....,n_years_model-1
  */
  using namespace density; // necessary to use AR1, SCALE, SEPARABLE
  int n_Ecov = Ecov_model.size();
  int n_y = years_use.size();
  validate_Ecov_re(Ecov_model, Ecov_process_pars, Ecov_re, Ecov_use_re);
  validate_Ecov_years(Ecov_model, Ecov_re, Ecov_use_re, years_use);
  //int n_y = Ecov_re.rows();
  matrix<Type> re_sim = Ecov_re;

  for(int i = 0; i < n_Ecov; i++){ // loop over Ecovs
    if(n_y == 0 || (Ecov_model(i) == 1 && n_y < 2)) continue;
    if(Ecov_use_re(i)==1){
      // Ecov model option 1: RW
      if((Ecov_model(i) == 1) & (Ecov_use_re(i)==1)){
        Type Ecov1 = Ecov_process_pars(0,i);; // Ecov_x in year 1 (fixed effect)
        Type Ecov_sig = exp(Ecov_process_pars(1,i)); // sd (sig_x in Eq1, pg 1262, Miller et al. 2016)
        re_sim(years_use(1),i) = rnorm(Ecov1, Ecov_sig);
        for(int y = 2; y < n_y; y++){
          re_sim(years_use(y),i) = rnorm(re_sim(years_use(y)-1,i), Ecov_sig);
        }
      }
      // Ecov model option 2: AR1
      if((Ecov_model(i) == 2) & (Ecov_use_re(i)==1)){
        Type Ecov_phi = geninvlogit(Ecov_process_pars(2,i), Type(-1), Type(1), Type(1)); // autocorrelation
        //Type Ecov_phi = -1 + 2/(1 + exp(-Ecov_process_pars(2,i))); // autocorrelation
        Type Ecov_sig = exp(Ecov_process_pars(1,i)) *pow(1 - pow(Ecov_phi,2),-0.5); // marginal sd
        //vector<Type> re_i = re_sim.col(i);
        vector<Type> re_i(n_y);
        AR1(Ecov_phi).simulate(re_i);
        re_i *= Ecov_sig;
        for(int y = 0; y < n_y; y++) re_sim(years_use(y),i) = re_i(y);
      }
    }
  }
  return(re_sim);
}

template<class Type>
matrix<Type> get_Ecov_out(matrix<Type> Ecov_x, int n_years_model, int n_years_proj, 
  vector<int> ind_Ecov_out_start, vector<int> ind_Ecov_out_end, vector<int> proj_Ecov_opt, vector<int> avg_years, matrix<Type> Ecov_use_proj){
  /* 
     get Ecov time series for years pertinent to specific effects on population
                    Ecov_x: Ecov time series spanning years of environmental data and population model 
             n_years_model: number of years for population model
              n_years_proj: number of years population model is projected
        ind_Ecov_out_start: which Ecov_x year is first year needed for effect on population
          ind_Ecov_out_end: which Ecov_x years is last year needed for effect on population
             proj_Ecov_opt: how to project each Ecov
                 avg_years: if proj_Ecov_opt = 2, which Ecov years to average over
             Ecov_use_proj: if proj_Ecov_opt = 4, (n_years_proj x n_Ecov) user-supplied fixed values to project Ecov at.
  */
  //int n_effects = Ecov_beta.dim(0); // 2 + n_indices (recruitment, mortality and any catchabilities)
  int n_Ecov = Ecov_x.cols();
  if(n_years_model <= 0 || n_years_proj < 0 ||
     n_years_proj > std::numeric_limits<int>::max() - n_years_model)
    Rf_error("Ecov output requires positive model years and nonnegative projection years");
  if(ind_Ecov_out_start.size() != n_Ecov || ind_Ecov_out_end.size() != n_Ecov ||
     (n_years_proj > 0 && proj_Ecov_opt.size() != n_Ecov))
    Rf_error("Ecov output index and projection vectors must have one entry per Ecov");
  for(int i = 0; i < n_Ecov; i++){
    int start = ind_Ecov_out_start(i), end = ind_Ecov_out_end(i);
    if(start < 0 || end < start || end >= Ecov_x.rows() ||
       end - start != n_years_model - 1)
      Rf_error("Ecov output endpoints must be in bounds and span n_years_model years");
    if(n_years_proj > 0){
      int opt = proj_Ecov_opt(i);
      if(opt < 1 || opt > 4) Rf_error("Invalid Ecov projection option");
      if(opt == 1 && n_years_proj > Ecov_x.rows() - 1 - end)
        Rf_error("Ecov process does not cover projection years");
      if(opt == 2){
        if(avg_years.size() == 0) Rf_error("Ecov projection averaging requires years");
        for(int y = 0; y < avg_years.size(); y++)
          if(avg_years(y) < 0 || avg_years(y) >= Ecov_x.rows())
            Rf_error("Ecov projection averaging year is out of bounds");
      }
      if(opt == 4 && (Ecov_use_proj.rows() < n_years_proj || Ecov_use_proj.cols() <= i))
        Rf_error("Ecov_use_proj does not cover requested projection years and covariates");
    }
  }
  matrix<Type> Ecov_out(n_years_model + n_years_proj, n_Ecov); // Pop model uses Ecov_out(t) for processes in year t (Ecov_x shifted by lag and padded)
  Ecov_out.setZero(); // set Ecov_out = 0
  for(int i = 0; i < n_Ecov; i++){
      int ct = 0;
      for(ct = 0; ct < n_years_model + n_years_proj; ct++){
        if(ct < n_years_model) {
          Ecov_out(ct,i) = Ecov_x(ind_Ecov_out_start(i) + ct,i);
        } else { //projection year
          if(proj_Ecov_opt(i) == 1) Ecov_out(ct,i) = Ecov_x(ind_Ecov_out_start(i) + ct,i); //continue Ecov process
          if(proj_Ecov_opt(i) == 2) { //average over avg_years
            // for(int ind = 0; ind < avg_years.size(); ind++) Ecov_out(ct,i) += Ecov_out(avg_years(ind),i)/Type(avg_years.size());
            for(int ind = 0; ind < avg_years.size(); ind++) Ecov_out(ct,i) += Ecov_x(avg_years(ind),i)/Type(avg_years.size());
          }
          if(proj_Ecov_opt(i) == 3) Ecov_out(ct,i) = Ecov_out(n_years_model-1,i); //use terminal year
          if(proj_Ecov_opt(i) == 4) Ecov_out(ct,i) = Ecov_use_proj(ct-n_years_model,i); //user-defined ecov
        }
      }
  }
  return Ecov_out;
}

template<class Type>
matrix<Type> get_Ecov_lm(matrix<Type>Ecov_beta, matrix<Type>Ecov_out, int n_years_model, int n_years_proj){
  /*
    Calculate ecov link model (possibly othogonal polynomial: b1*ecov + b2*ecov^2 + ...) --------------------
    to be added to linear model for particular parameter effects (R, q, M, mu, etc.).
    Ecov_beta: matrix (a sub-component of ecov_beta_*): n_ecov x n_poly
    Ecov_out: matrix returned by get_Ecov_out with Ecov time series spanning pertinent years
  */
  int n_Ecov = Ecov_out.cols();
  int n_poly = Ecov_beta.cols();
  // Ecov_lm stores the linear models for each Ecov each year.
  int ny = Ecov_out.rows(); //n_years_model + n_years_proj 
  matrix<Type> Ecov_lm(ny,n_Ecov); 
  Ecov_lm.setZero();
  for(int i = 0; i < n_Ecov; i++){
    vector<Type> thecol = Ecov_out.col(i);
    matrix<Type> X_poly(ny, n_poly);
    X_poly.setZero();
    if(n_poly == 1){ // n_poly = 1 if ecov effect is none or linear
      X_poly = thecol.matrix();
    } else { // n_poly > 1, get poly transformation for ith ecov
      X_poly = poly_trans(thecol, n_poly, n_years_model, n_years_proj); //need to treat X in projections differently
    }
    for(int y = 0; y < ny; y++){
      for(int j = 0; j < n_poly; j++){
        Ecov_lm(y,i) += Ecov_beta(i,j) * X_poly(y,j); // poly transformation returns design matrix, don't need to take powers
      }
    }
  }
  return Ecov_lm;
}

template<class Type>
matrix<Type> get_Ecov_lm(matrix<Type>Ecov_beta, matrix<Type>Ecov_out, int n_years_model, int n_years_proj, vector<int> n_poly_Ecov){
  /*
    Calculate ecov link model (possibly othogonal polynomial: b1*ecov + b2*ecov^2 + ...) --------------------
    to be added to linear model for particular parameter effects (R, q, M, mu, etc.).
    Ecov_beta: matrix (a sub-component of ecov_beta_*): n_ecov x n_poly
    Ecov_out: matrix returned by get_Ecov_out with Ecov time series spanning pertinent years
    n_poly_Ecov: vector (n_Ecov) for the polynomial order for each Ecov
  */
  int n_Ecov = Ecov_out.cols();
  //int n_poly = Ecov_beta.cols();
  // Ecov_lm stores the linear models for each Ecov each year.
  int ny = Ecov_out.rows(); //n_years_model + n_years_proj 
  matrix<Type> Ecov_lm(ny,n_Ecov); 
  Ecov_lm.setZero();
  for(int i = 0; i < n_Ecov; i++){
    vector<Type> thecol = Ecov_out.col(i);
    matrix<Type> X_poly(ny, n_poly_Ecov(i));
    //see(n_poly_Ecov(i));
    //see(ny);
    //see(n_years_model);
    //see(n_years_proj);
    X_poly.setZero();
    if(n_poly_Ecov(i) == 1){ // n_poly = 1 if ecov effect is none or linear
      X_poly = thecol.matrix();
    } else { // n_poly > 1, get poly transformation for ith ecov
      X_poly = poly_trans(thecol, n_poly_Ecov(i), n_years_model, n_years_proj); //need to treat X in projections differently
    }
    //see("Xpoly_done");
    //see(Ecov_beta.rows());
    //see(Ecov_beta.cols());
    //see(X_poly);
    for(int y = 0; y < ny; y++){
      for(int j = 0; j < n_poly_Ecov(i); j++){
        Ecov_lm(y,i) += Ecov_beta(i,j) * X_poly(y,j); // poly transformation returns design matrix, don't need to take powers
      }
    }
  }
  return Ecov_lm;
}

