//=========================================================================
// Dynare replication code for den Haan, Freund and Rendahl (2021)
// Extension to risk-aversion
// Last updated: September 2020
// Please email lukas.beat.freund@gmail.com with any questions
//=========================================================================

//=========================================================================
//User Options
//=========================================================================

% Pre-processor variable to enable computation of ergodic moments (slows down algorithm)
@# define OptionMoments = 0

//--------------------------------------------------------------------------

% Pre-processor variable to switch b/w Nash (1) and linear (0)
@# define OptionWNash = 1

% Pre-processor variable to switch to remove SDF from Nash wage (first term)
@# define OptionNashSDFConstant = 0

% Pre-processor variable to switch to hold MU at SS in Nash wage (second term)
@# define OptionNashMUConstant = 0

//--------------------------------------------------------------------------

% Pre-processor variable to switch b/w baseline timing (1) and LL alternative (0) 
@# define OptionTiming = 1

//=========================================================================
//DECLARATION OF ENDOGENOUS VARIABLES
//=========================================================================
var
    c       ${c}$               (long_name='Consumption')
    y       ${y}$               (long_name='Aggregate output')
    m       ${m}$               (long_name='Number of matches')                
    u       ${u}$               (long_name='Number of unemployed')
    us      ${u^s}$             (long_name='Number of searching workers')
    v       ${v}$               (long_name='Number of vacancies')
    f       ${f}$               (long_name='Job finding rate')
    h       ${h}$               (long_name='Vacancy filling rate')
    n       ${n}$               (long_name='Employment')
    J       ${J}$               (long_name='Value of a firm (with match)')
    w       ${w}$               (long_name='Actual wage')
    wN      ${w^N}$             (long_name='Nash wage')
    wL      ${w^{L}}$           (long_name='Linear wage')
    theta   ${\theta}$          (long_name='Labor market tightness')   
    z       ${z}$               (long_name='Productivity')
    sigma_z ${\sigma_z}$        (long_name='Technology shock volatility')

// The following are added for the extension to risk aversion...
    lambda  ${\lambda}$         (long_name='Marginal utility')
    RReal   ${R^r}$             (long_name='Real interest rate (gross)')

                   
// The following are auxililiary variables
    n_Obs 
    z_Obs 
    sigma_z_Obs
;
 
//=========================================================================
//DECLARATION OF EXOGENOUS VARIABLES
//=========================================================================
varexo
    eps_z   ${\epsilon_z}$          (long_name='Technology innovation')
    eps_s   ${\epsilon_{\sigma_z}}$ (long_name='Technology uncertainty innovation')
;
 
//=========================================================================
//DECLARATION OF PARAMETERS       
//=========================================================================
parameters
    betta        ${\beta}$               (long_name='Subjective discount factor')
    eta         ${\eta}$                (long_name='Elasticity of substitution b/w differentiated products')
    alpha       ${\alpha}$              (long_name='Elasticity of matching w.r.t. u')
    delta       ${\rho}$                (long_name='Exogenous job destruction rate')      
    kappa       ${\kappa}$              (long_name='Vacancy posting fixed cost')
    omegaL      ${\omega^L}$            (long_name='Bargaining weight for workers, linear wage')
    omegaN      ${\omega^N}$            (long_name='Bargaining weight for workers, Nash')
    rho_z       ${\rho_z}$              (long_name='Technology shock, persistence')
    rho_sigma_z ${\rho_{\sigma_z}}$     (long_name='Technoloy uncertainty shock, persistence')
    sigma_zbar  ${\bar{\sigma}_z}$      (long_name='Technology uncertainty shock, mean value')
    sigma_sigma_z ${\sigma_{\sigma_z}}$ (long_name='Technology uncertainty shock, sd of innovation')
    chiN         ${\bar{\chiN}}$        (long_name='Disutility of working')
    chiL         ${\chi}$               (long_name='Strike value')
    psi          ${\mu}$                (long_name='Matching efficiency')
    xss          ${\bar{x}}$             (long_name='Relative price of intermediates')

// The following are added for the extension to risk aversion
    xi          ${\xi}$                 (long_name='Inverse elasticity of intertemporal substitution')
    phi         ${\phi}$                (long_name='Flow unemployment benefits')
    lambdass

// The following are calibrated steady-state values of variables
    IndN        ${\mathbb{I}^N}$        (long_name='Indicator variable for wage process')
    zBar        ${\bar{Z}}$             (long_name='Technology shock, mean value')
;

//=========================================================================
//PARAMETER VALUES                
//=========================================================================


load inputParameters

set_param_value('xi',sPar.xi);
set_param_value('phi',sPar.phi);
set_param_value('zBar',sPar.zBar);
set_param_value('xss',sPar.xss);
set_param_value('betta',sPar.beta);
set_param_value('eta',sPar.eta);
set_param_value('alpha',sPar.alpha);
set_param_value('delta',sPar.delta);
set_param_value('kappa',sPar.kappa);
set_param_value('psi',sPar.psi);
set_param_value('lambdass',sPar.ss.lambda);

% Productivity shocks
set_param_value('rho_z',sPar.rho_z);
set_param_value('sigma_zbar',sPar.sigma_zbar);
set_param_value('sigma_sigma_z',sPar.sigma_sigma_z);
set_param_value('rho_sigma_z',sPar.rho_sigma_z);

@#if OptionWNash > 0  
IndN = 1;
@#else
IndN = 0;
@#endif 

//--------------------------------------------------------------------------
//    SS Relationships
//--------------------------------------------------------------------------

% Wage parameters are set in main file
set_param_value('omegaN',sPar.omegaN);
set_param_value('omegaL',sPar.omegaL);
set_param_value('chiN',sPar.chiN);
set_param_value('chiL',sPar.chiL);

//=========================================================================
//MODEL EQUATIONS                 
//=========================================================================
model;

//--------------------------------------------------------------------------
//   Representative household
//--------------------------------------------------------------------------
[name = 'Marginal utility']
lambda = (c)^(-xi);

[name = 'Bond Euler']
1 = betta*RReal*(lambda(+1)/lambda);

//--------------------------------------------------------------------------
//   Labor market
//--------------------------------------------------------------------------

[name = 'Searching workers']
us = 1-(1-delta)*n(-1);

[name = 'Matching function']
m = psi*us^alpha*v^(1-alpha);

[name = 'Job finding rate ']
f= m/us;

[name = 'Vacancy filling rate']
h = m/v;

[name = 'Labor market tightness']
theta = v/us;

[name = 'Employment law of motion']
n = (1-delta)*n(-1)+f*us;   

[name = 'Unemployment rate']
u = 1-n;
 
[name = 'Wage, linear']
wL = omegaL*(xss*z)+(1-omegaL)*chiL; 

[name = 'Wage, Nash']
@#if OptionNashSDFConstant >0
@#if OptionNashMUConstant >0
wN = omegaN*(xss*z+(1-delta)*betta*(kappa*theta(+1)))+(1-omegaN)*(phi + chiN/lambdass);
@#else
wN = omegaN*(xss*z+(1-delta)*betta*(kappa*theta(+1)))+(1-omegaN)*(phi + chiN/lambda);
@#endif
@#else
@#if OptionNashMUConstant >0
wN = omegaN*(xss*z+(1-delta)*betta*(lambda(+1)/lambda)*(kappa*theta(+1)))+(1-omegaN)*(phi + chiN/lambdass);
@#else
wN = omegaN*(xss*z+(1-delta)*betta*(lambda(+1)/lambda)*(kappa*theta(+1)))+(1-omegaN)*(phi + chiN/lambda); % default with both effects
@#endif
@#endif

[name = 'Actual wage']
 w = IndN*wN+(1-IndN)*wL;
//--------------------------------------------------------------------------
//  Intermediate goods sector
//--------------------------------------------------------------------------

[name = 'Production Function']
y = z*n;

[name = 'Match value']
kappa = h*J;

[name = 'Equity Euler']
J = xss*z-w+(1-delta)*betta*(lambda(+1)/lambda)*J(+1); % default with both effects

//--------------------------------------------------------------------------
//  Aggregation
//--------------------------------------------------------------------------

[name = 'Resource constraint']
y = c;  % To keep things reasonably transparent, I immediately imposed that vancancy 
% posting resource costs are rebated to the household


//--------------------------------------------------------------------------
//   Shock Processes
//--------------------------------------------------------------------------
 
[name = 'Productivity shock']
@#if OptionTiming > 0
z = (1-rho_z)* zBar+rho_z*z(-1)+sigma_z(-1)*eps_z;
@#else 
z = (1-rho_z)* zBar+rho_z*z(-1)+sigma_z*eps_z;
@#endif

[name = 'Productivity uncertainty shock']
log(sigma_z)=(1-rho_sigma_z)*log(sigma_zbar)+rho_sigma_z*log(sigma_z(-1))+sigma_sigma_z*eps_s;
//sigma_z=(1-rho_sigma_z)*sigma_zbar+rho_sigma_z*sigma_z(-1)+sigma_sigma_z*eps_s;

//--------------------------------------------------------------------------
//   Additional Variables
//--------------------------------------------------------------------------

// To trick Dynare into thinking that the following aren't state variables (for plotting/Andreasen method)
n_Obs = n;
z_Obs = z;
sigma_z_Obs = sigma_z;

end;

//=========================================================================
//SHOCK VARIANCES                 
//=========================================================================
shocks;
var eps_z = 1;
var eps_s = 1;      
end;

//=========================================================================
//STEADY STATE
//=========================================================================
% Provided in external steady-state file 

//--------------------------------------------------------------------------
// Check the starting values for the (det.) steady state
//--------------------------------------------------------------------------
resid;
 
//--------------------------------------------------------------------------
// Compute (det.) steady-state given the starting values
//--------------------------------------------------------------------------
steady;
 
//--------------------------------------------------------------------------
// Check Blanchard-Kahn-conditions
//--------------------------------------------------------------------------
check;

//=========================================================================
// POLICY FUNCTIONS
//=========================================================================
@#if OptionMoments > 0
stoch_simul(order=3,periods = 100000, drop = 50000,pruning,k_order_solver,irf=0,hp_filter=1600);
@#else
stoch_simul(order=3,pruning,k_order_solver,irf=0); 
@#endif

// Store names
vNames = M_.endo_names; 

//=========================================================================
// IRFs
//=========================================================================

//--------------------------------------------------------------------------
// Pure Uncertainty Shock IRF
//--------------------------------------------------------------------------
% Specifications
IrfPeriods = 20;                
BurninPeriods = 10000;   % periods for convergence from DSS to EMAS
 
% Compute EMAS
mEpsZero = zeros(BurninPeriods+IrfPeriods,M_.exo_nbr);         % shocks set to 0 to simulate without uncertainty
mIrfZero = simult_(oo_.dr.ys,oo_.dr,mEpsZero,options_.order)'; % simulate series
vEMAS = mIrfZero(1+BurninPeriods,:);                           % EMAS is any of the final points after burn-in
 
% Now simulate with an added impulse after the burn-in period
mEps = zeros(BurninPeriods+IrfPeriods,M_.exo_nbr);
mEps(1+BurninPeriods,strmatch('eps_s',M_.exo_names,'exact')) = 1; % select productivity uncertainty

% Simulate model: starting from EMWS, burn-in to get us to EMAS, then impulse
mIrf = simult_(oo_.dr.ys,oo_.dr,mEps,options_.order)';
 
% Compute proportional deviation from EMAS
mIrfProp = (mIrf(1+BurninPeriods+1:1+BurninPeriods+IrfPeriods,:)-mIrfZero(1+BurninPeriods+1:1+BurninPeriods+IrfPeriods,:))./repmat(vEMAS,IrfPeriods,1); % only valid for variables not yet logged
 
% Adjust values that are virtually zero (not exactly due to use of a solver) to exactly zero to avoid confusing graphs
mIrfProp(abs(mIrfProp)<1e-12)=0;
 
mIRFProp_zUncertainty_EMAS = mIrfProp;

//--------------------------------------------------------------------------
// GIRFs a la Andreasen et al. (2018)
//--------------------------------------------------------------------------
tic
addpath('simAndMoments3order');

// Default step: a first-order approximation needed for RunDynarePruning
stoch_simul(order = 1, noprint, nomoments, irf = 0);
f_11 = [oo_.dr.ghx oo_.dr.ghu];
vNames = M_.endo_names;
vNames_cell = cellstr(vNames);
 
// Standard Dynare command 
stoch_simul(order = 3, noprint, nomoments, irf = 0);
 
// Options for running the pruning codes
optPruning.numSim         = 2000;                   
optPruning.seedNum        = 1; 
optPruning.orderApp       = options_.order;                                     
optPruning.computeIRF     = 1;
optPruning.ySelect     = vNames;
optPruning.plotIRF = 0;
 
% Get a positive (+1) shock
optPruning.shockSize = +1;
outDynare = RunDynarePruning(optPruning,oo_,M_,f_11);
vNames_GIRF = outDynare.label_y;
vEMWS_GIRF = outDynare.Mean_y';

% Recover IRFs and scale by EMAS
Order = 3;

mIRFProp_GIRF_z_EMAS = outDynare.IRFy(:,:,1,1)'./repmat(vEMAS,IrfPeriods,1); % impose that look at 1st order
mIRFProp_GIRF_zVol_EMAS = outDynare.IRFy(:,:,3,2)'./repmat(vEMAS,IrfPeriods,1);

% Careful: the above only loads the IRFs for (what Dynare takes to be) controls. If you want to look at states
% you either have to declare them separately (as I've done e.g. for z, sigma_z, n) or use the following
% mIRFProp_GIRF_zVol_state = outDynare.IRFv(:,:,3,2)';

disp(['Run time for computation of IRFs a la Andreasen was ',num2str(toc),' seconds']);

//=========================================================================
// SAVE STUFF
//=========================================================================

save('.\Output\IRFs\IRFs_Spc', 'vEMAS', 'vNames', 'vNames_GIRF', 'mIRFProp_zUncertainty_EMAS', 'mIRFProp_GIRF_zVol_EMAS') 
