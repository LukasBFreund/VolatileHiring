//=========================================================================
// SaM model with option-value effect due to heterogeneous productivity 
// draws and finite mass of potential entrepeneurs.
// Version with additive production function, entrepreneur 'death'
// Holds and runs only for sigma_a > 0
//=========================================================================

//=========================================================================
//DECLARATION OF ENDOGENOUS VARIABLES
//=========================================================================
var
    m       ${m}$                    (long_name='Matches')
    n       ${n}$                    (long_name='Employment')
    u       ${n}$                    (long_name='Unemployment')
    us      ${u^s}$                  (long_name='Number of searching workers')
    v       ${v}$                    (long_name='Number of vacancies')
    h       ${h}$                    (long_name='Vacancy filling rate')
    z       ${z}$                    (long_name='Aggregate productivity')
    f       ${f}$                    (long_name='Job finding rate')
    y       ${y}$                    (long_name='Output')
    sigma_z ${\sigma_z}$             (long_name='Aggregate productivity shock volatility')
    p       ${p}$                    (long_name='Entry probability')
    aHat    ${\hat{a}}$              (long_name='Cutoff')
    aStar   ${a^*}$                  (long_name='Cond. expected productivity draw')
    JaHat   ${J(\hat{a})}$           (long_name='Match value at cutoff')
    JStar   ${J^*(a)}$               (long_name='Conditional expected match value')
    JU      ${J^U}$                  (long_name='Value of unmatched entrepreneur')
    JaHat_aHatSS ${\bar{J}(\hat{a})}$ (long_name='Match value at steady-state cutoff')

// The following are auxililiary variables
    JaHatPTilde
    JStarPTilde
    JU_d
    JaHatPTilde_aHatSS 
    n_Obs 
    z_Obs 
    sigma_z_Obs
;
 
//=========================================================================
//DECLARATION OF EXOGENOUS VARIABLES
//=========================================================================
varexo
    eps_z   ${\epsilon_z}$          (long_name='Aggregate productivity innovation')
    eps_s   ${\epsilon_{\sigma_z}}$ (long_name='Aggregate productivity uncertainty innovation')
;
 
//=========================================================================
//DECLARATION OF PARAMETERS       
//=========================================================================
parameters
    betta         ${\beta}$               (long_name='Subjective discount factor')
    eta           ${\eta}$                (long_name='Elasticity of substitution b/w differentiated products')
    alpha         ${\alpha}$              (long_name='Elasticity of matching w.r.t. u')
    delta         ${\delta}$              (long_name='Exogenous job destruction rate')      
    omega         ${\omega}$              (long_name='Bargaining weight for workers')
    kappa         %{\kappa}$              (long_name='Fixed vacancy posting cost')
    chi           ${\chi}$                (long_name='Strike value')
    psi           ${\mu}$                 (long_name='Matching efficiency')
    xss           ${\bar{x}}$             (long_name='Inverse markup')
    zBar          ${\bar{z}}$             (long_name='Aggregate productivity, steady-state value')
    rho_z         ${\rho_z}$              (long_name='Aggregate productivity shock, persistence')
    rho_sigma_z   ${\rho_{\sigma_z}}$     (long_name='Aggregate productivity uncertainty shock, persistence')
    sigma_zbar    ${\bar{\sigma}_z}$      (long_name='Aggregate productivity uncertainty shock, mean value')
    sigma_sigma_z ${\sigma_{\sigma_z}}$   (long_name='Aggregate productivity uncertainty shock, sd of innovation')

    sigma_a       ${\sigma_a}$             (long_name='Idiosyncratic productivity shock, sd')
    Upsilon       ${\Upsilon}$             (long_name='Fixed mass of entrepreneurs')
    aMean         ${\mu_a}$                (long_name='Idiosyncratic productivity shock, mean')
    aH            ${a_L}$                  (long_name='Idiosyncratic productivity shock, upper limit')
    aL            ${a^H}$                  (long_name='Idiosyncratic productivity shock, lowerl limit')
;

//=========================================================================
//PARAMETER VALUES                
//=========================================================================

load inputParameters

set_param_value('zBar', sPar.zBar);
set_param_value('betta',sPar.beta);
set_param_value('eta',sPar.eta);
set_param_value('alpha',sPar.alpha);
set_param_value('delta',sPar.delta);
set_param_value('kappa',sPar.kappa);
set_param_value('psi',sPar.psi);
set_param_value('xss',sPar.xss);

% Entrepreneurs
set_param_value('sigma_a',sPar.sigma_a);
set_param_value('Upsilon',sPar.Upsilon);
set_param_value('aMean',sPar.aMean);
set_param_value('aH',sPar.aH);
set_param_value('aL',sPar.aL);

% Wage
set_param_value('omega',sPar.omega);
set_param_value('chi',sPar.chi);

% Stochastic processes
set_param_value('rho_z',sPar.rho_z);
set_param_value('sigma_zbar',sPar.sigma_zbar);
set_param_value('sigma_sigma_z',sPar.sigma_sigma_z);
set_param_value('rho_sigma_z',sPar.rho_sigma_z);

//=========================================================================
//MODEL EQUATIONS                 
//=========================================================================
model;

[name = 'Searching workers']
us = 1-(1-delta)*n(-1);

[name = 'Matching function']
m = psi*us^alpha*v^(1-alpha);

[name = 'Job finding rate ']
f= m/us;

[name = 'Vacancy filling rate']
h = m/v;
%h = STEADY_STATE(h);

[name = 'Employment law of motion']
n = (1-delta)*n(-1)+m;        

[name = 'Unemployment rate']
u = 1-n;

[name = 'Match value at cutoff']
JaHat = ((1-omega)*(xss*aHat))/(1-betta*(1-delta)) + JaHatPTilde;
//JaHatPTilde = (1-omega)*(xss*z-chi)+betta*((1-delta)*JaHatPTilde(+1)); % entrepreneur death
JaHatPTilde = (1-omega)*(xss*z-chi)+betta*((1-delta)*JaHatPTilde(+1)+delta*JU(+1));

[name = 'Expected value of newly matched firm']
JStar = ((1-omega)*(xss*aStar))/(1-betta*(1-delta)) + JStarPTilde;
JStarPTilde = (1-omega)*(xss*z-chi)+betta*((1-delta)*JStarPTilde(+1)+delta*JU(+1));

[name = 'Value of unmatched entrepreneur']
JU = betta*JU(+1)-p*kappa+p*h*(JStar-betta*JU(+1));

[name = 'Indifference condition']
kappa = h*(JaHat-betta*JU(+1));

[name = 'Definition of p']
p = 1 - (aHat-aL)/(aH-aL);  %1-F(aHat)

[name = 'Definition of a-Star']
aStar = 0.5*(aH+aHat);  % E(a|a>aHat)

[name = 'Vacancy creation']
v = p*(Upsilon-(1-delta)*n(-1));

//--------------------------------------------------------------------------
//   Stochastic processes
//--------------------------------------------------------------------------

[name = 'Aggregate productivity']
z = (1-rho_z)* STEADY_STATE(z)+rho_z*z(-1)+sigma_z(-1)*eps_z;

[name = 'Aggregate productivity volatility']
//log(sigma_z)=(1-rho_sigma_z)*log(sigma_zbar)+rho_sigma_z*log(sigma_z(-1))+sigma_sigma_z*eps_s;
sigma_z=(1-rho_sigma_z)*sigma_zbar+rho_sigma_z*sigma_z(-1)+sigma_sigma_z*eps_s;

//--------------------------------------------------------------------------
//   Additional Variables
//--------------------------------------------------------------------------

[name = 'Output']
y = n*(z+aStar);

[name = 'Value of matched entrepreneur at steady state cut-off level']
JaHat_aHatSS = ((1-omega)*(xss*STEADY_STATE(aHat)))/(1-betta*(1-delta)) + JaHatPTilde_aHatSS;
JaHatPTilde_aHatSS = (1-omega)*(xss*z-chi)+betta*((1-delta)*JaHatPTilde_aHatSS(+1)+delta*JU(+1));

[name = 'Value of entrepreneur leaving without match']
JU_d = betta*JU(+1);

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
steady(solve_algo=4);
 
//--------------------------------------------------------------------------
// Check Blanchard-Kahn-conditions
//--------------------------------------------------------------------------
check;
 
//=========================================================================
// POLICY FUNCTIONS
//=========================================================================
//@#if OptionMoments > 0
//stoch_simul(order=3,periods = 100000, drop = 50000,pruning,k_order_solver,irf=0,hp_filter=1600);
//@#else
stoch_simul(order=3,pruning,k_order_solver,irf=0); 
//@#endif

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
mEps(1+BurninPeriods,strmatch('eps_s',M_.exo_names,'exact')) = 1; % select uncertainty

% Simulate model: starting from EMWS, burn-in to get us to EMAS, then impulse
mIrf = simult_(oo_.dr.ys,oo_.dr,mEps,options_.order)';
 
% Compute proportional deviation from EMAS
mIrfProp = (mIrf(1+BurninPeriods+1:1+BurninPeriods+IrfPeriods,:)-mIrfZero(1+BurninPeriods+1:1+BurninPeriods+IrfPeriods,:))./repmat(vEMAS,IrfPeriods,1); % only valid for variables not yet logged
 
% Adjust values that are virtually zero (not exactly due to use of a solver) to exactly zero to avoid confusing graphs
mIrfProp(abs(mIrfProp)<1e-12)=0;
 
mIRFProp_zUncertainty_EMAS = mIrfProp;

//--------------------------------------------------------------------------
// Standard productivity shock IRF
//--------------------------------------------------------------------------

% Method 1: manual simulation (at some order k)
stoch_simul(order = 1, noprint, nomoments, irf = 0);
mEps = zeros(BurninPeriods+IrfPeriods,M_.exo_nbr);
mEps(1+BurninPeriods,strmatch('eps_z',M_.exo_names,'exact')) = 1; % select productivity shock
mIrf = simult_(oo_.dr.ys,oo_.dr,mEps,options_.order)';
mIrfProp = (mIrf(1+BurninPeriods+1:1+BurninPeriods+IrfPeriods,:)-mIrfZero(1+BurninPeriods+1:1+BurninPeriods+IrfPeriods,:))./repmat(vEMAS,IrfPeriods,1); % only valid for variables not yet logged
mIrfProp(abs(mIrfProp)<1e-12)=0;
mIRFProp_Order1_z_EMAS = mIrfProp;

% Method 2: dynare internal (yields virtually identical results to method 1)
%stoch_simul(order=1,noprint, nomoments,irf=20,nograph); 
%sIRFDynare = oo_.irfs;

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
vNames_Andreasen = outDynare.label_y;
 vEMWS_Andreasen = outDynare.Mean_y';

mIRFProp_Andreasen_z_EMAS = outDynare.IRFy(:,:,options_.order,1)'./repmat(vEMAS,IrfPeriods,1);
mIRFProp_Andreasen_Order1_z_EMAS = outDynare.IRFy(:,:,1,1)'./repmat(vEMAS,IrfPeriods,1);
mIRFProp_Andreasen_Order3_z_EMAS = outDynare.IRFy(:,:,3,1)'./repmat(vEMAS,IrfPeriods,1);
mIRFProp_Andreasen_zVol_EMAS = outDynare.IRFy(:,:,options_.order,2)'./repmat(vEMAS,IrfPeriods,1);

% Careful: the above only loads the IRFs for (what Dynare takes to be) controls. If you want to look at states
% you either have to declare them separately (as I've done e.g. for z, sigma_z, n) or use the following
% mIRFProp_Andreasen_zVol_state = outDynare.IRFv(:,:,3,2)';

disp(['Run time for computation of IRFs a la Andreasen was ',num2str(toc),' seconds']);

//=========================================================================
// SAVE STUFF
//=========================================================================

save('Output\IRFs\IRFs_HetFirm_Gamma0Model_Spc', 'vEMAS', 'vNames', 'vNames_Andreasen', 'mIRFProp_zUncertainty_EMAS', 'mIRFProp_Andreasen_zVol_EMAS','mIRFProp_Andreasen_z_EMAS','mIRFProp_Andreasen_Order1_z_EMAS','mIRFProp_Order1_z_EMAS','mIRFProp_Andreasen_Order3_z_EMAS'); %,'sIRFDynare')

