function [ys,check] = dynareSAMOptionValue_NoDispersion_steadystate(ys,exo)

global M_ 
% read out parameters to access them with their name
NumberOfParameters = M_.param_nbr;
for ii = 1:NumberOfParameters
  paramname = deblank(M_.param_names(ii,:));
  eval([ paramname ' = M_.params(' int2str(ii) ');']);
end
% initialize indicator
check = 0;

%% Enter model equations here
    load inputParameters;  
    m           = sPar.ss.m;
    n           = sPar.ss.n;
    u           = 1-n;
    us          = sPar.ss.us;
    f           = m/us;
    v       	= sPar.ss.v;
    h           = sPar.ss.h;
    z           = sPar.ss.z;
    sigma_z     = sPar.ss.sigma_z;
    JaHat       = sPar.ss.JaHat;
    JStar       = sPar.ss.JStar;
    JaHatPTilde = ((1-sPar.omega)*(sPar.xss*sPar.ss.z-sPar.chi))/(1-sPar.beta*(1-sPar.delta));
    JStarPTilde = JaHatPTilde;
    JU          = sPar.ss.JU;
    aHat        = 0;
    aStar       = 0;
    p           = sPar.ss.p;  
    y           = n*(z+aStar);
    JaHat_aHatSS = JaHat;
    JaHatPTilde_aHatSS = JaHatPTilde;
    JU_d        = sPar.beta*JU;
    n_Obs       = n;
    z_Obs       = z;
    sigma_z_Obs = sigma_z;  
    
%% End own model equations

for iter = 1:length(M_.params) % Update parameters set in the file
  eval([ 'M_.params(' num2str(iter) ') = ' M_.param_names(iter,:) ';' ])
end

NumberOfEndogenousVariables = M_.orig_endo_nbr; % Auxiliary variables are set automatically
for ii = 1:NumberOfEndogenousVariables
  varname = deblank(M_.endo_names(ii,:));
  eval(['ys(' int2str(ii) ') = ' varname ';']);                % Get the steady state value of this variable
end
end