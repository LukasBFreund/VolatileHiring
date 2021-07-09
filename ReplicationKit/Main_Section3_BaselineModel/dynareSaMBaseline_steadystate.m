function [ys,check] = dynareSaMBaseline_steadystate(ys,exo)

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

    z = sPar.ss.z;
    sigma_z = sPar.ss.sigma_z;
    h = sPar.hBar; 
    u = sPar.ss.u;
    n = 1-u;
    m = sPar.ss.m;
    us = sPar.ss.us;
    f = sPar.ss.f;
    v = sPar.ss.v;
    theta = sPar.ss.theta;
    y = sPar.ss.y;
    c = css; % depends on option in mod file
    J = sPar.ss.J;
    wL = sPar.ss.w;      
    wN = sPar.ss.w;  
    w = sPar.ss.w;     
    n_Obs = n;
    z_Obs = z;
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