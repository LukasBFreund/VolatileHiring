function [ys,check] = dynareFRUncertainty_steadystate(ys,exo)

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
    z = zss;
    sigma_z = sigma_zbar;
    h = hss; 
    u = uss;
    n = nss;
    m = mss;
    us = usss;
    f = fss;
    v = vss;
    theta = v/us;
    y = yss;
    c = css;
    x = xss;
    J = Jss;
    wL = wss;      
    wN = wss;  
    w = wss;     
    n_Obs = n;
    z_Obs = z;
    sigma_z_Obs = sigma_z;
    
    lambda = css^(-xi);
    RReal = 1/betta;
    
    EJ1 = J;
    ELambdaratio = 1;
    CovJLambdaratio = 0;
    JCheck = J;
    JNoCov = J;
    JNoSDF = J;
    JNoNonlin = J;
  
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