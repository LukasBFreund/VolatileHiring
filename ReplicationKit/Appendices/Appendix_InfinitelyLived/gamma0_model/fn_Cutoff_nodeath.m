%==========================================================================
% Parent file: SaM model with option value due to het. prod. draws.
% This file: computes the errors for cutoff aHat and hiring rate h from the  
% key steady-state cutoff equation, taking into account endogeneity of h.
% This version: assumes a~U(aL,aH), so can solve symbolically; case 
% where y = z+a...

%...but replace assumption that entrepreneurs 'die' upon separation
% with non-zero continuation value

% Structure: 
%           Inputs: 
%                   - vGuess: 2x1 vector of aHat and h guesses
%                   - sPar:   structure of parameter values
%
%           Functions called:
%                  - none (unlike with normal)                 
%
%           Output:  
%                   - vError: 2x1 vector of errors
%
% Last updated: June 2021 
%==========================================================================

function vError = fn_Cutoff_nodeath(vGuess,sPar)
    %% Prepare
    % Unpack guess
    aHat  = vGuess(1);
    h     = vGuess(2);

    % Compute implied p and a*
    p = 1 - (aHat-sPar.aL)/(sPar.aH-sPar.aL);  %1-F(aHat)
    aStar = 0.5*(sPar.aH+aHat);  % E(a|a>aHat)
    
    % Note that if we wanted to use the normal distribution...
    % p  = 1-normcdf(aHat,sPar.aMean,sPar.sigma_a);
    % normpdfEval = ((1/sqrt(2*pi))*exp(-0.5*((aHat-sPar.aMean)/sPar.sigma_a).^2)); 
    % aStar = sPar.aMean+sPar.sigma_a*normpdfEval/p;

    %% Compute error for aHat
    LHS   = sPar.kappa/h;
    
    % Version with entrepreneur death
    
%     LambdaaStar =  ((1-sPar.omega)*(sPar.xss*(aStar+sPar.ss.z)-sPar.chi));
%     LambdaaHat =  ((1-sPar.omega)*(sPar.xss*(aHat+sPar.ss.z)-sPar.chi));
%     JaStar =   LambdaaStar/(1-sPar.beta*(1-sPar.delta));
%     JU  = (p*(h*JaStar-sPar.kappa))/(1-sPar.beta*(1-p*h));
%     JUTilde = sPar.beta*JU;
%     JaHat = LambdaaHat/(1-sPar.beta*(1-sPar.delta));
%    Error_aHat = -LHS + JaHat - JUTilde;

    % Version without entrepreneur death
    RHS1Num   = (1-sPar.omega)*(sPar.xss*(aStar+sPar.ss.z)-sPar.chi)+sPar.beta*(1-sPar.delta)*p*sPar.kappa;
    RHS1Den   = 1-sPar.beta*(1-sPar.delta)*(1-p*h);
    RHS1  = RHS1Num/RHS1Den;
    RHS2  = ((1-sPar.omega)*sPar.xss*(aStar-aHat))/(1-sPar.beta*(1-sPar.delta));    
    Error_aHat = -LHS + RHS1 - RHS2;

    % Alternative approach that should yield the same outcome: involves
    % numerical approximation but is robust to generalization
%     eqn = @(J) [-J(1) + sPar.beta*J(1) - p*sPar.kappa+p*h*(J(2)-sPar.beta*J(1));    
%     -J(2) + (1-sPar.omega)*(sPar.xss*(aStar+sPar.ss.z)-sPar.chi)+sPar.beta*((1-sPar.delta)*(J(2)-J(1))+J(1))];
%     J0 = [0.1, 0.3];
%     options = optimoptions('fsolve','Display','iter');
%     [JSol,~] = fsolve(eqn,J0,options);
%     JU = JSol(1);
%     JUTilde = sPar.beta*JU;
%     %JaStar = JSol(2);
%     JaHat = ((1-sPar.omega)*(sPar.xss*(aHat+sPar.ss.z)-sPar.chi)+sPar.beta*sPar.delta*JU)/(1-sPar.beta*(1-sPar.delta));
%     Error_aHat = -LHS + JaHat - JUTilde;

    
    %% Compute error for hiring rate h, 
    % first compute implied vacancies and employment, using 
    % v   = prob*(upsilon - (1-delta)*n);
    % n = (h*v)/delta;
    % and substitute the 2nd into the first to solve for v

    v = p*sPar.Upsilon/(1+p*(1-sPar.delta)*h/sPar.delta);
    n = h*v/sPar.delta;

    % Error between implied (from definition) hiring rate h and our guess 
    hImplied = sPar.psi*((1-(1-sPar.delta)*n)/v)^sPar.alpha;
    Error_h = -h + hImplied;

    vError(1) = Error_aHat;
    vError(2) = Error_h;
end
