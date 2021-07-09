%==========================================================================
% Parent file: SaM model with option value due to het. prod. draws.
% This file: computes the errors for cutoff aHat and hiring rate h from the  
% key steady-state cutoff equation, taking into account endogeneity of h.
% 
% This version: 'gamma model'

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
% Last updated: July 2021
%==========================================================================

function vError = fn_Cutoff_gamma(vGuess,sPar)
    %% Prepare
    % Unpack guess
    aHat  = vGuess(1);
    h     = vGuess(2);

    % Compute implied p and a*
    p = 1 - (aHat-sPar.aL)/(sPar.aH-sPar.aL);  %1-F(aHat)
    aStar = 0.5*(sPar.aH+aHat);  % E(a|a>aHat)

    %% Compute error for aHat
    LHS   = sPar.kappa/h;
  
    eqn = @(J) [-J(1) + sPar.beta*J(1) - p*sPar.kappa+p*h*(J(3)-sPar.beta*J(1));  
            -J(2) + (1-sPar.gamma)*J(1) + sPar.gamma*sPar.beta*J(2);
            -J(3) + (1-sPar.omega)*(sPar.xss*(aStar+sPar.ss.z)-sPar.chi)+sPar.beta*((1-sPar.delta)*J(3)+sPar.delta*J(2))];

    J0 = [0.1, 0.1, 0.3];
    options = optimoptions('fsolve','Display','none','Algorithm','trust-region',...
        'MaxFunctionEvaluations',100000,'StepTolerance',1e-12,'OptimalityTolerance',1e-12);
    [JSol,~] = fsolve(eqn,J0,options);
    JU = JSol(1);
    JUtilde = JSol(2);
    JaHat = ((1-sPar.omega)*(sPar.xss*(aHat+sPar.ss.z)-sPar.chi)+sPar.beta*sPar.delta*JUtilde)/(1-sPar.beta*(1-sPar.delta));

    Error_aHat = -LHS + JaHat - sPar.beta*JU; 

    %% Compute error for hiring rate h -- accounting for gamma
    % First compute implied vacancies and employment
    % Old equations (gamma = 0):
    % v   = prob*(upsilon - (1-delta)*n);
    % n = (h*v)/delta; % from LoM
    %
    % New equations:
    % v = 
    % n = (h*v)/delta;
    % 
    % Substituting the second into the first and solving for v yields the
    % below
    v = p*sPar.Upsilon/(1+p*(((1-sPar.delta)+(sPar.delta*sPar.gamma)/(1-sPar.gamma))*(h/sPar.delta)));  
    n = h*v/sPar.delta;

    % Error between implied (from matching fn) hiring rate h and our guess 
    hImplied = sPar.psi*((1-(1-sPar.delta)*n)/v)^sPar.alpha;
    Error_h = -h + hImplied;

    vError(1) = Error_aHat;
    vError(2) = Error_h;
end
