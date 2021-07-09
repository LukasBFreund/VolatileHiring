%==========================================================================
% Parent file: SaM model with option value due to het. prod. draws.
% This file: solves for ahat and chi that satisfies cut-off equation and elaticity
%
% This version: assumes a~U(aL,aH); z+a and entrepreneurs have 0 value upon separation
%
% Structure: 
%           Inputs: 
%                   - vGuess: 2x1 vector of omega and chi guesses
%                   - sPar:   structure of parameter values
%                   - target:  target value of ss-elasticity of h w.r.t. z
%
%           Functions called:
%                  - fnElasticity: computes ss-elasticity of h w.r.t. z              
%
%           Output:  
%                   - vError: 2x1 vector of errors
%
% Last updated: 22 October 2020
%==========================================================================

function vError = fn_omega_chi(vGuess,sPar,target)

    % Unpack guess and some parameters to ease notation
    omega = vGuess(1,1);
    chi   = vGuess(1,2);

    aHat  = sPar.ss.aHat;
    p      = sPar.ss.p;
    aStar = sPar.ss.aStar;
    h     = sPar.ss.h;
    n     = sPar.ss.n;
    v     = sPar.ss.v;
    lnh   = log(h);
    lambdaastar = (1-omega)*(sPar.xss*(sPar.ss.z+aStar) - chi);
    lambdaahat  = (1-omega)*(sPar.xss*(sPar.ss.z+aHat)  - chi);

    % Compute error 

    vError(1) = -sPar.kappa/h + lambdaahat/(1-sPar.beta*(1-sPar.delta))...
        - ((sPar.beta*p*h*lambdaastar)/((1-sPar.beta*(1-sPar.delta))*(1-sPar.beta*(1-p*h)))...
        - (sPar.beta*p*sPar.kappa)/(1-sPar.beta*(1-p*h)));

    % Now let's do the fun part -> or rather, do it in fnElasticity

    elasticity = fn_Elasticity([omega, chi],sPar);
    vError(2) = abs(target-elasticity); 

    %disp('omega chi elasticity target'); 
    %disp([omega chi elasticity target])
    %disp([vError(1) vError(2)]);
end
