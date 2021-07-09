%==========================================================================
% Parent file: SaM model with option value due to het. prod. draws,
% gamma version
%
% Structure: 
%           Inputs: 
%                   - vGuess: 2x1 vector of omega and chi guesses
%                   - sPar:   structure of parameter values
%                   - target:  target value of ss-elasticity of h w.r.t. z
%
%           Functions called:
%                  - fn_Elasticity_nodeath: computes ss-elasticity of h w.r.t. z   
%                    (adjusted s.t. no entrepreneur death)
%
%           Output:  
%                   - vError: 2x1 vector of errors
%
% Last updated: June 2021
%==========================================================================


function vError = fn_omega_chi_gamma(vGuess,sPar,target)

    % Unpack guess and some parameters to ease notation
    omega = vGuess(1,1);
    chi   = vGuess(1,2);

    aHat  = sPar.ss.aHat;
    p      = sPar.ss.p;
    aStar = sPar.ss.aStar;
    h     = sPar.ss.h;
%     n     = sPar.ss.n;
%     v     = sPar.ss.v;
%     lnh   = log(h);
%     lambdaastar = (1-omega)*(sPar.xss*(sPar.ss.z+aStar) - chi);
%     lambdaahat  = (1-omega)*(sPar.xss*(sPar.ss.z+aHat)  - chi);

   % Compute error 
    LHS   = sPar.kappa/h;
  
    eqn = @(J) [-J(1) + sPar.beta*J(1) - p*sPar.kappa+p*h*(J(3)-sPar.beta*J(1));  
            -J(2) + (1-sPar.gamma)*J(1) + sPar.gamma*sPar.beta*J(2);
            -J(3) + (1-omega)*(sPar.xss*(aStar+sPar.ss.z)-chi)+sPar.beta*((1-sPar.delta)*J(3)+sPar.delta*J(2))];

    J0 = [0.1, 0.1, 0.3];
    %options = optimoptions('fsolve','Display','none');
    options = optimoptions('fsolve','Display','none','Algorithm','trust-region',...
        'MaxFunctionEvaluations',100000,'StepTolerance',1e-12,'OptimalityTolerance',1e-12);
    [JSol,~] = fsolve(eqn,J0,options);
    JU = JSol(1);
    JUtilde = JSol(2);
    JaHat = ((1-omega)*(sPar.xss*(aHat+sPar.ss.z)-chi)+sPar.beta*sPar.delta*JUtilde)/(1-sPar.beta*(1-sPar.delta));

    vError(1) = -LHS + JaHat - sPar.beta*JU; % careful, this is a different JUtilde than before


    % Compute elasticity
    elasticity = fn_Elasticity_gamma([omega, chi],sPar);  
    
    vError(2) = abs(target-elasticity); 

    %disp('omega chi elasticity target'); 
    %disp([omega chi elasticity target])
    %disp([vError(1) vError(2)]);
end
