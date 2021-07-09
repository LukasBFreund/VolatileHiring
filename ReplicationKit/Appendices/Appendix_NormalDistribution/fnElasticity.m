%==========================================================================
% Parent file: SaM model with option value due to het. prod. draws.
% This file: compute steady-state elasticity of hiring rate wrt agg prod.
%
% This version: a~Normal; numerical approach (validated for uniform dist.)
%
% Case where y = z+a and entrepreneurs have 0 value upon separation
%
% Structure: 
%           Inputs: 
%                   - x: 2x1 vector of omega and chi guesses
%                   - sPar:   structure of parameter values
%
%           Functions called:
%                  - fn_Cutoff                
%
%           Output:  
%                   - elasticity: 1x1 scalar elasticity
%
% Last updated: 4 November 2020
%==========================================================================

function elasticity = fnElasticity(x,sPar)
    
    sParGuess = sPar;

    % Impose conjectured values of omega and chi 
    sParGuess.omega = x(1);
    sParGuess.chi = x(2);

    % Set up *numerical* computation of elasticity 
    difz     = 0.0000001;
    %dlnz     = log((sParGuess.ss.z+difz)/sParGuess.ss.z);
  
    % Explicitly compute h before change in z (b/c influenced by para change)
    vInit = [0.00,0.7];  
    options = optimoptions('fsolve','Display','none',...
                'MaxFunctionEvaluations',100000,'StepTolerance',1e-16,'OptimalityTolerance',1e-16);
    
     vSol = fsolve(@(y) fn_Cutoff(y,sParGuess),vInit,options);
     vSol = real(vSol);
     h0 = vSol(2);
            
     % Compute h after change in z   
     sParGuess.ss.z = sParGuess.ss.z + difz;
     vSol = fsolve(@(y) fn_Cutoff(y,sParGuess),vInit,options);
     vSol = real(vSol);
     h = vSol(2);
     elasticity = ((h-h0)/h0)/((difz)/sParGuess.ss.z);

end
