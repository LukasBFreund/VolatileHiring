%==========================================================================
% Parent file: SaM model with option value due to het. prod. draws
% *and* no entrepreneur death upon separation
%
% This file: compute steady-state elasticity of hiring rate wrt agg prod.
%
% This version: assumes a~U(aL,aH), z+a
%
% Structure: 
%           Inputs: 
%                   - x: 2x1 vector of omega and chi guesses
%                   - sPar:   structure of parameter values
%
%           Functions called:
%                  - none (unlike with normal)                 
%
%           Output:  
%                   - elasticity: 1x1 scalar elasticity
%
% Last updated: June 2021
%==========================================================================

function elasticity = fn_Elasticity_nodeath(x,sPar)

    omega = x(1);
    chi = x(2);

    aHat  = sPar.ss.aHat;
    p      = sPar.ss.p;
    aStar = sPar.ss.aStar;
    h     = sPar.ss.h;
    n     = sPar.ss.n;
    v     = sPar.ss.v;
    lnh   = log(h);
    lambdaastar = (1-omega)*(sPar.xss*(sPar.ss.z+aStar) - chi);
    lambdaahat  = (1-omega)*(sPar.xss*(sPar.ss.z+aHat)  - chi);

    syms Qh Qv Qn Qahat Qastar Qprob  Qlnh Qlambdaastar Qlambdaahat Qavez 

    % Set up the equations

% Entrepreneur death
    
%     eq   = [ -Qlambdaastar  + (1-omega)*(sPar.xss*(Qastar+Qavez) - chi);  
%     -Qlambdaahat   + (1-omega)*(sPar.xss*(Qahat+Qavez) - chi); 
%     -Qastar        + (Qahat +  sPar.aMean+sPar.uni)/2;   
%     -Qprob         + 1-((Qahat - (sPar.aMean-sPar.uni))/(2*sPar.uni));
%     -Qh            + sPar.psi*((1-(1-sPar.delta)*Qn)/Qv)^sPar.alpha;
%     -Qv            + Qprob*(sPar.Upsilon - (1-sPar.delta)*Qn);
%     -Qn            + Qh*Qv/sPar.delta;         
%     -sPar.kappa/Qh + Qlambdaahat/(1-sPar.beta*(1-sPar.delta))...
%       - ((sPar.beta*Qprob*Qh*Qlambdaastar)/((1-sPar.beta*(1-sPar.delta))*(1-sPar.beta*(1-Qprob*Qh)))...
%       - (sPar.beta*Qprob*sPar.kappa)/(1-sPar.beta*(1-Qprob*Qh)));
%     -Qlnh          + log(Qh);
%     ];

% No entrepreneur death
    eq   = [ -Qlambdaastar  + (1-omega)*(sPar.xss*(Qastar+Qavez) - chi);  
    -Qlambdaahat   + (1-omega)*(sPar.xss*(Qahat+Qavez) - chi); 
    -Qastar        + (Qahat +  sPar.aMean+sPar.uni)/2;   
    -Qprob         + 1-((Qahat - (sPar.aMean-sPar.uni))/(2*sPar.uni));
    -Qh            + sPar.psi*((1-(1-sPar.delta)*Qn)/Qv)^sPar.alpha;
    -Qv            + Qprob*(sPar.Upsilon - (1-sPar.delta)*Qn);
    -Qn            + Qh*Qv/sPar.delta;         
    -sPar.kappa/Qh  + (Qlambdaastar + sPar.beta*(1-sPar.delta)*Qprob*sPar.kappa) ...
                / (1-sPar.beta*(1-sPar.delta)*(1-Qprob*Qh))                 ...
                - (Qlambdaastar - Qlambdaahat)                            ...
                / (1-sPar.beta*(1-sPar.delta));
    -Qlnh          + log(Qh);
    ];


    % calculate elasticity
    Jacob = jacobian(eq,[Qh   Qv Qn Qahat Qastar Qprob Qlnh     Qlambdaastar Qlambdaahat Qavez]);
    Jacob = subs(Jacob, [Qh   Qv Qn Qahat Qastar Qprob Qlnh     Qlambdaastar Qlambdaahat Qavez], ...
    [ h    v  n  aHat  aStar  p  lnh      lambdaastar  lambdaahat sPar.ss.z]);
    Jacob  = double(Jacob);
    JacobY = Jacob(:,1:9);
    Jacobx = Jacob(:,10);
    output = -JacobY\Jacobx;
    elasticity = output(7,1)*sPar.ss.z; 
   
end