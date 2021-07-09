%==========================================================================
% Parent file: SaM model with option value due to het. prod. draws.
%
% Gamma version -> attempt at symbolic
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
% Last updated: November 2020
%==========================================================================

function elasticity = fn_Elasticity_gamma(x,sPar)
    
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
    s = n/(1-sPar.gamma);
    Gammma = sPar.Upsilon-(1-sPar.delta)*n - sPar.delta*sPar.gamma*s;
    
    %  Value function evaluated at ss given omega and chi
    eqn = @(J) [-J(1) + sPar.beta*J(1) - p*sPar.kappa+p*h*(J(3)-sPar.beta*J(1));  
    -J(2) + (1-sPar.gamma)*J(1) + sPar.gamma*sPar.beta*J(2);
    -J(3) + (1-omega)*(sPar.xss*(aStar+sPar.ss.z)-chi)+sPar.beta*((1-sPar.delta)*J(3)+sPar.delta*J(2))];

    J0 = [0.1, 0.1, 0.3];
    options =  optimoptions('fsolve','Display','none',...
        'MaxFunctionEvaluations',100000,'StepTolerance',1e-14,'OptimalityTolerance',1e-14);
    [JSol,~] = fsolve(eqn,J0,options);
    JU = JSol(1);
    JUtilde = JSol(2);
    JaStar = JSol(3);
    JaHat = ((1-omega)*(sPar.xss*(aHat+sPar.ss.z)-chi)+sPar.beta*sPar.delta*JUtilde)/(1-sPar.beta*(1-sPar.delta));

   
    % Symbolic system of equations
    syms Qh Qv Qn Qahat Qastar Qprob  Qlnh Qlambdaastar Qlambdaahat QJaStar QJaHat QJU QJUtilde QGamma Qs Qavez
 
    eq   = [ -Qlambdaastar  + (1-omega)*(sPar.xss*(Qastar+Qavez) - chi);  
    -Qlambdaahat   + (1-omega)*(sPar.xss*(Qahat+Qavez) - chi); 
    -Qastar        + (Qahat +  sPar.aMean+sPar.uni)/2;   
    -Qprob         + 1-((Qahat - (sPar.aMean-sPar.uni))/(2*sPar.uni));
    -Qh            + sPar.psi*((1-(1-sPar.delta)*Qn)/Qv)^sPar.alpha;
    -Qv            + Qprob*QGamma;    
    -Qn            + Qh*Qv/sPar.delta;         
    -sPar.kappa/Qh + QJaHat - sPar.beta*QJU;
    -QJaStar       + Qlambdaastar + sPar.beta*((1-sPar.delta)*QJaStar+sPar.delta*QJUtilde); 
    -QJaHat        + (Qlambdaahat +sPar.beta*sPar.delta*QJUtilde)/(1-sPar.beta*(1-sPar.delta));
    -QJU           + (Qprob*(Qh*QJaStar-sPar.kappa))/(1-sPar.beta*(1-Qprob*Qh));  
    -QJUtilde      + (1-sPar.gamma)*QJU + sPar.gamma*sPar.beta*QJUtilde;
    -QGamma        + sPar.Upsilon-(1-sPar.delta)*Qn - sPar.delta*sPar.gamma*Qs;
    -Qs            + Qn/(1-sPar.gamma);
    -Qlnh          + log(Qh);      
      ];
    
    % Calculate elasticity
    Jacob = jacobian(eq,[Qh Qv Qn Qahat Qastar Qprob  Qlnh Qlambdaastar Qlambdaahat QJaStar QJaHat QJU QJUtilde QGamma Qs Qavez]);
    Jacob = subs(Jacob, [Qh Qv Qn Qahat Qastar Qprob  Qlnh Qlambdaastar Qlambdaahat QJaStar QJaHat QJU QJUtilde QGamma Qs Qavez], ...
    [ h    v  n  aHat  aStar  p  lnh      lambdaastar  lambdaahat JaStar JaHat JU JUtilde Gammma s sPar.ss.z]);
    Jacob  = double(Jacob);
    JacobY = Jacob(:,1:15); % symmetric
    Jacobx = Jacob(:,16); % z is last
    output = -JacobY\Jacobx;
    elasticity = output(7,1)*sPar.ss.z; % logh is in 7thr
   
end
