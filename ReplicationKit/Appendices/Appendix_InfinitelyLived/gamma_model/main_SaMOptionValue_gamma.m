%%=========================================================================
% Replication code for den Haan, Freund and Rendahl (2021)
%
% Infinitely lived firms, but with stochastic draws for separated 
% entrepreneurs: with probability gamma, they remain idle for one
% period after separation (i.e., 1-gammma Poisson arrival rate of new draw).
% 
% Where does gamma > 0 lead to adjustments vs. gamma = 0 version,
% specifically when allowing for recalibration?
% - cutoff function
% - elasticity function
% - steady-state computation
%
% Run on Matlab R2019b, Dynare 4.4.3
% Last updated: July 2021
%
% For any questions, please email lukas.beat.freund@gmail.com
%%=========================================================================

%% Housekeeping
%--------------------------------------------------------------------------
clear; close all; clc;
TimeStart = tic;
load inputSettings

%% Key user choices
%--------------------------------------------------------------------------

% Choose whether to apply recalibration procedure or not 
optionRecalibration = 'yes';        % 'yes' or 'no'
optionRecalibration_type = 'standard'; % 'standard' or 'fixed' (imposes omega and chi from death-baseline with sigma_a = 0.003)

% Set probability of separated  getting new draw right away
sPar.gamma = 0;          
p_gamma = sPar.gamma;   

% Idiosyncratic dipersion
sPar.sigma_a = 0.025;  

% Choose how many different models (or IRFs) to show
sPlotting.options.IRFs = 'Pure uncertainty'; % choose b/w 'Pure uncertainty', 'Pure uncertainty and total volatility'

% Plotting
OptionPrint = 0;        
TargetPath ='Output\Figures\';
FigName = 'fig_HetFirm_GammaModel_sigma0025_recalib_gamma0';
sPlotting.IrfPeriodsNum = 10;
sPlotting.numSubplotV = 4;
sPlotting.numSubplotH = 2;
sPlotting.variables = {'JaHat_aHatSS','JU_d','p','aStar','f','h','y','u'}; 
sPlotting.names = {'Value of entrepreneur at steady-state cutoff level','Value of unmatched entrepreneur',...
         'Entry probability','Cond. expected. productivity draw','Job finding rate',...
         'Hiring rate','Output','Unemployment rate'};
     
     
% 2x2
% sPlotting.numSubplotV = 2;
% sPlotting.numSubplotH = 2;
% sPlotting.variables = {'JaHat_aHatSS','JU_d','p','u'}; 
% sPlotting.names = {'Value of entrepreneur at steady-state cutoff level','Value of unmatched entrepreneur',...
%          'Entry probability','Unemployment rate'};
     
%% Parameterization 
%---------------------------------------------------------------------------
% Shared parameters
pTarget = 0.5;                      % targeted ss entry-probability

% Usual parameters
sPar.zBar = 1;   
sPar.aMean = 0; 	
sPar.hBar = 0.7;    
sPar.uBar = 0.064;  
sPar.beta = 0.99;   
sPar.eta = 10;    
sPar.xss = (sPar.eta-1)/sPar.eta;
sPar.alpha = 0.5; 
sPar.delta = 0.1;   
sPar.TZElasticity = 7.0514; 
target = -sPar.alpha*sPar.TZElasticity; % we're working with the ss-elasticity of h w.r.t. z
sPar.rho_z = 0.95;           
sPar.sigma_zbar = 0.01;    
sPar.rho_sigma_z = 0.76;       
sPar.sigma_sigma_z = 0.392*sPar.sigma_zbar; % level specification
sPar.ss.sigma_z = sPar.sigma_zbar;
  
sPar.uni = sPar.sigma_a/(sqrt(1/3));
sPar.aH = sPar.aMean + sPar.uni;
sPar.aL = sPar.aMean - sPar.uni;  

switch optionRecalibration 
    case 'no'
        
        sPar.ss.z = sPar.zBar;
        hss = sPar.hBar;  
        nss = 1-sPar.uBar;
        mss = sPar.delta*nss;
        usss = 1-(1-sPar.delta)*nss;
        fss = mss/usss;
        vss = mss/hss;
        thetass= vss/usss;
        sPar.psi = mss/(usss^sPar.alpha*vss^(1-sPar.alpha));
        sPar.kappa = 0.14;
        sPar.wss = sPar.xss*sPar.zBar-(1-sPar.beta*(1-sPar.delta))*(sPar.kappa/hss);  
        FS = (sPar.xss*1)/(sPar.alpha*sPar.TZElasticity);
        sPar.chi = sPar.xss*sPar.zBar-FS; 
        sPar.omega = (sPar.wss-sPar.chi)/(sPar.xss*sPar.zBar-sPar.chi); 
  
        % Finite mass of firms: set such that no discontinuity when moving from zero to non-zero sigma_a
        sPar.Pss = pTarget;
        % sPar.Upsilon = vss/sPar.Pss + (1-sPar.delta)*nss;
        sPar.ss.s =  nss/(1-sPar.gamma);
        sPar.Upsilon = vss/sPar.Pss + (1-sPar.delta)*nss + sPar.delta*sPar.gamma*sPar.ss.s; 
        
    case 'yes' 
        
        % We can use the following values, e.g. sPar.hBar for sPar.ss.h b/c
        % our recalibration ensures that we're hitting these targets even if
        % sigma_a>0
        sPar.ss.p = pTarget;   
        sPar.ss.h = sPar.hBar;
        sPar.ss.n = 1-sPar.uBar;
        sPar.ss.m = sPar.delta*sPar.ss.n;
        sPar.ss.us = 1-(1-sPar.delta)*sPar.ss.n;
        sPar.ss.f = sPar.ss.m/sPar.ss.us;
        sPar.ss.v = sPar.ss.m/sPar.ss.h;
        sPar.ss.theta = sPar.ss.v/sPar.ss.us;
        sPar.psi = sPar.ss.m/(sPar.ss.us^sPar.alpha*sPar.ss.v^(1-sPar.alpha));          
        
        % Given our target for p, we get aHat and aStar analytically given
        % uniform distribution formulas
        sPar.ss.aHat = sPar.aMean + sPar.uni - (sPar.ss.p*2*sPar.uni);   % p = 1-F(aHat)       
        sPar.ss.aStar = 0.5*(sPar.aH+sPar.ss.aHat);                      % a* = E[a|a=>aHat)
        
        % We recalibrate aggregate productivity such that mean total
        % productivity is still equal to unity...
        sPar.ss.z = 1-sPar.ss.aStar;     
        
        % ...which is neat because that way vacancy posting costs will be
        % unchanged
        sPar.kappa =  sPar.ss.n*(sPar.ss.z+sPar.ss.aStar)*0.02/sPar.ss.v; 
        
        % Mass of firms follows from targeted vacancies and targeted entry
        % probability...
        % ... and in this version is also affected by gamma!
        sPar.ss.s =  sPar.ss.n/(1-sPar.gamma);
        sPar.Upsilon = sPar.ss.v/sPar.ss.p + (1-sPar.delta)*sPar.ss.n + sPar.delta*sPar.gamma*sPar.ss.s; 

        % We still have to solve for omega and chi to hit the
        % ... targeted elasticity w.r.t. z (target value comes from free-entry model)

        % To initialize omega and chi, use values from sigma_a = 0
        w0 = sPar.xss*(sPar.aMean+sPar.ss.z)-(1-sPar.beta*(1-sPar.delta))*(sPar.kappa/sPar.ss.h);  
        FS = (sPar.xss*(sPar.aMean+sPar.ss.z))/(sPar.alpha*sPar.TZElasticity);
        chi0 = sPar.xss*sPar.ss.z-FS; 
        omega0 = (w0-chi0)/(sPar.xss*sPar.ss.z-chi0); 

        % Allow for different methods of recalibration
        switch optionRecalibration_type
            case 'standard'              
                % Solve for targeted omega and chi
                % NB: the function called below is adjusted for the 'gamma' formulation
                vInitOmegaChi = [omega0 chi0];                
                options = optimoptions('fsolve','Display','iter','Algorithm','trust-region',...
                'MaxFunctionEvaluations',100000,'StepTolerance',1e-6,'OptimalityTolerance',1e-6);
                [vSolOmegaChi0,~,~] = fsolve(@(x) fn_omega_chi_gamma(x,sPar,target),vInitOmegaChi,options);    
                % Second run to refine
                options = optimoptions('fsolve','Display','iter','Algorithm','trust-region',...
                'MaxFunctionEvaluations',100000,'StepTolerance',1e-12,'OptimalityTolerance',1e-12);
                [vSolOmegaChi,~,exitflag] = fsolve(@(x) fn_omega_chi_gamma(x,sPar,target),vSolOmegaChi0,options);
                sPar.omega = vSolOmegaChi(1);
                sPar.chi = vSolOmegaChi(2);
            
            case 'fixed'
                % values for omega and chi here come from 'death' model with sigma_a = 0.003
                sPar.omega = 0.635860132254258;
                sPar.chi = 0.756773618926686;
        end
end

%% Steady-state 
%---------------------------------------------------------------------------

% Solve jointly for aHat and h (cutoff function accounts for gamma)
vInit = [0.001,0.7];  
options = optimoptions('fsolve','Display','none',...
        'MaxFunctionEvaluations',100000,'StepTolerance',1e-12,'OptimalityTolerance',1e-12);  
vSol = fsolve(@(x) fn_Cutoff_gamma(x,sPar),vInit,options);
vSol = real(vSol);
aHat = vSol(1);   h = vSol(2);

% Compute p and a* implied by aHat
p = 1 - (aHat-sPar.aL)/(sPar.aH-sPar.aL);  % 1-F(aHat)
aStar = 0.5*(sPar.aH+aHat);                % E(a|a>aHat)

% Compute other relevant variables
v     = p*sPar.Upsilon/(1+p*(((1-sPar.delta)+(sPar.delta*sPar.gamma)/(1-sPar.gamma))*(h/sPar.delta))); % account for gamma
n     = h*v/sPar.delta;
f     = sPar.psi*((1-(1-sPar.delta)*n)/v)^(sPar.alpha-1);
y     = (aStar+sPar.ss.z)*n;
LambdaaStar =  ((1-sPar.omega)*(sPar.xss*(aStar+sPar.ss.z)-sPar.chi));
LambdaaHat =  ((1-sPar.omega)*(sPar.xss*(aHat+sPar.ss.z)-sPar.chi));

% Implied value functions need to be solved for jointly in this variant
eqn = @(J) [-J(1) + sPar.beta*J(1) - p*sPar.kappa+p*h*(J(3)-sPar.beta*J(1));  %  J^U
            -J(2) + (1-sPar.gamma)*J(1) + sPar.gamma*sPar.beta*J(2);          %  J^U^tilde (!= beta*JU)
            -J(3) + (1-sPar.omega)*(sPar.xss*(aStar+sPar.ss.z)-sPar.chi)+sPar.beta*((1-sPar.delta)*J(3)+sPar.delta*J(2))]; % J(a^*)
    
J0 = [0.1, 0.1, 0.3];
options = optimoptions('fsolve','Display','iter');
[JSol,fval] = fsolve(eqn,J0,options);
JU = JSol(1);
JUtilde = JSol(2);
JaStar = JSol(3);
JaHat = ((1-sPar.omega)*(sPar.xss*(aHat+sPar.ss.z)-sPar.chi)+sPar.beta*sPar.delta*JUtilde)/(1-sPar.beta*(1-sPar.delta));

% Quick check -- this should be 0 given recalibration
JDiff = JaHat-sPar.beta*JU;     
hImplied = sPar.kappa./JDiff; 
disp('Check hiring rate -- difference between solver and implied should be 0')
disp(hImplied-h)

% Put everything together into a sub-element of the parameter structure 
% (called by Dynare)
sPar.ss.n = n;
sPar.ss.m = sPar.delta*sPar.ss.n;
sPar.ss.us = 1-(1-sPar.delta)*sPar.ss.n;
sPar.ss.v = v;
sPar.ss.h = h;
sPar.ss.aHat = aHat;
sPar.ss.aStar = aStar;
sPar.ss.p = p;
sPar.ss.JaHat = JaHat;
sPar.ss.JStar = JaStar;
sPar.ss.JU = JU;

save inputParameters sPar

%% Run Dynare
%---------------------------------------------------------------------------
dynare dynareSAMOptionValue noclearall

%% Get IRFs, make relevant manipulations for visualization, store results
%---------------------------------------------------------------------------

% Select relevant IRFs
switch sPlotting.options.IRFs  
    case 'Pure uncertainty'
          mIRF_1 = mIRFProp_zUncertainty_EMAS;
          sPlotting.numModels = 1;
    case 'Pure uncertainty and total volatility'
          mIRF_1 = mIRFProp_zUncertainty_EMAS;
          mIRF_2 = mIRFProp_Andreasen_zVol_EMAS;
          sPlotting.numModels = 2;
          sPlotting.legend.labels = {'Pure uncertainty effect';'Total volatility effect'};
end

sPlotting.options.legend.show = min(sPlotting.numModels-1,1); % only show if >1 models

% Put together
if sPlotting.numModels == 1
    sResults.aIRF = reshape([mIRF_1],[size(mIRF_1,1),size(mIRF_1,2),sPlotting.numModels]);
elseif sPlotting.numModels == 2
     sResults.aIRF = reshape([mIRF_1,mIRF_2],[size(mIRF_1,1),size(mIRF_1,2),sPlotting.numModels]);
end

% For rates (u, h, f, p), switch to ppt instead of percent
uPos = strmatch('u',vNames,'exact');
hPos = strmatch('h',vNames,'exact');
fPos = strmatch('f',vNames,'exact');
pPos = strmatch('p',vNames,'exact');

sResults.aIRF(:,uPos,:) = sResults.aIRF(:,uPos,:)*vEMAS(uPos);
sResults.aIRF(:,hPos,:) = sResults.aIRF(:,hPos,:)*vEMAS(hPos);
sResults.aIRF(:,fPos,:) = sResults.aIRF(:,fPos,:)*vEMAS(fPos);
sResults.aIRF(:,pPos,:) = sResults.aIRF(:,pPos,:)*vEMAS(pPos);

% Look at absolute deviations from relevant value functions
aHatPos = strmatch('aHat',vNames,'exact');   % 0 for p = 1/2
JU_dPos = strmatch('JU_d',vNames,'exact');
JaHatPos = strmatch('JaHat',vNames,'exact');
JaHat_aHatSSPos = strmatch('JaHat_aHatSS',vNames,'exact');
sResults.aIRF(:,JU_dPos,:) = sResults.aIRF(:,JU_dPos,:)*vEMAS(JU_dPos);
sResults.aIRF(:,JaHatPos,:) = sResults.aIRF(:,JaHatPos,:)*vEMAS(JaHatPos);
sResults.aIRF(:,JaHat_aHatSSPos,:) = sResults.aIRF(:,JaHat_aHatSSPos,:)*vEMAS(JaHat_aHatSSPos);
sResults.aIRF(:,aHatPos,:) = sResults.aIRF(:,aHatPos,:)*vEMAS(aHatPos);

% Adjust values that are virtually zero (not exactly due to use of a solver) to exactly zero to avoid confusing graphs
sResults.aIRF(abs(sResults.aIRF)<1e-10)=0;

sResults.names = vNames;

%% Plotting 
%---------------------------------------------------------------------------

fig1=figure;
for iV = 1:numel(sPlotting.variables)
    subplot(sPlotting.numSubplotV,sPlotting.numSubplotH,iV)
    for iN = 1:sPlotting.numModels
        p1=plot(0:sPlotting.IrfPeriodsNum,100*sResults.aIRF(1:sPlotting.IrfPeriodsNum+1,strmatch(sPlotting.variables{iV},vNames,'exact'),iN),sSettings.lines.list{iN},'LineWidth',sSettings.lines.width,'Color',sSettings.colors.list{iN});
        hold on
        plot(0:sPlotting.IrfPeriodsNum,zeros(sPlotting.IrfPeriodsNum+1,1),sSettings.lines.list{3},'HandleVisibility','off','Color','k','Linewidth',0.5);
        xlim([0 sPlotting.IrfPeriodsNum]);
        set(gca,'XTick',[0:2:sPlotting.IrfPeriodsNum],'FontSize',sSettings.font.size.axisticks,'fontname',sSettings.font.name);
        ax = gca;
        ax.YAxis.Exponent = 0; 
        box on
    end
    title(sPlotting.names{iV},'FontSize',sSettings.font.size.default,'fontname',sSettings.font.name,'FontWeight','normal');

   if ismember(sPlotting.variables(iV),{'u','f','h','p'}) == 1
        ylabel('Deviation (ppts.)','FontSize',sSettings.font.size.axis,'fontname',sSettings.font.name);
    elseif ismember(sPlotting.variables(iV),{'aHat','JU_d','JaHat_aHatSS','JU1'}) == 1
        ylabel('Deviation (abs.)','FontSize',sSettings.font.size.axis,'fontname',sSettings.font.name);
    else        
       ylabel('Deviation (pct.)','FontSize',sSettings.font.size.axis,'fontname',sSettings.font.name);
    end    

if iV>=numel(sPlotting.variables)-1
    xlabel('Time (quarters)','FontSize',sSettings.font.size.axis,'fontname',sSettings.font.name);
end

end

if sPlotting.options.legend.show == 1
    subplot(sPlotting.numSubplotV,sPlotting.numSubplotH,1)
    legend1 = legend(sPlotting.legend.labels);
    set(legend1,'fontname','times','Location','best','FontSize',sSettings.font.size.legend);
end
   
%% Printing
%---------------------------------------------------------------------------

sSettings.plots.xSize = sPlotting.numSubplotH*8.75; 
sSettings.plots.ySize = sPlotting.numSubplotV*6.25; 
sSettings.plots.xCut = 1; sSettings.plots.yCut = 1;

set(gcf,'Units','centimeters','Position',[0 0 sSettings.plots.xSize sSettings.plots.ySize],'PaperUnits','centimeters' ...
      ,'PaperPosition',[0 0 sSettings.plots.xSize sSettings.plots.ySize],'PaperSize',[sSettings.plots.xSize-sSettings.plots.xCut sSettings.plots.ySize-sSettings.plots.yCut],'PaperPositionMode','auto')
  
if OptionPrint == 1      
    FigNamepdf =horzcat(horzcat(TargetPath,FigName),'.pdf');
    print(fig1,FigNamepdf,'-dpdf','-painters')
end

%% Done!
%---------------------------------------------------------------------------
