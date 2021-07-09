%%=========================================================================
% Replication code for den Haan, Freund and Rendahl (2021)
%
% This file: run version with risk aversion and plot results
% for a single model
%
% Run on Matlab R2019b, Dynare 4.4.3
% Last updated: September 2020
% For any questions, please email lukas.beat.freund@gmail.com
%%=========================================================================

%% Housekeeping
%--------------------------------------------------------------------------
clear; close all; clc;

%% Main plotting choices
%--------------------------------------------------------------------------
% Load standard design settings
load inputSettings

% Choose how many different models (or IRFs) to show
sPlotting.options.IRFs = 'Pure uncertainty'; % choose b/w 'Pure uncertainty', 'Pure uncertainty and total volatility', 'Productivity'

% Choose whether to print and, if so, set figure name and target folder to save figure in
sPlotting.options.print = 0;        
sPlotting.options.targetPath ='.\Output\Figures\';
sPlotting.options.figName = 'fig_SaM_Baseline_Specification';
sPlotting.options.color = 0;

% Variables to show     
sPlotting.vVariables = {'J','w','u','c'};
sPlotting.vVNames = {'Match value', 'Wage','Unemployment rate','Consumption'};
sPlotting.numSubplotV = 2;
sPlotting.numSubplotH = 2;
sPlotting.IRFPeriods = 10;
save inputPlotting sPlotting

%% Model specification
%---------------------------------------------------------------------------
% NB: choice directly in mod file: wage specification (Nash vs linear), 
% stochastic processes, treatment of vacancy posting costs

% Let's set the parameter values here for sake of transparency
sPar.xi    = 1; % = 0 implies we're back to the baseline case
sPar.zBar  = 1;  
sPar.hBar  = 0.7;                         
sPar.uBar  = 0.064; 
sPar.beta  = 0.99;   
sPar.eta   = 10;    
sPar.alpha = 0.5; 
sPar.delta = 0.1;                      
sPar.kappa = 0.14;     

% Productivity process
sPar.rho_z         = 0.95;           
sPar.sigma_zbar    = 0.01;    
sPar.rho_sigma_z   = 0.76;       
sPar.sigma_sigma_z = 0.392; % if use stochastic process in levels: 0.392*sPar.sigma_zbar; 

% This implies...
sPar.ss.z       = sPar.zBar;
sPar.ss.sigma_z = sPar.sigma_zbar;
sPar.ss.u       = sPar.uBar;
sPar.ss.y       = sPar.ss.z*(1-sPar.ss.u);
sPar.ss.m       = sPar.delta*(1-sPar.ss.u);
sPar.ss.us      = 1-(1-sPar.delta)*((1-sPar.uBar));
sPar.ss.f       = sPar.ss.m/sPar.ss.us;
sPar.ss.v       = sPar.ss.m/sPar.hBar;
sPar.ss.theta   = sPar.ss.v/sPar.ss.us;
sPar.psi        = sPar.ss.m/(sPar.ss.us^sPar.alpha*sPar.ss.v^(1-sPar.alpha));
sPar.xss        = (sPar.eta-1)/sPar.eta;
sPar.ss.J       = (sPar.kappa/sPar.hBar);
sPar.ss.w       = sPar.xss*sPar.zBar-(1-sPar.beta*(1-sPar.delta))*sPar.ss.J; 
sPar.ss.c       = 1-sPar.uBar;
sPar.ss.lambda  = sPar.ss.c^(-sPar.xi); 

% Wage specification: Nash
% NB: this differs from the baseline model for two reasons
% (i) integrates flow unemployment benefits phi in Nash case
sPar.phi = 0.25;
% (ii) allows for risk aversion, which alters the calibration 

% We follow Leduc-Liu and set the bargaining weight to 1/2;
% from the wage eqn we then get the outside option
sPar.omegaN = 1/2;
sPar.chiN = sPar.ss.lambda*((1/(1-sPar.omegaN))*sPar.ss.w...
    - (sPar.omegaN/(1-sPar.omegaN))*(sPar.xss*sPar.zBar+sPar.beta*(1-sPar.delta)*sPar.kappa*sPar.ss.theta)-sPar.phi);
FSN = sPar.xss*sPar.zBar-(sPar.phi+sPar.chiN/sPar.ss.lambda);  % notice that (vPar.phi+vPar.chiN/lambdass) corresponds
% to vPar.chiN in the risk-neutral baseline model that abstracts from phi
TZElasticityN =(sPar.xss*sPar.zBar/FSN)*((1-sPar.beta*(1-sPar.delta)...
                +sPar.omegaN*sPar.beta*(1-sPar.delta)*sPar.ss.theta*sPar.hBar)/(sPar.alpha*(1-sPar.beta*(1-sPar.delta))...
                +sPar.omegaN*sPar.beta*(1-sPar.delta)*sPar.ss.theta*sPar.hBar));
 
% Linear wage: target the same elasticity of labor
% market tightness w.r.t. productivity as implied by Nash bargaining under the above choices and compute the
% implied outside option and bargaining weight
% For the linear wage, we just target the same steady-state elasticity of $\theta$ w.r.t. z and compute the
% implied outside option and bargaining weight
TZElasticityL   = TZElasticityN;
FSL             = (sPar.xss*sPar.zBar)/(sPar.alpha*TZElasticityL);
sPar.chiL       = sPar.xss*sPar.zBar-FSL;
sPar.omegaL     = (sPar.ss.w-sPar.chiL)/(sPar.xss*sPar.zBar-sPar.chiL);

% Save parameters and clean up
save inputParameters sPar
clear FSL FSN TZElasticityL TZElasticityN

%% Run dynare 
%--------------------------------------------------------------------------

% Call dynare
dynare dynareSaM_riskAversion 

% To keep things clean, I find it convenient to clear the workspace and reload
% only what's needed
clear 
load inputSettings;
load inputPlotting;
load inputParameters;
load('Output/IRFs/IRFs_Spc.mat'); % loads the most recently created IRFs

% Select relevant IRFs
switch sPlotting.options.IRFs  
    case 'Pure uncertainty'
          mIRF_1 = mIRFProp_zUncertainty_EMAS;
          sPlotting.numModels = 1;
    case 'Pure uncertainty and total volatility'
          mIRF_1 = mIRFProp_zUncertainty_EMAS;
          mIRF_2 = mIRFProp_GIRF_zVol_EMAS;
          sPlotting.numModels = 2;
          sPlotting.legend.labels = {'Pure uncertainty effect';'Total volatility effect'};
    case 'Productivity'
          mIRF_1 = mIRFProp_GIRF_z_EMAS; % first order
          sPlotting.numModels = 1;
end

sPlotting.options.legend.show = min(sPlotting.numModels-1,1); % only show if >1 models

% Put together
if sPlotting.numModels == 1
sResults.aIRF = reshape([mIRF_1],[size(mIRF_1,1),size(mIRF_1,2),sPlotting.numModels]);
elseif sPlotting.numModels == 2
sResults.aIRF = reshape([mIRF_1,mIRF_2],[size(mIRF_1,1),size(mIRF_1,2),sPlotting.numModels]);
end
sResults.aIRF (abs(sResults.aIRF )<1e-10)=0; % adjust values that are virtually zero (not exactly due to use of a solver) to exactly zero to avoid confusing graphs

% For rates (u, h, f), switch to ppt instead of percent
uPos = strmatch('u',vNames,'exact');
fPos = strmatch('f',vNames,'exact');
hPos = strmatch('h',vNames,'exact');
sResults.aIRF(:,uPos,:) = sResults.aIRF(:,uPos,:)*vEMAS(uPos);
sResults.aIRF(:,fPos,:) = sResults.aIRF(:,fPos,:)*vEMAS(fPos);
sResults.aIRF(:,hPos,:) = sResults.aIRF(:,hPos,:)*vEMAS(hPos);

%% Plotting 
%--------------------------------------------------------------------------

figure;
for iV = 1:numel(sPlotting.vVariables)
    subplot(sPlotting.numSubplotV,sPlotting.numSubplotH,iV)
for iN = 1:sPlotting.numModels
    box on 
    hold on
    plot(0:sPlotting.IRFPeriods,100*sResults.aIRF(1:sPlotting.IRFPeriods+1,strmatch(sPlotting.vVariables{iV},vNames,'exact'),iN),sSettings.lines.list{iN},'LineWidth',sSettings.lines.width,'Color',sSettings.colors.list{iN});
    hold on
    plot(0:sPlotting.IRFPeriods,zeros(sPlotting.IRFPeriods+1,1),sSettings.lines.zeros,'HandleVisibility','off','Color','k','Linewidth',0.1);
    xlim([0 sPlotting.IRFPeriods]);
    set(gca,'XTick',[0:2:sPlotting.IRFPeriods],'FontSize',sSettings.font.size.axisticks,'fontname',sSettings.font.name);
    ax = gca;     ax.YAxis.Exponent = 0; 
end
    title(sPlotting.vVNames{iV},'FontSize',sSettings.font.size.default,'fontname',sSettings.font.name,'FontWeight','normal');
    
    if ismember(sPlotting.vVariables(iV),{'u','f','h'}) == 1
        ylabel('Percentage poins','FontSize',sSettings.font.size.axis,'fontname',sSettings.font.name);
    else
       ylabel('Percent','FontSize',sSettings.font.size.axis,'fontname',sSettings.font.name);
    end    
    if iV >= (sPlotting.numSubplotV*sPlotting.numSubplotH)-1
       xlabel('Time (quarters)','FontSize',sSettings.font.size.axis,'fontname',sSettings.font.name);
    end
end

if sPlotting.options.legend.show == 1
subplot(sPlotting.numSubplotV,sPlotting.numSubplotH,1)
legend1 = legend(sPlotting.legend.labels);
set(legend1,'fontname',sSettings.font.name,'Location','best','FontSize',sSettings.font.size.legend,'interpreter','none')
end

clear hfig ax fPos hPos iN legend1 mIRF_1 mIRF_2  p1 uPos vEMAS  vNames_GIRF mIRFProp_zUncertainty_EMAS mIRFProp_GIRF_zVol_EMAS mIRFProp_GIRF_Order1_z_EMAS

%% Print 
%--------------------------------------------------------------------------

set(gcf,'Units','centimeters','Position',[0 0 sSettings.plots.xSize sSettings.plots.ySize],'PaperUnits','centimeters' ...
      ,'PaperPosition',[0 0 sSettings.plots.xSize sSettings.plots.ySize],'PaperSize',[sSettings.plots.xSize-sSettings.plots.xCut sSettings.plots.ySize-sSettings.plots.yCut],'PaperPositionMode','auto')

if sPlotting.options.print == 1      
    if sPlotting.options.color == 1
         figNamepdf =horzcat(horzcat(sPlotting.options.targetPath,horzcat(sPlotting.options.figName,'_CS')),'.pdf');
    else
         figNamepdf =horzcat(horzcat(sPlotting.options.targetPath,sPlotting.options.figName),'.pdf');
    end
  print(gcf,figNamepdf,'-dpdf','-painters')
end

%% Done!
%----------------------------------------------------------------------------
% delete a few files to keep the folder clean
clear figNamepdf iV vNames
!del dynareSaM_riskAversion.jnl dynareSaM_riskAversion.log dynareSaM_riskAversion.m dynareSaM_riskAversion_dynamic.m 
!del dynareSaM_riskAversion_set_auxiliary_variables.m dynareSaM_riskAversion_static.m main_SaM_riskAversion.asv
!del dynareSaM_riskAversion_results.mat

