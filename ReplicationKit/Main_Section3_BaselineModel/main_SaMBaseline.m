%%=========================================================================
% Replication code for den Haan, Freund and Rendahl (2021)
%
% Baseline model
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
% Choose how many different models (or IRFs) to show
sPlotting.options.IRFs = 'Pure uncertainty and total volatility'; % choose b/w 'Pure uncertainty', 'Pure uncertainty and total volatility', 'Productivity',

% Choose whether to print and, if so, set figure name and target folder to save figure in
sPlotting.options.print = 0;        
sPlotting.options.targetPath ='.\Output\Figures\';
sPlotting.options.figName = 'fig_Baseline_WNash';
sPlotting.options.color = 0;

%% Model specification
%---------------------------------------------------------------------------
% NB: choice directly in mod file: wage specification (Nash vs linear), 
% stochastic processes, treatment of vacancy posting costs

% Let's set the parameter values here for sake of transparency
sPar.zBar  = 1;  
sPar.hBar  = 0.7;                         
sPar.uBar  = 0.064; 
sPar.beta  = 0.99;   
sPar.eta   = 10;    
sPar.alpha = 0.5; % change to 0.2 (or some similar value) for Fig 3 
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

% For the two wage specifications, we proceed as follows.
% For Nash, we just follow Leduc-Liu and set the bargaining weight to 1/2;
% from the wage eqn we then get the outside option
sPar.omegaN     = 1/2;
sPar.chiN       = (1/(1-sPar.omegaN))*(sPar.ss.w-sPar.omegaN*(sPar.xss*sPar.zBar+sPar.beta*(1-sPar.delta)*sPar.kappa*sPar.ss.theta));
FSN             = sPar.xss*sPar.zBar-sPar.chiN; % fundamental surplus
TZElasticityN   =(sPar.xss*sPar.zBar/FSN)*((1-sPar.beta*(1-sPar.delta)+sPar.omegaN*sPar.beta*(1-sPar.delta)*sPar.ss.theta*sPar.hBar)/(sPar.alpha*(1-sPar.beta*(1-sPar.delta))+sPar.omegaN*sPar.beta*(1-sPar.delta)*sPar.ss.theta*sPar.hBar));
 
% For the linear wage, we just target the same steady-state elasticity of $\theta$ w.r.t. z and compute the
% implied outside option and bargaining weight
TZElasticityL   = TZElasticityN;
FSL             = (sPar.xss*sPar.zBar)/(sPar.alpha*TZElasticityL);
sPar.chiL       = sPar.xss*sPar.zBar-FSL;
sPar.omegaL     = (sPar.ss.w-sPar.chiL)/(sPar.xss*sPar.zBar-sPar.chiL);

% Save parameters and clean up
save inputParameters sPar
clear FSL FSN TZElasticityL TZElasticityN

%% Generic design
%---------------------------------------------------------------------------
% Set these explicitly ones; same settings used in the rest of the paper
% and just loaded via inputSettings.map

% Design
sSettings.plots.color = sPlotting.options.color;
if sSettings.plots.color == 1
    sSettings.colors.list = {[0,42,91]/255,[255 127 0]/255,[128,128,128]/255,[114,47,55]/255};
else
    sSettings.colors.list = {[0.2,0.2,0.2],[0.33 0.33 0.33],[0.46,0.46,0.46],[0.6,0.6,0.6]};
end
sSettings.lines.list        = {'-';'--';':';'-.'};
sSettings.lines.width       = 1.6;
sSettings.font.name         = 'times';
sSettings.lines.width       = 1.6;
sSettings.lines.zeros       = '-';
sSettings.font.size.default = 10;
sSettings.font.size.legend  = 8;
sSettings.font.size.axis    = 10;
sSettings.font.size.axisticks = 8;
sSettings.font.name         = 'times';
sSettings.plots.xSize       = 2*8.75; 
sSettings.plots.ySize       = 3*6.25; 
sSettings.plots.xCut        = 2;
sSettings.plots.yCut        = 0.5;

save inputSettings sSettings
save inputPlotting sPlotting

%% Run dynare 
%--------------------------------------------------------------------------

% Call dynare
dynare dynareSaMBaseline 

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
          sPlotting.options.legend.labels = {'Pure uncertainty effect';'Total volatility effect'};
    case 'Productivity'
          mIRF_1 = mIRFProp_GIRF_Order1_z_EMAS;
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

%% Manual plotting 
%--------------------------------------------------------------------------
% Here manually chosen
sPlotting.IRFPeriods = 10; 
sPlotting.numSubplotV = 3;
sPlotting.numSubplotH = 2;

figure;
for iN = 1:sPlotting.numModels
    
subplot(sPlotting.numSubplotV,sPlotting.numSubplotH,1)
p1=plot(0:sPlotting.IRFPeriods,100*sResults.aIRF(1:sPlotting.IRFPeriods+1,strmatch('J',vNames,'exact'),iN),sSettings.lines.list{iN},'LineWidth',sSettings.lines.width,'Color',sSettings.colors.list{iN});
hold on
plot(0:sPlotting.IRFPeriods,zeros(sPlotting.IRFPeriods+1,1),sSettings.lines.zeros,'HandleVisibility','off','Color','k','Linewidth',0.1);
xlim([0 sPlotting.IRFPeriods]);
ylim([-0.125 0.125]);
set(gca,'XTick',[0:2:sPlotting.IRFPeriods],'FontSize',sSettings.font.size.axisticks,'fontname',sSettings.font.name);
ax = gca; ax.YAxis.Exponent = 0; 
box on
title('Match value','FontSize',sSettings.font.size.default,'fontname',sSettings.font.name,'FontWeight','normal');
ylabel('Deviation (pct.)','FontSize',sSettings.font.size.axis,'fontname',sSettings.font.name);

subplot(sPlotting.numSubplotV,sPlotting.numSubplotH,2)
p1=plot(0:sPlotting.IRFPeriods,100*sResults.aIRF(1:sPlotting.IRFPeriods+1,strmatch('w',vNames,'exact'),iN),sSettings.lines.list{iN},'LineWidth',sSettings.lines.width,'Color',sSettings.colors.list{iN});
hold on
plot(0:sPlotting.IRFPeriods,zeros(sPlotting.IRFPeriods+1,1),sSettings.lines.zeros,'HandleVisibility','off','Color','k','Linewidth',0.1);
xlim([0 sPlotting.IRFPeriods]);
ylim([-0.00575 0.00575]);
set(gca,'XTick',[0:2:sPlotting.IRFPeriods],'FontSize',sSettings.font.size.axisticks,'fontname',sSettings.font.name);
ax = gca; ax.YAxis.Exponent = 0; 
box on
title('Wage','FontSize',sSettings.font.size.default,'fontname',sSettings.font.name,'FontWeight','normal');
ylabel('Deviation (pct.)','FontSize',sSettings.font.size.axis,'fontname',sSettings.font.name);

subplot(sPlotting.numSubplotV,sPlotting.numSubplotH,3)
p1=plot(0:sPlotting.IRFPeriods,100*sResults.aIRF(1:sPlotting.IRFPeriods+1,strmatch('theta',vNames,'exact'),iN),sSettings.lines.list{iN},'LineWidth',sSettings.lines.width,'Color',sSettings.colors.list{iN});
hold on
plot(0:sPlotting.IRFPeriods,zeros(sPlotting.IRFPeriods+1,1),sSettings.lines.zeros,'HandleVisibility','off','Color','k','Linewidth',0.1);
xlim([0 sPlotting.IRFPeriods]);
ylim([-0.25 0.25]);
set(gca,'XTick',[0:2:sPlotting.IRFPeriods],'FontSize',sSettings.font.size.axisticks,'fontname',sSettings.font.name);
ax = gca; ax.YAxis.Exponent = 0; 
box on
title('Labor market tightness','FontSize',sSettings.font.size.default,'fontname',sSettings.font.name,'FontWeight','normal');
ylabel('Deviation (pct.)','FontSize',sSettings.font.size.axis,'fontname',sSettings.font.name);

subplot(sPlotting.numSubplotV,sPlotting.numSubplotH,4)
p1=plot(0:sPlotting.IRFPeriods,100*sResults.aIRF(1:sPlotting.IRFPeriods+1,strmatch('h',vNames,'exact'),iN),sSettings.lines.list{iN},'LineWidth',sSettings.lines.width,'Color',sSettings.colors.list{iN});
hold on
plot(0:sPlotting.IRFPeriods,zeros(sPlotting.IRFPeriods+1,1),sSettings.lines.zeros,'HandleVisibility','off','Color','k','Linewidth',0.1);
xlim([0 sPlotting.IRFPeriods]);
ylim([-0.25 0.25]);  % figure 1
%ylim([-0.125 0.125]); % figure 2
%ylim([-0.025 0.025]); % figure 3
set(gca,'XTick',[0:2:sPlotting.IRFPeriods],'FontSize',sSettings.font.size.axisticks,'fontname',sSettings.font.name);
ax = gca; ax.YAxis.Exponent = 0; 
box on
title('Hiring rate','FontSize',sSettings.font.size.default,'fontname',sSettings.font.name,'FontWeight','normal');
ylabel('Deviation (ppts.)','FontSize',sSettings.font.size.axis,'fontname',sSettings.font.name);

subplot(sPlotting.numSubplotV,sPlotting.numSubplotH,5)
p1=plot(0:sPlotting.IRFPeriods,100*sResults.aIRF(1:sPlotting.IRFPeriods+1,strmatch('f',vNames,'exact'),iN),sSettings.lines.list{iN},'LineWidth',sSettings.lines.width,'Color',sSettings.colors.list{iN});
hold on
plot(0:sPlotting.IRFPeriods,zeros(sPlotting.IRFPeriods+1,1),sSettings.lines.zeros,'HandleVisibility','off','Color','k','Linewidth',0.1);
xlim([0 sPlotting.IRFPeriods]);
ylim([-0.125 0.125]);
set(gca,'XTick',[0:2:sPlotting.IRFPeriods],'FontSize',sSettings.font.size.axisticks,'fontname',sSettings.font.name);
ax = gca; ax.YAxis.Exponent = 0; 
box on
title('Job finding rate','FontSize',sSettings.font.size.default,'fontname',sSettings.font.name,'FontWeight','normal');
ylabel('Deviation (ppts.)','FontSize',sSettings.font.size.axis,'fontname',sSettings.font.name);
xlabel('Time (quarters)','FontSize',sSettings.font.size.axis,'fontname',sSettings.font.name);

subplot(sPlotting.numSubplotV,sPlotting.numSubplotH,6)
p1=plot(0:sPlotting.IRFPeriods,100*sResults.aIRF(1:sPlotting.IRFPeriods+1,strmatch('u',vNames,'exact'),iN),sSettings.lines.list{iN},'LineWidth',sSettings.lines.width,'Color',sSettings.colors.list{iN});
hold on
plot(0:sPlotting.IRFPeriods,zeros(sPlotting.IRFPeriods+1,1),sSettings.lines.zeros,'HandleVisibility','off','Color','k','Linewidth',0.1);
xlim([0 sPlotting.IRFPeriods]);
ylim([-0.055 0.055]); % figure 1
%ylim([-0.025 0.025]); % figures 2 and 3
set(gca,'XTick',[0:2:sPlotting.IRFPeriods],'FontSize',sSettings.font.size.axisticks,'fontname',sSettings.font.name);
ax = gca; ax.YAxis.Exponent = 0; 
box on
title('Unemployment rate','FontSize',sSettings.font.size.default,'fontname',sSettings.font.name,'FontWeight','normal');
ylabel('Deviation (ppts.)','FontSize',sSettings.font.size.axis,'fontname',sSettings.font.name);
xlabel('Time (quarters)','FontSize',sSettings.font.size.axis,'fontname',sSettings.font.name);

end

if sPlotting.options.legend.show == 1
subplot(sPlotting.numSubplotV,sPlotting.numSubplotH,1)
legend1 = legend(sPlotting.options.legend.labels);
set(legend1,'fontname',sSettings.font.name,'Location','north','FontSize',sSettings.font.size.legend,'interpreter','none')
end

clear hfig ax fPos hPos iN legend1 mIRF_1 mIRF_2  p1 uPos vEMAS  vNames_GIRF mIRFProp_zUncertainty_EMAS mIRFProp_GIRF_zVol_EMAS mIRFProp_GIRF_Order1_z_EMAS
%  vNames
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
clear figNamepdf
!del dynareSaMBaseline.jnl dynareSaMBaseline.log dynareSaMBaseline.m dynareSaMBaseline_dynamic.m 
!del dynareSaMBaseline_set_auxiliary_variables.m dynareSaMBaseline_static.m main_SaM_Baseline.asv
!del dynareSaMBaseline_results.mat

%% 'Appendix'
%----------------------------------------------------------------------------

%% Remark on the (G)IRF for productivity level
% Why did I use "z_Obs", "sigma_z_Obs" and "n_Obs" above? Purely out of laziness and to "trick"
% Dynare into believing that these are control rather than state variables. 
% Else you need to read in a separate Andreasen-method-generated IRFs for state variables.
% Can compare this to:
% Test1 = mIRFProp_Andreasen_zVol_state(:,2); % positions; outDynare.label_v
% Test2 = mIRFProp_Andreasen_zVol_EMAS(:,strmatch('z_Obs',vNames))

%% Alternative plotting method
% If users want to explore variables other than the ones preselected above...

%{
sPlotting.vVariables = {'J','w','theta','h','f','u'};
sPlotting.vVNames = {'Firm value', 'Wage','Labor market tightness','Hiring rate','Job finding rate','Unemployment rate'};
sPlotting.numSubplotV = 3;
sPlotting.numSubplotH = 2;
sPlotting.IRFPeriods = 10;

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
        ylabel('Deviation (ppts.)','FontSize',sSettings.font.size.axis,'fontname',sSettings.font.name);
    else
       ylabel('Deviation (pct.)','FontSize',sSettings.font.size.axis,'fontname',sSettings.font.name);
    end    
    if iV >= (sPlotting.numSubplotV*sPlotting.numSubplotH)-1
       xlabel('Time (quarters)','FontSize',sSettings.font.size.axis,'fontname',sSettings.font.name);
    end
end
%}