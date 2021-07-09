%%=========================================================================
% Replication code for den Haan, Freund and Rendahl (2020)
% Extension to risk-aversion (used in Appendix B)
% Run on Matlab R2019b, Dynare 4.4.3
% Last updated: August 2020
% For any questions please email me at lukas.beat.freund@gmail.com
%%=========================================================================
 
%% Housekeeping
%--------------------------------------------------------------------------
clear;
close all;
clc;
TimeStart = tic;

%% User choices
%---------------------------------------------------------------------------
% Choose models (1: pure uncertainty effect only; 2: also total vol effect)
NumModels = 2;

% Choose whether to print and, if so, set figure name and target folder to save figure in
OptionPrint = 0;        
TargetPath ='.\Output\';
FigName = 'fig_Specification';

%% Variables
vVariables = {'J','w','u','c'};
vVNames = {'Firm value', 'Wage','Unemployment rate','Consumption'};

NumSubplotV = 2;
NumSubplotH = 2;
OptionLegend = NumModels-1; % only use if 2 models
Labels{1}  = 'Pure uncertainty effect';
Labels{2} =  'Total volatility effect';

%% Parameters

% Choose between wage specification (Nash vs linear), stochastic processes,
% treatment of vacancy posting costs: directly in .mod file

% Let's set the parameter values here for sake of transparency
vPar.xi    = 1; % = 0 implies we're back to the baseline case
vPar.zBar  = 1;  
vPar.hBar  = 0.7;                         
vPar.uBar  = 0.064; 
vPar.beta  = 0.99;   
vPar.eta   = 10;    
vPar.alpha = 0.5; 
vPar.delta = 0.1;                      
vPar.kappa = 0.14;     

% This implies...
mss = vPar.delta*(1-vPar.uBar);
usss = 1-(1-vPar.delta)*((1-vPar.uBar));
fss = mss/usss;
vss = mss/vPar.hBar;
thetass= vss/usss;
xss = (vPar.eta-1)/vPar.eta;
Jss = (1/vPar.hBar)*(vPar.kappa);
wss = xss*vPar.zBar-(1-vPar.beta*(1-vPar.delta))*Jss; 
css = 1-vPar.uBar;
lambdass = css^(-vPar.xi); 

% Wage specification: Nash
% NB: this differs from the baseline model for two reasons
% (i) integrates flow unemployment benefits phi in Nash case
vPar.phi = 0.25;
% (ii) allows for risk aversion, which alters the calibration 

% We follow Leduc-Liu and set the bargaining weight to 1/2;
% from the wage eqn we then get the outside option
vPar.omegaN = 1/2;
vPar.chiN = lambdass*((1/(1-vPar.omegaN))*wss...
    - (vPar.omegaN/(1-vPar.omegaN))*(xss*vPar.zBar+vPar.beta*(1-vPar.delta)*vPar.kappa*thetass)-vPar.phi);
FSN = xss*vPar.zBar-(vPar.phi+vPar.chiN/lambdass);  % notice that (vPar.phi+vPar.chiN/lambdass) corresponds
% to vPar.chiN in the risk-neutral baseline model that abstracts from phi
TZElasticityN =(xss*vPar.zBar/FSN)*((1-vPar.beta*(1-vPar.delta)...
                +vPar.omegaN*vPar.beta*(1-vPar.delta)*thetass*vPar.hBar)/(vPar.alpha*(1-vPar.beta*(1-vPar.delta))...
                +vPar.omegaN*vPar.beta*(1-vPar.delta)*thetass*vPar.hBar));
 
% Linear wage: target the same elasticity of labor
% market tightness w.r.t. productivity as implied by Nash bargaining under the above choices and compute the
% implied outside option and bargaining weight
TZElasticityL = TZElasticityN;
FSL = (xss*vPar.zBar)/(vPar.alpha*TZElasticityL);
vPar.chiL = xss*vPar.zBar-FSL;                              % effectively includes phi
vPar.omegaL = (wss-vPar.chiL)/(xss*vPar.zBar-vPar.chiL);

% Productivity process
vPar.rho_z = 0.95;           
vPar.sigma_zbar = 0.01;    
vPar.rho_sigma_z = 0.76;       
vPar.sigma_sigma_z = 0.392; % levels: 0.392*vPar.sigma_zbar; 

% Save parameters 
save Parameters vPar

%% Design choices
%---------------------------------------------------------------------------
IRFPeriodsNum = 10;
OptionGreycolor = 1;    
vLinestyle = {'-','--','-.',':'};
LinestyleExp = '-'; 
FontsizeDefault = 10;
FontsizeAxis = 10;
FontsizeAxisticks = 8;
FontSizeLegend = 8;
Fonttype = 'times';
LinewidthDefault = 1.6;
LinewidthAlt = 1;

if OptionGreycolor == 0
vColors = {[0,42,91]/255,[255 127 0]/255,[128,128,128]/255,[114,47,55]/255};
elseif OptionGreycolor == 1
vColors = {[0.2,0.2,0.2],[0.33 0.33 0.33],[0.46,0.46,0.46],[0.6,0.6,0.6]};
end
ColorZeros = 'k';
StyleZeros = ':';

% Adjust some style options
set(groot, 'DefaultTextInterpreter', 'none');
set(groot, 'defaultAxesTickLabelInterpreter','none');
set(groot, 'defaultLegendInterpreter','none');

%% Run dynare 
%--------------------------------------------------------------------------

% Call dynare
dynare dynareDHFRUncertainty_riskAversion noclearall

if NumModels == 1
mIRF_1 = mIRFProp_zUncertainty_EMAS;
elseif NumModels == 2
mIRF_1 = mIRFProp_zUncertainty_EMAS;
mIRF_2 = mIRFProp_Andreasen_zVol_EMAS; % responses of control variables only!
end

% Put together
if NumModels == 1
aIRF = reshape([mIRF_1],[size(mIRF_1,1),size(mIRF_1,2),NumModels]);
elseif NumModels == 2
aIRF = reshape([mIRF_1,mIRF_2],[size(mIRF_1,1),size(mIRF_1,2),NumModels]);
end

% For rates (u, h, f), switch to ppt instead of percent
uPos = strmatch('u',vNames,'exact');
fPos = strmatch('f',vNames,'exact');
hPos = strmatch('h',vNames,'exact');
aIRF(:,uPos,:) = aIRF(:,uPos,:)*vEMAS(uPos);
aIRF(:,fPos,:) = aIRF(:,fPos,:)*vEMAS(fPos);
aIRF(:,hPos,:) = aIRF(:,hPos,:)*vEMAS(hPos);


%% Plotting via loop
%-------------------------------------------------------------------------

hfig=figure;
for iV = 1:numel(vVariables)
subplot(NumSubplotV,NumSubplotH,iV)
for iN = 1:NumModels
p1=plot(0:IRFPeriodsNum,100*aIRF(1:IRFPeriodsNum+1,strmatch(vVariables{iV},vNames,'exact'),iN),vLinestyle{iN},'LineWidth',LinewidthDefault,'Color',vColors{iN});
hold on
plot(0:IRFPeriodsNum,zeros(IRFPeriodsNum+1,1),StyleZeros,'HandleVisibility','off','Color',ColorZeros,'Linewidth',0.5);
xlim([0 IRFPeriodsNum]);
set(gca,'XTick',[0:2:IRFPeriodsNum],'FontSize',FontsizeAxisticks,'fontname',Fonttype);
ax = gca;
ax.YAxis.Exponent = 0; % avoid scientific notation
box on
end
title(vVNames{iV},'FontSize',FontsizeDefault,'fontname',Fonttype,'FontWeight','normal');

if ismember(iV,[4,5,6])==1
   ylabel('Percentage points','FontSize',FontsizeAxis,'fontname',Fonttype); 
else
 ylabel('Percent','FontSize',FontsizeAxis,'fontname',Fonttype);
end

%if iV >= 3
xlabel('Time (quarters)','FontSize',FontsizeAxis,'fontname',Fonttype);
%end
%end

end

if OptionLegend == 1
subplot(NumSubplotV,NumSubplotH,1)
legend1 = legend(Labels);
set(legend1,'fontname','times','Location','best','FontSize',FontSizeLegend,'interpreter','latex')
end


%% Manual plotting 
%--------------------------------------------------------------------------
%{
hfig=figure;
for iN = 1:NumModels
    
subplot(3,2,1)
p1=plot(0:IRFPeriodsNum,100*aIRF(1:IRFPeriodsNum+1,strmatch('J',vNames,'exact'),iN),vLinestyle{iN},'LineWidth',LinewidthDefault,'Color',vColors{iN});
hold on
plot(0:IRFPeriodsNum,zeros(IRFPeriodsNum+1,1),StyleZeros,'HandleVisibility','off','Color',ColorZeros,'Linewidth',0.5);
xlim([0 IRFPeriodsNum]);
ylim([-0.125 0.125]);
set(gca,'XTick',[0:2:IRFPeriodsNum],'FontSize',FontsizeAxisticks,'fontname',Fonttype);
ax = gca;
ax.YAxis.Exponent = 0; % avoid scientific notation
box on
title('Firm value','FontSize',FontsizeDefault,'fontname',Fonttype,'FontWeight','normal');
ylabel('Percent deviation','FontSize',FontsizeAxis,'fontname',Fonttype);

subplot(3,2,2)
p1=plot(0:IRFPeriodsNum,100*aIRF(1:IRFPeriodsNum+1,strmatch('w',vNames,'exact'),iN),vLinestyle{iN},'LineWidth',LinewidthDefault,'Color',vColors{iN});
hold on
plot(0:IRFPeriodsNum,zeros(IRFPeriodsNum+1,1),StyleZeros,'HandleVisibility','off','Color',ColorZeros,'Linewidth',0.5);
xlim([0 IRFPeriodsNum]);
ylim([-0.00575 0.00575]);
set(gca,'XTick',[0:2:IRFPeriodsNum],'FontSize',FontsizeAxisticks,'fontname',Fonttype);
ax = gca;
ax.YAxis.Exponent = 0; % avoid scientific notation
box on
title('Wage','FontSize',FontsizeDefault,'fontname',Fonttype,'FontWeight','normal');
ylabel('Percent deviation','FontSize',FontsizeAxis,'fontname',Fonttype);

subplot(3,2,3)
p1=plot(0:IRFPeriodsNum,100*aIRF(1:IRFPeriodsNum+1,strmatch('theta',vNames,'exact'),iN),vLinestyle{iN},'LineWidth',LinewidthDefault,'Color',vColors{iN});
hold on
plot(0:IRFPeriodsNum,zeros(IRFPeriodsNum+1,1),StyleZeros,'HandleVisibility','off','Color',ColorZeros,'Linewidth',0.5);
xlim([0 IRFPeriodsNum]);
ylim([-0.25 0.25]);
set(gca,'XTick',[0:2:IRFPeriodsNum],'FontSize',FontsizeAxisticks,'fontname',Fonttype);
ax = gca;
ax.YAxis.Exponent = 0; % avoid scientific notation
box on
title('Labor market tightness','FontSize',FontsizeDefault,'fontname',Fonttype,'FontWeight','normal');
ylabel('Percent deviation','FontSize',FontsizeAxis,'fontname',Fonttype);

subplot(3,2,4)
p1=plot(0:IRFPeriodsNum,100*aIRF(1:IRFPeriodsNum+1,strmatch('h',vNames,'exact'),iN),vLinestyle{iN},'LineWidth',LinewidthDefault,'Color',vColors{iN});
hold on
plot(0:IRFPeriodsNum,zeros(IRFPeriodsNum+1,1),StyleZeros,'HandleVisibility','off','Color',ColorZeros,'Linewidth',0.5);
xlim([0 IRFPeriodsNum]);
ylim([-0.25 0.25]);
%ylim([-0.125 0.125]);
%ylim([-0.025 0.025]);
set(gca,'XTick',[0:2:IRFPeriodsNum],'FontSize',FontsizeAxisticks,'fontname',Fonttype);
ax = gca;
ax.YAxis.Exponent = 0; % avoid scientific notation
box on
title('Hiring rate','FontSize',FontsizeDefault,'fontname',Fonttype,'FontWeight','normal');
ylabel('Percentage point deviation','FontSize',FontsizeAxis,'fontname',Fonttype);

subplot(3,2,5)
p1=plot(0:IRFPeriodsNum,100*aIRF(1:IRFPeriodsNum+1,strmatch('f',vNames,'exact'),iN),vLinestyle{iN},'LineWidth',LinewidthDefault,'Color',vColors{iN});
hold on
plot(0:IRFPeriodsNum,zeros(IRFPeriodsNum+1,1),StyleZeros,'HandleVisibility','off','Color',ColorZeros,'Linewidth',0.5);
xlim([0 IRFPeriodsNum]);
ylim([-0.125 0.125]);
set(gca,'XTick',[0:2:IRFPeriodsNum],'FontSize',FontsizeAxisticks,'fontname',Fonttype);
ax = gca;
ax.YAxis.Exponent = 0; % avoid scientific notation
box on
title('Job finding rate','FontSize',FontsizeDefault,'fontname',Fonttype,'FontWeight','normal');
ylabel('Percentage point deviation','FontSize',FontsizeAxis,'fontname',Fonttype);
xlabel('Time (quarters)','FontSize',FontsizeAxis,'fontname',Fonttype);

subplot(3,2,6)
p1=plot(0:IRFPeriodsNum,100*aIRF(1:IRFPeriodsNum+1,strmatch('u',vNames,'exact'),iN),vLinestyle{iN},'LineWidth',LinewidthDefault,'Color',vColors{iN});
hold on
plot(0:IRFPeriodsNum,zeros(IRFPeriodsNum+1,1),StyleZeros,'HandleVisibility','off','Color',ColorZeros,'Linewidth',0.5);
xlim([0 IRFPeriodsNum]);
ylim([-0.055 0.055]);
%ylim([-0.025 0.025]);
set(gca,'XTick',[0:2:IRFPeriodsNum],'FontSize',FontsizeAxisticks,'fontname',Fonttype);
ax = gca;
ax.YAxis.Exponent = 0; % avoid scientific notation
box on
title('Unemployment rate','FontSize',FontsizeDefault,'fontname',Fonttype,'FontWeight','normal');
ylabel('Percentage point deviation','FontSize',FontsizeAxis,'fontname',Fonttype);
xlabel('Time (quarters)','FontSize',FontsizeAxis,'fontname',Fonttype);

end

if OptionLegend == 1
subplot(3,2,1)
legend1 = legend(Labels);
set(legend1,'fontname','times','Location','best','FontSize',FontSizeLegend,'interpreter','latex')
end
%}

%% Print 
%--------------------------------------------------------------------------
xSize = 2*8.75; 
ySize = 3*6.25; 
xCut = 2;
yCut = 0.5;

set(gcf,'Units','centimeters','Position',[0 0 xSize ySize],'PaperUnits','centimeters' ...
      ,'PaperPosition',[0 0 xSize ySize],'PaperSize',[xSize-xCut ySize-yCut],'PaperPositionMode','auto')

if OptionPrint == 1      
    if OptionGreycolor == 0
         FigNamepdf =horzcat(horzcat(TargetPath,horzcat(FigName,'_CS')),'.pdf');
    else
         FigNamepdf =horzcat(horzcat(TargetPath,FigName),'.pdf');
    end
  print(hfig,FigNamepdf,'-dpdf','-painters')
end

%% Done!
%----------------------------------------------------------------------------
TimeEnd = toc(TimeStart);
disp(['Total run time was ',num2str(TimeEnd),' seconds']);

%% 'Appendix'
%----------------------------------------------------------------------------
%% Remark on the (G)IRF for productivity level
% Why did I use "z_Obs", "sigma_z_Obs" and "n_Obs" above? Purely out of laziness and to "trick"
% Dynare into believing that these are control rather than state variables. 
% Else you need to read in a separate Andreasen-method-generated IRFs for state variables.
% Can compare this to:
% Test1 = mIRFProp_Andreasen_zVol_state(:,2); % positions; outDynare.label_v
% Test2 = mIRFProp_Andreasen_zVol_EMAS(:,strmatch('z_Obs',vNames))
