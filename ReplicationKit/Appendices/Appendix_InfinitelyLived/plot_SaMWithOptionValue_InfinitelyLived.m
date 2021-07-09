%%=========================================================================
% Replication code for den Haan, Freund and Rendahl (2021)
%
% This file: compare plots for baseline model vs. alternative
% with infinitely-lived entrepreneurs.

% Run on Matlab R2019b, Dynare 4.4.3
% Last updated: July 2021
% For any questions please email lukas.beat.freund@gmail.com
%%=========================================================================

%% Housekeeping
%--------------------------------------------------------------------------
clear; close all; clc;

%% User options
%--------------------------------------------------------------------------
% Choose whether to print 
sPlotting.options.print = 1;        

% Load general design settings
load inputSettings

%% What variables to plot
sPlotting.vVariables = {'JaHat_aHatSS','JU_d','p','u'}; 
sPlotting.vVNames = {'Value of entrepreneur at steady-state cutoff level','Value of unmatched entrepreneur',...
         'Entry probability','Unemployment rate'};
sPlotting.IRFPeriods = 10;
sPlotting.numSubplotV = 2;
sPlotting.numSubplotH = 2;
sPlotting.legend.show = 1;

%% Load data from different models
%--------------------------------------------------------------------------
sPlotting.numModels = 3;
sPlotting.legend.labels{1} = 'Baseline';
sPlotting.legend.labels{2} = 'Infinitely-lived, $\gamma = 0$';
sPlotting.legend.labels{3} = 'Infinitely-lived, $\gamma \rightarrow 1$';
sPlotting.legend.position = 'northeast';

%% Model 1 - Baseline: 'death' (long-lived); 0.003 (~max)
load(fullfile('.', 'Inputs\', 'IRFs_Baseline_sigmaa003_recalib'));
mIRFProp_1=mIRFProp_zUncertainty_EMAS; 

% Adjust units on y axis for select variables
uPos = strmatch('u',vNames,'exact');
hPos = strmatch('h',vNames,'exact');
pPos = strmatch('p',vNames,'exact');
JU_dPos = strmatch('JU_d',vNames,'exact'); % and here use abs value
JaHat_aHatSSPos = strmatch('JaHat_aHatSS',vNames,'exact');

mIRFProp_1(:,uPos,:) = mIRFProp_1(:,uPos,:)*vEMAS(uPos);
mIRFProp_1(:,hPos,:) = mIRFProp_1(:,hPos,:)*vEMAS(hPos);
mIRFProp_1(:,pPos,:) = mIRFProp_1(:,pPos,:)*vEMAS(pPos);
mIRFProp_1(:,JU_dPos,:) = mIRFProp_1(:,JU_dPos,:)*vEMAS(JU_dPos);
mIRFProp_1(:,JaHat_aHatSSPos,:) =mIRFProp_1(:,JaHat_aHatSSPos,:)*vEMAS(JaHat_aHatSSPos);

%% Model 2 - 'no death' (infinitely-lived); 0.025 (~max)
load(fullfile('.', 'Inputs\', 'IRFs_Gamma0Model_sigmaa025_recalib')); 
mIRFProp_2=mIRFProp_zUncertainty_EMAS; 
mIRFProp_2 = mIRFProp_2(:,1:24); % has more endo variables

% Adjust units on y axis for select variables
uPos = strmatch('u',vNames,'exact');
hPos = strmatch('h',vNames,'exact');
pPos = strmatch('p',vNames,'exact');
JU_dPos = strmatch('JU_d',vNames,'exact'); % and here use abs value
JaHat_aHatSSPos = strmatch('JaHat_aHatSS',vNames,'exact');
mIRFProp_2(:,uPos,:) = mIRFProp_2(:,uPos,:)*vEMAS(uPos);
mIRFProp_2(:,hPos,:) = mIRFProp_2(:,hPos,:)*vEMAS(hPos);
mIRFProp_2(:,pPos,:) = mIRFProp_2(:,pPos,:)*vEMAS(pPos);
mIRFProp_2(:,JU_dPos,:) = mIRFProp_2(:,JU_dPos,:)*vEMAS(JU_dPos);
mIRFProp_2(:,JaHat_aHatSSPos,:) =mIRFProp_2(:,JaHat_aHatSSPos,:)*vEMAS(JaHat_aHatSSPos);

%% Model 3 - gamma = 0.999, sigma_a = 0.003, death calib
load(fullfile('.', 'Inputs\', 'IRFs_GammaModel_sigmaa003_BaselineCalib_gamma0999')); 
mIRFProp_3=mIRFProp_zUncertainty_EMAS; 
mIRFProp_3 = mIRFProp_3(:,1:24); % has more endo variables

% Adjust units on y axis for select variables
uPos = strmatch('u',vNames,'exact');
hPos = strmatch('h',vNames,'exact');
pPos = strmatch('p',vNames,'exact');
JU_dPos = strmatch('JU_d',vNames,'exact'); % and here use abs value
JaHat_aHatSSPos = strmatch('JaHat_aHatSS',vNames,'exact');
mIRFProp_3(:,uPos,:) = mIRFProp_3(:,uPos,:)*vEMAS(uPos);
mIRFProp_3(:,hPos,:) = mIRFProp_3(:,hPos,:)*vEMAS(hPos);
mIRFProp_3(:,pPos,:) = mIRFProp_3(:,pPos,:)*vEMAS(pPos);
mIRFProp_3(:,JU_dPos,:) = mIRFProp_3(:,JU_dPos,:)*vEMAS(JU_dPos);
mIRFProp_3(:,JaHat_aHatSSPos,:) =mIRFProp_3(:,JaHat_aHatSSPos,:)*vEMAS(JaHat_aHatSSPos);

%% Put together
if sPlotting.numModels == 1
sResults.aIRF = reshape([mIRFProp_1],[size(mIRFProp_1,1),size(mIRFProp_1,2),sPlotting.numModels]);
elseif sPlotting.numModels == 2
sResults.aIRF = reshape([mIRFProp_1,mIRFProp_2],[size(mIRFProp_1,1),size(mIRFProp_1,2),sPlotting.numModels]);
elseif sPlotting.numModels == 3
sResults.aIRF = reshape([mIRFProp_1,mIRFProp_2,mIRFProp_3],[size(mIRFProp_1,1),size(mIRFProp_1,2),sPlotting.numModels]);
elseif sPlotting.numModels == 4
sResults.aIRF = reshape([mIRFProp_1,mIRFProp_2,mIRFProp_3,mIRFProp_4],[size(mIRFProp_1,1),size(mIRFProp_1,2),sPlotting.numModels]);
end

% Adjust values that are virtually zero (not exactly due to use of a solver) to exactly zero to avoid confusing graphs
sResults.aIRF(abs(sResults.aIRF)<1e-10)=0;
 
%% Plot 
fig=figure;
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
    
    if ismember(sPlotting.vVariables(iV),{'u','f','h','p'}) == 1
        ylabel('Deviation (ppts.)','FontSize',sSettings.font.size.axis,'fontname',sSettings.font.name);
    elseif ismember(sPlotting.vVariables(iV),{'aHat','JU_d','JaHat_aHatSS','JU1'}) == 1
        ylabel('Deviation (abs.)','FontSize',sSettings.font.size.axis,'fontname',sSettings.font.name);
    else        
       ylabel('Deviation (pct.)','FontSize',sSettings.font.size.axis,'fontname',sSettings.font.name);
    end    
    if iV >= (sPlotting.numSubplotV*sPlotting.numSubplotH)-1
       xlabel('Time (quarters)','FontSize',sSettings.font.size.axis,'fontname',sSettings.font.name);
    end
    if ismember(sPlotting.vVariables(iV),{'JU_d'}) == 1
          ylim([0 0.05])
       end
end

if sPlotting.legend.show == 1
subplot(sPlotting.numSubplotV,sPlotting.numSubplotH,1)
legend1 = legend(sPlotting.legend.labels);
set(legend1,'fontname',sSettings.font.name,'Location',sPlotting.legend.position,'FontSize',sSettings.font.size.legend,'interpreter','latex')
end

%% Print 
sSettings.plots.xSize = 2*8.75; sSettings.plots.ySize = 2*6.25;  sSettings.plots.xCut = 1; sSettings.plots.yCut = 0.5;

set(gcf,'Units','centimeters','Position',[0 0 sSettings.plots.xSize sSettings.plots.ySize],'PaperUnits','centimeters' ...
      ,'PaperPosition',[0 0 sSettings.plots.xSize sSettings.plots.ySize],'PaperSize',[sSettings.plots.xSize-sSettings.plots.xCut sSettings.plots.ySize-sSettings.plots.yCut],'PaperPositionMode','auto')
  
if sPlotting.options.print == 1      
    print(fig,'.\Output\Fig_App_HetFirm_Comparison_InfinitelyLived','-dpdf','-painters')
end
