%%=========================================================================
% Replication code for den Haan, Freund and Rendahl (2021)
%
% This file: plot results for extended version of heterogeneous-firm
% version that allows for endogenous Upsilon

% Run on Matlab R2019b, Dynare 4.4.3
% Last updated: July 2021
% For any questions please email lukas.beat.freund@gmail.com
%%=========================================================================

%% Housekeeping
%--------------------------------------------------------------------------
clear; close all; clc;

%% User Options
%--------------------------------------------------------------------------

% Choose whether to print 
sPlotting.options.print = 1;        
sPlotting.options.targetPath ='.\Output\';

% Load general design settings
load inputSettings

%% What to plot and how
sPlotting.IRFPeriods = 10;
sPlotting.vVariables = {'JU_d','u','g','Upsilon'};
sPlotting.vVNames = {'Value of unmatched entrepreneur','Unemployment rate', 'Entry probability: entrepreneurs','Mass of potential entrepreneurs'};
sPlotting.numSubplotV = 2;
sPlotting.numSubplotH = 2;
sPlotting.legend.show = 1;
sPlotting.legend.labels{3} = '$\iota = 0$';
sPlotting.legend.labels{2} = '$\iota = 10$';
sPlotting.legend.labels{1} = '$\iota = 100$';
sPlotting.legend.position = 'northeast';

%% Load files and recover IRFs
%--------------------------------------------------------------------------
sPlotting.numModels = 3;

%% iota = 100
load(fullfile('.', 'Inputs\', 'IRFs_SaMOptionValue_Uniform_EndoUps_isoelastic_iota100')); 
mIRFProp_1=mIRFProp_zUncertainty_EMAS; 

% For rates (u, h, f, p), switch to ppt 
uPos = strmatch('u',vNames,'exact');
hPos = strmatch('h',vNames,'exact');
pPos = strmatch('p',vNames,'exact');
gPos = strmatch('g',vNames,'exact');
JU_dPos = strmatch('JU_d',vNames,'exact'); % and here use abs value

mIRFProp_1(:,uPos,:) = mIRFProp_1(:,uPos,:)*vEMAS(uPos);
mIRFProp_1(:,hPos,:) = mIRFProp_1(:,hPos,:)*vEMAS(hPos);
mIRFProp_1(:,pPos,:) = mIRFProp_1(:,pPos,:)*vEMAS(pPos);
mIRFProp_1(:,gPos,:) = mIRFProp_1(:,gPos,:)*vEMAS(gPos);
mIRFProp_1(:,JU_dPos,:) = mIRFProp_1(:,JU_dPos,:)*vEMAS(JU_dPos);

%% iota = 10
load(fullfile('.', 'Inputs\', 'IRFs_SaMOptionValue_Uniform_EndoUps_isoelastic_iota10')); 
mIRFProp_2=mIRFProp_zUncertainty_EMAS; 

% Transform as well
uPos = strmatch('u',vNames,'exact');
hPos = strmatch('h',vNames,'exact');
pPos = strmatch('p',vNames,'exact');
gPos = strmatch('g',vNames,'exact');
JU_dPos = strmatch('JU_d',vNames,'exact'); % and here use abs value
mIRFProp_2(:,uPos,:) = mIRFProp_2(:,uPos,:)*vEMAS(uPos);
mIRFProp_2(:,hPos,:) = mIRFProp_2(:,hPos,:)*vEMAS(hPos);
mIRFProp_2(:,pPos,:) = mIRFProp_2(:,pPos,:)*vEMAS(pPos);
mIRFProp_2(:,gPos,:) = mIRFProp_2(:,gPos,:)*vEMAS(gPos);
mIRFProp_2(:,JU_dPos,:) = mIRFProp_2(:,JU_dPos,:)*vEMAS(JU_dPos);

%% iota = 0
load(fullfile('.', 'Inputs\', 'IRFs_SaMOptionValue_Uniform_EndoUps_isoelastic_iota0')); 
mIRFProp_3=mIRFProp_zUncertainty_EMAS; 

% Transform as well
uPos = strmatch('u',vNames,'exact');
hPos = strmatch('h',vNames,'exact');
pPos = strmatch('p',vNames,'exact');
gPos = strmatch('g',vNames,'exact');
JU_dPos = strmatch('JU_d',vNames,'exact'); % and here use abs value
mIRFProp_3(:,uPos,:) = mIRFProp_3(:,uPos,:)*vEMAS(uPos);
mIRFProp_3(:,hPos,:) = mIRFProp_3(:,hPos,:)*vEMAS(hPos);
mIRFProp_3(:,pPos,:) = mIRFProp_3(:,pPos,:)*vEMAS(pPos);
mIRFProp_3(:,gPos,:) = mIRFProp_3(:,gPos,:)*vEMAS(gPos);
mIRFProp_3(:,JU_dPos,:) = mIRFProp_3(:,JU_dPos,:)*vEMAS(JU_dPos);

% Put together
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
%--------------------------------------------------------------------------
fig1=figure;
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
    
    if ismember(sPlotting.vVariables(iV),{'u','f','h','p','g'}) == 1
        ylabel('Deviation (ppts.)','FontSize',sSettings.font.size.axis,'fontname',sSettings.font.name);
    elseif ismember(sPlotting.vVariables(iV),{'aHat','JU_d','JaHat_aHatSS','JU1'}) == 1
        ylabel('Deviation (abs.)','FontSize',sSettings.font.size.axis,'fontname',sSettings.font.name);
    else        
       ylabel('Deviation (pct.)','FontSize',sSettings.font.size.axis,'fontname',sSettings.font.name);
    end    
    if iV >= (sPlotting.numSubplotV*sPlotting.numSubplotH)-1
       xlabel('Time (quarters)','FontSize',sSettings.font.size.axis,'fontname',sSettings.font.name);
    end
end

if sPlotting.legend.show == 1
subplot(sPlotting.numSubplotV,sPlotting.numSubplotH,1)
legend1 = legend(sPlotting.legend.labels);
set(legend1,'fontname',sSettings.font.name,'Location',sPlotting.legend.position,'FontSize',sSettings.font.size.legend,'interpreter','latex')
end

%% Print 
%--------------------------------------------------------------------------

xSize = 2*8.75; ySize = 2*6.25;  xCut = 1; yCut = 0.5;
set(gcf,'Units','centimeters','Position',[0 0 xSize ySize],'PaperUnits','centimeters' ...
      ,'PaperPosition',[0 0 xSize ySize],'PaperSize',[xSize-xCut ySize-yCut],'PaperPositionMode','auto')

  
if sPlotting.options.print == 1      
    print(fig1,'Output\fig_App_SaM_HetFirm_EndoMass_Isoelastic_Comparison','-dpdf','-painters')
end
