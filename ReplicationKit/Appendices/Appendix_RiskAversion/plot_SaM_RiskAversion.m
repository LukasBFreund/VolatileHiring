%%=========================================================================
% Replication code for den Haan, Freund and Rendahl (2021)
%
% This file: plot results for baseline SaM model, under risk aversion
%
% Run on Matlab R2019b, Dynare 4.4.3
% Last updated: November 2020
% For any questions please email lukas.beat.freund@gmail.com
%%=========================================================================

%% Housekeeping
%--------------------------------------------------------------------------
clear; close all; clc;

%% User Options
%--------------------------------------------------------------------------
% Select from figures
figureName = 'fig_SaM_RiskAversion_WageLin';
% 'fig_SaM_RiskAversion_WageLin'
% 'fig_SaM_RiskAversion_WageComparison'

% Choose whether to print and, if so, set figure name and target folder to save figure in
sPlotting.options.print = 1;        
sPlotting.options.targetPath ='.\Output\Figures\';

% Load general design settings
load inputSettings

%% Load files and recover IRFs
%--------------------------------------------------------------------------

switch figureName
%--------------------------------------------------------------------------   
    case 'fig_SaM_RiskAversion_WageLin'
    sPlotting.options.figName = figureName;
    sPlotting.numModels = 2; 

    % Baseline case with risk aversion
    load(fullfile('.', 'Output\IRFs\', 'IRFs_LogU_WageLin_Baseline')); mIRFProp_1=mIRFProp_zUncertainty_EMAS; 
    uPos = strmatch('u',vNames,'exact'); % Show unemployment in ppt, rest in percent
    mIRFProp_1(:,uPos,:) = mIRFProp_1(:,uPos,:)*vEMAS(uPos);

    % Imposing risk-neutral pricing
    load(fullfile('.', 'Output\IRFs\', 'IRFs_LogU_WageLin_RiskNeutralPricing')); mIRFProp_2=mIRFProp_zUncertainty_EMAS; 
    mIRFProp_2(:,uPos,:) = mIRFProp_2(:,uPos,:)*vEMAS(uPos);

    sPlotting.vVariables = {'J','w','u','c'};
    sPlotting.vVNames = {'Match value', 'Wage','Unemployment rate','Consumption'};
    sPlotting.numSubplotV = 2;
    sPlotting.numSubplotH = 2;
    sPlotting.legend.show = 1;
    sPlotting.legend.labels{1} = 'Log utility';
    sPlotting.legend.labels{2} = 'Risk-neutral pricing';
    sPlotting.legend.position = 'southeast';

%--------------------------------------------------------------------------   
    case 'fig_SaM_RiskAversion_WageComparison'
    sPlotting.options.figName = figureName;
    sPlotting.numModels = 4; 

    load(fullfile('.', 'Output\IRFs\', 'IRFs_LogU_WageLin_Baseline')); mIRFProp_1=mIRFProp_zUncertainty_EMAS; 
    uPos = strmatch('u',vNames,'exact'); 
    mIRFProp_1(:,uPos,:) = mIRFProp_1(:,uPos,:)*vEMAS(uPos);

    load(fullfile('.', 'Output\IRFs\', 'IRFs_LogU_WageNash_Baseline')); mIRFProp_2=mIRFProp_zUncertainty_EMAS; 
    mIRFProp_2(:,uPos,:) = mIRFProp_2(:,uPos,:)*vEMAS(uPos);

    load(fullfile('.', 'Output\IRFs\', 'IRFs_LogU_WageNash_RiskNeutralWage')); mIRFProp_3=mIRFProp_zUncertainty_EMAS; 
    mIRFProp_3(:,uPos,:) = mIRFProp_3(:,uPos,:)*vEMAS(uPos);

    load(fullfile('.', 'Output\IRFs\', 'IRFs_LogU_WageNash_MUConstantInNash')); mIRFProp_4=mIRFProp_zUncertainty_EMAS; 
    mIRFProp_4(:,uPos,:) = mIRFProp_4(:,uPos,:)*vEMAS(uPos);

    
    sPlotting.vVariables = {'J','w','u','c'};
    sPlotting.vVNames = {'Match value', 'Wage','Unemployment rate','Consumption'};
    sPlotting.numSubplotV = 2;
    sPlotting.numSubplotH = 2;
    sPlotting.legend.show = 1;
    sPlotting.legend.labels{1} = 'Linear wage';
    sPlotting.legend.labels{2} = 'Nash wage';
    sPlotting.legend.labels{3} = 'Risk-neutral Nash wage';
    sPlotting.legend.labels{4} = 'Nash wage, constant leisure value';
    sPlotting.legend.position = 'southeast';

end
 
 % Some common settings
sPlotting.IRFPeriods = 10;
sPlotting.options.color = 0;

 %% 3. Manipulate
%--------------------------------------------------------------------------
          
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

if sPlotting.legend.show == 1
subplot(sPlotting.numSubplotV,sPlotting.numSubplotH,1)
legend1 = legend(sPlotting.legend.labels);
set(legend1,'fontname',sSettings.font.name,'Location',sPlotting.legend.position,'FontSize',sSettings.font.size.legend,'interpreter','none')
end

%% Print 
%--------------------------------------------------------------------------

numPlots = sPlotting.numSubplotV*sPlotting.numSubplotH;

switch numPlots
    case 1
    xSize = 17.5/2; 
    ySize = 6.25;
    xCut = -0.5;
    yCut = 0;
 case 2
    xSize = 17.5; 
    ySize = 6.25;
    xCut = 1;
    yCut = 0;
    case 4
    xSize = 17.5; 
    ySize = 12.5;
    xCut = 1.5;
    yCut = 0.5;
    case 6
    xSize = 17.5; 
    ySize = 12.5+6.25; 
    xCut = 2;
    yCut = 0.5;
    case 8
    xSize = 17.5; 
    ySize = 12.5;
    xCut = 1.8;
    yCut = 0.8;
end

set(gcf,'Units','centimeters','Position',[0 0 xSize ySize],'PaperUnits','centimeters' ...
      ,'PaperPosition',[0 0 xSize ySize],'PaperSize',[xSize-xCut ySize-yCut],'PaperPositionMode','auto')

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
clear ax iN iV uPos vEMAS vNames vNames_Andreasen legend1 mIRFProp_Andreasen_zVol_EMAS mIRFProp_zUncertainty_EMAS mIRFProp_1 mIRFProp_2 numPlots xCut xSize yCut ySize