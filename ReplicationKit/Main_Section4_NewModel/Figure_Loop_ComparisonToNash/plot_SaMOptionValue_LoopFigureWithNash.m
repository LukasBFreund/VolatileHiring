%%=========================================================================
% Replication code for den Haan, Freund and Rendahl (2021)
%
% This file: compare uncertainty effects in the new model, looping
% over sigma_a, to those in the baseline model but with Nash bargaining

% Run on Matlab R2019b, Dynare 4.4.3
% Last updated: October 2020

% For any questions, please email lukas.beat.freund@gmail.com
%%=========================================================================

%% Housekeeping
%--------------------------------------------------------------------------
clear; close all; clc;
load inputSettings

%% Load data
%--------------------------------------------------------------------------

% New model, with recalibration
load('Inputs\resultsOptionValue_p05_Recalib_sigmaa001To003')
vSigma_a = sResults.sigmaass;
uPos_1 = strmatch('u',sResults.names,'exact');
mResults = sResults.aStore(1,:,1,:);  % impact effect
mResults_Recalib = reshape(mResults,size(mResults,2),numel(vSigma_a)); %IRFs
sResults_Recalib = sResults;
mResults_TV = sResults.aStore(1,:,2,:);  
mResults_TV = reshape(mResults_TV,size(mResults_TV,2),numel(vSigma_a));
uIRF_TV_Recalib = sResults.aStore(:,uPos_1,2,:); 
uIRF_TV_Recalib = reshape(uIRF_TV_Recalib,20,20);
uIRF_TV_Recalib_Max = max(uIRF_TV_Recalib); 
clear mResults sResults

% Comparison: baseline model with Nash
load('Inputs\results_Nash_Baseline')
mResults_Nash_PU = sResults.aIRF(:,:,1);
mResults_Nash_TV = sResults.aIRF(:,:,2);
uIRF_Nash_Impact = mResults_Nash_PU(1,strmatch('u',sResults.vNames,'exact'));
uIRF_Nash_MaxTV = max(mResults_Nash_TV(:,strmatch('u',sResults.vNames,'exact')));

%% Plotting
%--------------------------------------------------------------------------

fig=figure;
subplot(1,2,1)
p1=plot(vSigma_a,100*mResults_Recalib(uPos_1,:),'LineWidth',sSettings.lines.width,'linestyle',sSettings.lines.list{1},'Color',sSettings.colors.list{1});
hold on 
p2 = plot(vSigma_a,ones(1,numel(vSigma_a))*100*uIRF_Nash_Impact,'LineWidth',sSettings.lines.width,'linestyle',sSettings.lines.list{2},'Color',sSettings.colors.list{2});
xlabel('$\sigma_a$','FontSize',sSettings.font.size.axis,'fontname',sSettings.font.name,'interpreter','latex');
set(gca,'XTick',[vSigma_a(1):((vSigma_a(end)-vSigma_a(1))/2):vSigma_a(end)],'FontSize',sSettings.font.size.axisticks,'fontname',sSettings.font.name);
ylabel('Effect on unemp. rate (ppts.)','FontSize',sSettings.font.size.axis,'fontname',sSettings.font.name,'interpreter','none');
box on
title('Impact effect','FontSize',sSettings.font.size.default,'fontname',sSettings.font.name,'FontWeight','normal');
ax = gca;
ax.XAxis.Exponent = 0; 
ax.YAxis.Exponent = 0; 

subplot(1,2,2)
p1=plot(vSigma_a,100*uIRF_TV_Recalib_Max,'LineWidth',sSettings.lines.width,'linestyle',sSettings.lines.list{1},'Color',sSettings.colors.list{1});
hold on
p2 = plot(vSigma_a,ones(1,numel(vSigma_a))*100*uIRF_Nash_MaxTV,'LineWidth',sSettings.lines.width,'linestyle',sSettings.lines.list{2},'Color',sSettings.colors.list{2});
xlabel('$\sigma_a$','FontSize',sSettings.font.size.axis,'fontname',sSettings.font.name,'interpreter','latex');
set(gca,'XTick',[vSigma_a(1):((vSigma_a(end)-vSigma_a(1))/2):vSigma_a(end)],'FontSize',sSettings.font.size.axisticks,'fontname',sSettings.font.name);
ylabel('Effect on unemp. rate (ppts.)','FontSize',sSettings.font.size.axis,'fontname',sSettings.font.name,'interpreter','none');
box on
title('Maximum total volatility effect ','FontSize',sSettings.font.size.default,'fontname',sSettings.font.name,'FontWeight','normal');
ax = gca;
ax.XAxis.Exponent = 0; 
ax.YAxis.Exponent = 0; 

% Print
xSize = 2*8.75; ySize = 1*6.25;  xCut = 1; yCut = 0;

set(gcf,'Units','centimeters','Position',[0 0 xSize ySize],'PaperUnits','centimeters' ...
      ,'PaperPosition',[0 0 xSize ySize],'PaperSize',[xSize-xCut ySize-yCut],'PaperPositionMode','auto')

print(fig,'Output\fig_SaM_OptionValue_ComparisonToNash','-dpdf','-painters')
