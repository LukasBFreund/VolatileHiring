%%=========================================================================
% Replication code for den Haan, Freund and Rendahl (2021)
%
% SaM model that features option-value effects due to dispersion in
% idiosyncratic productivity draws and a finite mass of potential 
% entrepreneurs.

% This file: plot 2x2 figure that loops over sigma_a values 
% without recalibration

% Run on Matlab R2019b, Dynare 4.4.3
% Last updated: October 2020

% For any questions, please email lukas.beat.freund@gmail.com
%%=========================================================================

%% Housekeeping
%--------------------------------------------------------------------------
clear; close all; clc;
TimeStart = tic;
load inputSettings

%% Load data
%--------------------------------------------------------------------------

% No recalibration, variable h
load('Inputs\resultsOptionValue_p05_NoRecalib_sigmaa0To001_Fine')
vSigma_a = sResults.sigmaass;
uPos = strmatch('u',sResults.names,'exact');
palpha = sResults.sPar(1,1).alpha;
mResults = sResults.aStore(1,:,1,:);  % IRFs
mResults_1 = reshape(mResults,size(mResults,2),numel(vSigma_a));
sResults_1 = sResults;

% No recalibration, fixed h
load('Inputs\resultsOptionValue_p05_NoRecalib_sigmaa00005To001_FixedH_Fine') % careful, has a different starting point
mResults = sResults.aStore(1,:,1,:); 
vSigma_a_fixed = sResults.sigmaass;
mResults_2 = reshape(mResults,size(mResults,2),numel(vSigma_a_fixed)); %IRFs
sResults_2 = sResults;
clear mResults sResults


%% Plotting
%--------------------------------------------------------------------------

fig=figure;
subplot(1,2,1)
p1=plot(vSigma_a,100*sResults_1.vuss,'LineWidth',sSettings.lines.width,'linestyle',sSettings.lines.list{1},'Color',sSettings.colors.list{1});
xlabel('$\sigma_a$','FontSize',sSettings.font.size.axis,'fontname',sSettings.font.name,'interpreter','latex');
set(gca,'XTick',[vSigma_a(1):((vSigma_a(end)-vSigma_a(1))/2):vSigma_a(end)],'FontSize',sSettings.font.size.axisticks,'fontname',sSettings.font.name);
ylabel('Percent','FontSize',sSettings.font.size.axis,'fontname',sSettings.font.name);
box on
title('Steady-state unemployment rate','FontSize',sSettings.font.size.default,'fontname',sSettings.font.name,'FontWeight','normal');
ax = gca;
ax.XAxis.Exponent = 0; 

subplot(1,2,2)
p1=plot(vSigma_a,-(1/palpha)*sResults_1.vElasticity,'LineWidth',sSettings.lines.width,'linestyle',sSettings.lines.list{1},'Color',sSettings.colors.list{1});
%  multiplication by -1/alpha is to go from hiring rate to labor market tightness (given alpha = 0.5)
hold on
xlabel('$\sigma_a$','FontSize',sSettings.font.size.axis,'fontname',sSettings.font.name,'interpreter','latex');
set(gca,'XTick',[vSigma_a(1):((vSigma_a(end)-vSigma_a(1))/2):vSigma_a(end)],'FontSize',sSettings.font.size.axisticks,'fontname',sSettings.font.name);
ylabel('Elasticity of tightness','FontSize',sSettings.font.size.axis,'fontname',sSettings.font.name,'interpreter','none');
box on
title('Steady-state change in agg. prod.','FontSize',sSettings.font.size.default,'fontname',sSettings.font.name,'FontWeight','normal');
ax = gca;
ax.XAxis.Exponent = 0; 

% Print
xSize = 2*8.75; ySize = 1*6.25;  xCut = 1; yCut = 0;
set(gcf,'Units','centimeters','Position',[0 0 xSize ySize],'PaperUnits','centimeters' ...
      ,'PaperPosition',[0 0 xSize ySize],'PaperSize',[xSize-xCut ySize-yCut],'PaperPositionMode','auto')
print(fig,'Output\fig_SaMOptionValue_NoRecalib_Properties_1x2','-dpdf','-painters')

%% Alternative figure that shows additional subplots 
% Bottom row of panels here shows the impact effect of an uncertainty
% shock with a variable hiring rate vs. fixed, keeping 
% the parameters at the 'unrecalibrated' levels

%{
fig=figure;
subplot(2,2,1)
p1=plot(vSigma_a,100*sResults_1.vuss,'LineWidth',sSettings.lines.width,'linestyle',sSettings.lines.list{1},'Color',sSettings.colors.list{1});
xlabel('$\sigma_a$','FontSize',sSettings.font.size.axis,'fontname',sSettings.font.name,'interpreter','latex');
set(gca,'XTick',[vSigma_a(1):((vSigma_a(end)-vSigma_a(1))/2):vSigma_a(end)],'FontSize',sSettings.font.size.axisticks,'fontname',sSettings.font.name);
ylabel('Percent','FontSize',sSettings.font.size.axis,'fontname',sSettings.font.name);
box on
title('Steady-state unemployment rate','FontSize',sSettings.font.size.default,'fontname',sSettings.font.name,'FontWeight','normal');
ax = gca;
ax.XAxis.Exponent = 0; 

subplot(2,2,2)
p1=plot(vSigma_a,-(1/palpha)*sResults_1.vElasticity,'LineWidth',sSettings.lines.width,'linestyle',sSettings.lines.list{1},'Color',sSettings.colors.list{1});
%  multiplication by -1/alpha is to go from hiring rate to labor market tightness (given alpha = 0.5)
hold on
xlabel('$\sigma_a$','FontSize',sSettings.font.size.axis,'fontname',sSettings.font.name,'interpreter','latex');
set(gca,'XTick',[vSigma_a(1):((vSigma_a(end)-vSigma_a(1))/2):vSigma_a(end)],'FontSize',sSettings.font.size.axisticks,'fontname',sSettings.font.name);
ylabel('Elasticity of tightness','FontSize',sSettings.font.size.axis,'fontname',sSettings.font.name,'interpreter','none');
box on
title('Steady-state change in agg. prod.','FontSize',sSettings.font.size.default,'fontname',sSettings.font.name,'FontWeight','normal');
ax = gca;
ax.XAxis.Exponent = 0; 

subplot(2,2,3)
% To keep the plot clean, we dont show sigma_a=0
% But one could compute manually the sigma_a case
% uss = sResults_2.sPar.uBar; uss=uss(:);
% nss = 1-uss;
% uSigmaA0 = 1-(1-0.1)*nss;
% IRFSigmaA0 = uSigmaA0-uss; % we start from steady-state, then there are 0 vacancies posted
%p1=plot(vSigma_a,100*[IRFSigmaA0 mResults_2(uPos,2:end)],'LineWidth',sSettings.lines.width,'linestyle',sSettings.lines.list{1},'Color',sSettings.colors.list{1});
p1=plot(vSigma_a_fixed,100*[mResults_2(uPos,1:end)],'LineWidth',sSettings.lines.width,'linestyle',sSettings.lines.list{1},'Color',sSettings.colors.list{1});
xlabel('$\sigma_a$','FontSize',sSettings.font.size.axis,'fontname',sSettings.font.name,'interpreter','latex');
set(gca,'XTick',[0,0.005, vSigma_a(end)],'FontSize',sSettings.font.size.axisticks,'fontname',sSettings.font.name);

ylabel('Impact effect on unemp. (ppts.)','FontSize',sSettings.font.size.axis,'fontname',sSettings.font.name,'interpreter','none');
xlim([0,0.01])
ylim([0,1])
box on
title('Uncertainty shock (fixed hiring prob.)','FontSize',sSettings.font.size.default,'fontname',sSettings.font.name,'FontWeight','normal');
ax = gca;
ax.XAxis.Exponent = 0; 
ax.YAxis.Exponent = 0; 

subplot(2,2,4)
p1=plot(vSigma_a,100*mResults_1(uPos,:),'LineWidth',sSettings.lines.width,'linestyle',sSettings.lines.list{1},'Color',sSettings.colors.list{1});
xlabel('$\sigma_a$','FontSize',sSettings.font.size.axis,'fontname',sSettings.font.name,'interpreter','latex');
set(gca,'XTick',[vSigma_a(1):((vSigma_a(end)-vSigma_a(1))/2):vSigma_a(end)],'FontSize',sSettings.font.size.axisticks,'fontname',sSettings.font.name);
ylabel('Impact effect on unemp. (ppts.)','FontSize',sSettings.font.size.axis,'fontname',sSettings.font.name,'interpreter','none');
box on
title('Uncertainty shock (variable hiring prob.)','FontSize',sSettings.font.size.default,'fontname',sSettings.font.name,'FontWeight','normal');
ax = gca;
ax.XAxis.Exponent = 0; 
ax.YAxis.Exponent = 0; 

% Print
xSize = 2*8.75; ySize = 2*6.25;  xCut = 1; yCut = 0.5;
set(gcf,'Units','centimeters','Position',[0 0 xSize ySize],'PaperUnits','centimeters' ...
      ,'PaperPosition',[0 0 xSize ySize],'PaperSize',[xSize-xCut ySize-yCut],'PaperPositionMode','auto')
print(fig,'Output\fig_SaMOptionValue_NoRecalib_Loop00005to001_2x2','-dpdf','-painters')
%}