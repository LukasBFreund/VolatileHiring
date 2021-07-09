%%=========================================================================
% Replication code for den Haan, Freund and Rendahl (2021)
%
% 2-period version of the SaM model with option value
%
% Plot value of waiting with fixed hiring rate
%
% Run on Matlab R2019b
% Last updated: November 2020
% For any questions, please email lukas.beat.freund@gmail.com
%%=========================================================================

%% Housekeeping
%--------------------------------------------------------------------------
clear; close all; clc;
TimeStart = tic;

% Printing?
optionPrint = 1;

% Adjust some style options
load inputSettings
sSettings.plots.xSize = 8.75;
sSettings.plots.ySize = 6.25;
sSettings.plots.xCut = 0; sSettings.plots.yCut = 0;

sSettings.colors.list{5}=[0.7 0.7 0.7];
sSettings.colors.list{5}=[0.9 0.9 0.9];

%% Parameterization
%--------------------------------------------------------------------------
sPar.zss = 1;

% By construction, beta = 1 and delta = 0; we don't target u-rate or 
% vacancies, so set  Upsilon = 1
sPar.Upsilon = 1;
sPar.alpha = 0;                     % exogenous hiring rate!

% As in the full model, we target p=1/2
pTarget = 1/2; 

% To facilitate comparison between stochastic and non-stochastic hiring
% specifications, target hBar=1
hTarget = 1;
sPar.psi = hTarget*(pTarget^(sPar.alpha));

% We target kappa such that cutoff = 0 in model with no risk
% h*(1+0) = kappa => this is just sPar.kappa = hTarget
sPar.kappa = hTarget;

% Dispersion
% Here we take a simple approach and just treat aBar and Delta as fixed parameters
aBar = 0.1;
Delta = 0.05;   


%% Period 1
%--------------------------------------------------------------------------
syms a
assume(-aBar<=a<=aBar)

p1 = (aBar-a)/(2*aBar);                    % prob(a>aHat); where symbolic a here stands for aHat1
h1 = sPar.psi*(p1*sPar.Upsilon).^(-sPar.alpha);

profits1 = h1*(sPar.zss+a)-sPar.kappa;     
profits1_Alt = sPar.zss+a-(sPar.kappa/h1); % non-stochastic hiring 

%% Period 2
%--------------------------------------------------------------------------

% Cutoff in period 2: depends on aggregate state
syms z2 

p2 = (aBar-a)./(2*aBar);                         % treat aHat as a
h2 = sPar.psi.*(p2*sPar.Upsilon).^(-sPar.alpha);
aHat2 = sPar.kappa./h2-z2; % cutoff eqn, as a fn of a

% For good or bad aggregate z2
tempG = subs(aHat2,z2,1+Delta);
tempB = subs(aHat2,z2,1-Delta);
if sPar.alpha>0
aHat2G = double((solve(tempG==a))); 
aHat2B = double((solve(tempB==a))); 
elseif sPar.alpha==0
    aHat2G = double(tempG);
    aHat2B = double(tempB);
end

% Implied entry and hiring probabilities
p2G = double(subs(p2,aHat2G));
h2G = double(subs(h2,aHat2G));
p2B = double(subs(p2,aHat2B)); 
h2B = double(subs(h2,aHat2B));

% The following should be approximately 0
CheckG = h2G*(aHat2G+1+Delta)-sPar.kappa;
CheckB = h2B*(aHat2B+1-Delta)-sPar.kappa;


%% Figure
%--------------------------------------------------------------------------

va = linspace(-aBar,aBar,100);
profits2No_cond = max(0,hTarget*(1+a)-sPar.kappa);
profits2G_cond = max(0,h2G*(1+Delta+a)-sPar.kappa);
profits2B_cond = max(0,h2B*(1-Delta+a)-sPar.kappa);

% To shade expected profits
Area1=[double(subs(profits2G_cond,va));double(subs(profits2No_cond,va))];
Area2=[double(subs(profits2B_cond,va));double(subs(profits2No_cond,va))];
x2 = [va, fliplr(va)];
inBetween1 = [Area1(1,:), fliplr(Area1(2,:))];
inBetween2 = [Area2(1,:), fliplr(Area2(2,:))];

% Figure 
fig1=figure;
PlotArea1 = fill(x2, inBetween1, sSettings.colors.list{3}); % [0.9 0.9 0.9]);
set(PlotArea1,'EdgeAlpha',0);
hold on
PlotArea2 = fill(x2, inBetween2, sSettings.colors.list{4}); % [0.7 0.7 0.7]);
set(PlotArea2,'EdgeAlpha',0);

p1 = plot(va,double(subs(profits2No_cond,va)),sSettings.lines.list{1},'color',sSettings.colors.list{1},'linewidth',sSettings.lines.width);
p2 = plot(va,double(subs(profits2G_cond,va)),sSettings.lines.list{2},'color',sSettings.colors.list{3},'linewidth',sSettings.lines.width);
p3 = plot(va,double(subs(profits2B_cond,va)),sSettings.lines.list{2},'color',sSettings.colors.list{4},'linewidth',sSettings.lines.width);
xlabel('Idiosyncratic productivity level','FontSize',sSettings.font.size.default,'fontname',sSettings.font.name,'interpreter','none');
ylabel('Period-2 profits','FontSize',sSettings.font.size.default,'fontname',sSettings.font.name,'interpreter','none');
set(gca,'XTick',[va(1):(va(end)-va(1))/4:va(end)],'FontSize',sSettings.font.size.axisticks,'fontname',sSettings.font.name);

% Add annotations
txt = '$\hat{a}_{2,+} \searrow$';
text(aHat2G-0.0275,+0.0076,txt,'Interpreter','latex','FontSize',sSettings.font.size.axisticks,'fontname',sSettings.font.name);
txt = '$\hat{a}_{2,-} \downarrow$';
text(aHat2B-0.02,+0.0075,txt,'Interpreter','latex','FontSize',sSettings.font.size.axisticks,'fontname',sSettings.font.name);

txt = sprintf('Increase in profits  \nin expansion');
txt=[txt,'$\rightarrow$'];
text(-0.05,+0.06,txt,'Interpreter','latex','color',sSettings.colors.list{1},'FontSize',6,'fontname',sSettings.font.name);
txt = sprintf(' Decrease \nin profits \nin recession');
txt=['$\leftarrow$',txt];
text(0.063,+0.01,txt,'Interpreter','latex','color',sSettings.colors.list{1},'FontSize',6,'fontname',sSettings.font.name);

legend1 = legend([p1 p2 p3],'$\max\{0,-\kappa+\bar{h}(1+a)\}$',...
    '$\max\{0,-\kappa+\bar{h}(1+\Delta+a)\}$','$\max\{0,-\kappa+\bar{h}(1-\Delta+a)\}$');
set(legend1,'fontname','times','Location','best','FontSize',sSettings.font.size.legend,'interpreter','latex');
set(gcf,'Units','centimeters','Position',[0 0 sSettings.plots.xSize sSettings.plots.ySize],'PaperUnits','centimeters' ...
      ,'PaperPosition',[0 0 sSettings.plots.xSize sSettings.plots.ySize],'PaperSize',[sSettings.plots.xSize-sSettings.plots.xCut sSettings.plots.ySize-sSettings.plots.yCut],'PaperPositionMode','auto')

xlim([-aBar,aBar+0.02])
vline(aBar,':k');
txt=('$\leftarrow \bar{a}$');
text(0.1,0.13,txt,'Interpreter','latex','color',sSettings.colors.list{1},'FontSize',6,'fontname',sSettings.font.name);

if optionPrint == 1
print(fig1,'Output\fig_App_2P_Period2Profits_hFixed_hBar1','-dpdf','-painters')
end



