%%=========================================================================
% Replication code for den Haan, Freund and Rendahl (2021)
%
% 2-period version of the SaM model with option value
%
% Plot cutoff determination with fixed hiring rate
%
% Run on Matlab R2019b
% Last updated: November 2020
% For any questions, please email lukas.beat.freund@gmail.com
%%=========================================================================

%% Housekeeping
%--------------------------------------------------------------------------
clear; close all; clc;
TimeStart = tic;

% Print?
optionPrint = 1;

% Adjust some style options
load inputSettings
sSettings.plots.xSize = 8.75;
sSettings.plots.ySize = 6.25;
sSettings.plots.xCut = 0; sSettings.plots.yCut = 0;

%% Parameterization
%--------------------------------------------------------------------------
sPar.zss = 1;

%% Search frictions
% By construction, beta = 1 and delta = 0; we don't target u-rate or 
% vacancies, so set  Upsilon = 1
sPar.Upsilon = 1;
sPar.alpha = 0;               % exogenous hiring rate!

% As in the full model, we target p=1/2
pTarget = 1/2; 

% For comparability reasons, target h = 1
hTarget = 1;
sPar.psi = hTarget*(pTarget^(sPar.alpha));
sPar.kappa = hTarget;

% Dispersion: now treat only aBar as fixed, Delta as symbolic
aBar = 0.1;
syms Delta
assume(0<Delta<1)

%% Period 1
%--------------------------------------------------------------------------
syms a
assume(-aBar<=a<=aBar)
va = linspace(-aBar,aBar,100);

% Entry probability
p1 = (aBar-a)/(2*aBar);                         % prob(a>aHat); where symbolic a here stands for aHat1
h1 = sPar.psi*(p1*sPar.Upsilon).^(-sPar.alpha);
profits1 =  h1*(sPar.zss+a)-sPar.kappa;

%% Period 2
%--------------------------------------------------------------------------
% Cutoff in period 2: depends on aggregate state
syms z2 
p2 = (aBar-a)/(2*aBar);                         % treat aHat as a
h2 = sPar.psi*(p2*sPar.Upsilon).^(-sPar.alpha);
aHat2 = sPar.kappa/h2-z2;                       

% For good or bad aggregate z2
tempG = subs(aHat2,z2,1+Delta);
tempB = subs(aHat2,z2,1-Delta);
if sPar.alpha>0
    aHat2G = ((solve(tempG==a))); 
    aHat2B = ((solve(tempB==a))); 
elseif sPar.alpha==0
    aHat2G = (tempG);
    aHat2B = (tempB);
end

% It's useful to pin down the following
p2G = (subs(p2,aHat2G));
h2G = (subs(h2,aHat2G));
p2B = (subs(p2,aHat2B)); 
h2B = (subs(h2,aHat2B));

% The following should be approx. 0
CheckG = h2G*(aHat2G+1+Delta)-sPar.kappa;
CheckB = h2B*(aHat2B+1-Delta)-sPar.kappa;

% Profits for individual firm, *conditional* on entry 
profits2G = h2G*(1+Delta+a)-sPar.kappa;
profits2B = h2B*(1-Delta+a)-sPar.kappa;
%profits2G_Alt =  (1+Delta+a)-sPar.kappa/h2G;
%profits2B_Alt =  (1-Delta+a)-sPar.kappa/h2B;

% To compute expected profits, need aStar (cond. on entry)
aStar2G = (aHat2G+aBar)/2;
aStar2B = (aHat2B+aBar)/2;

EProfits2 = 0.5*(p2G*subs(profits2G,a,aStar2G))...
            + 0.5*(p2B*subs(profits2B,a,aStar2B)); 


%% Cutoff determination in period 1
%--------------------------------------------------------------------------

% Solve h*J(aHat)-kappa = E[J2*] for particular values of Delta 
DeltaValue1 = 0;
DeltaValue2 = 0.05;
vDelta = linspace(DeltaValue1,DeltaValue2,20);

for iD = 1:numel(vDelta)
    DeltaValue = vDelta(iD);
    vaHat1_Sol(iD) = (vpasolve(profits1-subs(EProfits2,DeltaValue)==0)); 
end

% Compute implied hiring rate
%h1_Sol_1 = (subs(h1,vaHat1_Sol(1)));
%h1_Sol_2 = (subs(h1,vaHat1_Sol(end)));


%% Figure
%--------------------------------------------------------------------------
if sPar.alpha==0
    xLow = -aBar;
    xHigh = +aBar;
else
    xLow = double(vaHat1_Sol(1) - 0.01);
    xHigh = double(vaHat1_Sol(end) + 0.01);
end


% For period 2, compute expected profits and upper and lower case
Profits2G = (1+Delta-sPar.kappa+(aBar-Delta)/2);
Profits2B = (1-Delta-sPar.kappa+(aBar+Delta)/2);
mE2=[ones(1,size(va,2))*double(subs(Profits2G,DeltaValue2));ones(1,size(va,2))*double(subs(Profits2B,DeltaValue2))];

x2 = [va, fliplr(va)];
inBetween = [mE2(1,:), fliplr(mE2(2,:))];

fig1=figure;
p1 = fplot(max(profits1,0),[xLow xHigh],sSettings.lines.list{1},'color',sSettings.colors.list{1},'linewidth',sSettings.lines.width);
hold on
p2 = plot(va,ones(1,size(va,2))*double(subs(EProfits2,DeltaValue1)),sSettings.lines.list{1},'color',sSettings.colors.list{2},'linewidth',1);
p3 = plot(va,ones(1,size(va,2))*double(subs(EProfits2,DeltaValue2)),sSettings.lines.list{2},'color',sSettings.colors.list{3},'linewidth',sSettings.lines.width);
set(gca,'XTick',[xLow:(xHigh-xLow)/4:xHigh],'FontSize',sSettings.font.size.axisticks,'fontname',sSettings.font.name);
xlabel('Idiosyncratic productivity level, a','FontSize',sSettings.font.size.axis,'fontname',sSettings.font.name,'interpreter','none');
ylabel('Profits','FontSize',sSettings.font.size.axis,'fontname',sSettings.font.name,'interpreter','none');
xlim([xLow xHigh])
ylim([0 0.1])

% Annotation
vline(vaHat1_Sol(1),'k:');
vline(vaHat1_Sol(end),'k:');

txt = '$\hat{a}_1 (\Delta = 0) \searrow$';
text(vaHat1_Sol(1)-0.06,0.01,txt,'Interpreter','latex','FontSize',sSettings.font.size.axisticks,'fontname',sSettings.font.name);

txt = '$\swarrow \hat{a}_1 (\Delta > 0)$';
text(vaHat1_Sol(end)+0.001,+0.005,txt,'Interpreter','latex','FontSize',sSettings.font.size.axisticks,'fontname',sSettings.font.name);


legend1 = legend([p1 p2 p3],'Period-1 profits (at cutoff)',...
    'Exp. Period-2 profits, $\Delta = 0$','Exp. Period-2 profits, $\Delta > 0$','interpreter','latex');
set(legend1,'fontname','times','Location','northwest','FontSize',sSettings.font.size.legend);

set(gcf,'Units','centimeters','Position',[0 0 sSettings.plots.xSize sSettings.plots.ySize],'PaperUnits','centimeters' ...
      ,'PaperPosition',[0 0 sSettings.plots.xSize sSettings.plots.ySize],'PaperSize',...
          [sSettings.plots.xSize-sSettings.plots.xCut sSettings.plots.ySize-sSettings.plots.yCut],'PaperPositionMode','auto')

 xlim([-aBar,aBar+0.02])
vline(aBar,':k');
txt=('$\leftarrow \bar{a}$');
text(0.1,0.08,txt,'Interpreter','latex','color',sSettings.colors.list{1},'FontSize',6,'fontname',sSettings.font.name);
  ax = gca;
        ax.XAxis.Exponent = 0; 
     
if optionPrint == 1
print(fig1,'Output\fig_App_2P_CutoffDetermination_hFixed_hBar1_aBar01_Delta005','-dpdf','-painters')
end

