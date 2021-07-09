function [residual, g1, g2, g3] = dynareSAMOptionValue_dynamic(y, x, params, steady_state, it_)
%
% Status : Computes dynamic model for Dynare
%
% Inputs :
%   y         [#dynamic variables by 1] double    vector of endogenous variables in the order stored
%                                                 in M_.lead_lag_incidence; see the Manual
%   x         [M_.exo_nbr by nperiods] double     matrix of exogenous variables (in declaration order)
%                                                 for all simulation periods
%   params    [M_.param_nbr by 1] double          vector of parameter values in declaration order
%   it_       scalar double                       time period for exogenous variables for which to evaluate the model
%
% Outputs:
%   residual  [M_.endo_nbr by 1] double    vector of residuals of the dynamic model equations in order of 
%                                          declaration of the equations
%   g1        [M_.endo_nbr by #dynamic variables] double    Jacobian matrix of the dynamic model equations;
%                                                           rows: equations in order of declaration
%                                                           columns: variables in order stored in M_.lead_lag_incidence
%   g2        [M_.endo_nbr by (#dynamic variables)^2] double   Hessian matrix of the dynamic model equations;
%                                                              rows: equations in order of declaration
%                                                              columns: variables in order stored in M_.lead_lag_incidence
%   g3        [M_.endo_nbr by (#dynamic variables)^3] double   Third order derivative matrix of the dynamic model equations;
%                                                              rows: equations in order of declaration
%                                                              columns: variables in order stored in M_.lead_lag_incidence
%
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

%
% Model equations
%

residual = zeros(25, 1);
T19 = params(8)*y(7)^params(3);
T22 = y(8)^(1-params(3));
T92 = (y(15)-params(17))/params(15);
T100 = 1/sqrt(2*params(18));
T105 = exp((-0.5)*T92^2);
T174 = params(8)*getPowerDeriv(y(7),params(3),1);
T181 = getPowerDeriv(y(8),1-params(3),1);
T207 = 1/params(15);
T213 = exp((-(T92*T92))/2);
T220 = T105*(-0.5)*T207*2*T92;
T269 = T213*(-(T92*T207+T92*T207))/2;
lhs =y(7);
rhs =1-(1-params(4))*y(1);
residual(1)= lhs-rhs;
lhs =y(4);
rhs =T19*T22;
residual(2)= lhs-rhs;
lhs =y(11);
rhs =y(4)/y(7);
residual(3)= lhs-rhs;
lhs =y(9);
rhs =y(4)/y(8);
residual(4)= lhs-rhs;
lhs =y(5);
rhs =(1-params(4))*y(1)+y(4);
residual(5)= lhs-rhs;
lhs =y(6);
rhs =1-y(5);
residual(6)= lhs-rhs;
lhs =y(17);
rhs =(1-params(5))*params(9)*y(15)/(1-(1-params(4))*params(1))+y(22);
residual(7)= lhs-rhs;
lhs =y(22);
rhs =(1-params(5))*(params(9)*y(10)-params(7))+params(1)*(1-params(4))*y(30);
residual(8)= lhs-rhs;
lhs =y(18);
rhs =(1-params(5))*params(9)*y(16)/(1-(1-params(4))*params(1))+y(23);
residual(9)= lhs-rhs;
lhs =y(23);
rhs =(1-params(5))*(params(9)*y(10)-params(7))+params(1)*(1-params(4))*y(31);
residual(10)= lhs-rhs;
lhs =y(19);
rhs =params(1)*y(29)-y(14)*params(6)+y(9)*y(14)*(y(18)-params(1)*y(29));
residual(11)= lhs-rhs;
lhs =params(6);
rhs =y(9)*(y(17)-params(1)*y(29));
residual(12)= lhs-rhs;
lhs =y(14);
rhs =1-normcdf(T92,0,1);
residual(13)= lhs-rhs;
lhs =y(21);
rhs =T100*T105;
residual(14)= lhs-rhs;
lhs =y(16);
rhs =params(17)+params(15)*y(21)/y(14);
residual(15)= lhs-rhs;
lhs =y(8);
rhs =y(14)*(params(16)-(1-params(4))*y(1));
residual(16)= lhs-rhs;
lhs =y(10);
rhs =(1-params(11))*(steady_state(7))+params(11)*y(2)+y(3)*x(it_, 1);
residual(17)= lhs-rhs;
lhs =y(13);
rhs =(1-params(12))*params(13)+y(3)*params(12)+params(14)*x(it_, 2);
residual(18)= lhs-rhs;
lhs =y(12);
rhs =y(5)*(y(10)+y(16));
residual(19)= lhs-rhs;
lhs =y(20);
rhs =(1-params(5))*params(9)*(steady_state(12))/(1-(1-params(4))*params(1))+y(25);
residual(20)= lhs-rhs;
lhs =y(25);
rhs =(1-params(5))*(params(9)*y(10)-params(7))+params(1)*(1-params(4))*y(32);
residual(21)= lhs-rhs;
lhs =y(24);
rhs =params(1)*y(29);
residual(22)= lhs-rhs;
lhs =y(26);
rhs =y(5);
residual(23)= lhs-rhs;
lhs =y(27);
rhs =y(10);
residual(24)= lhs-rhs;
lhs =y(28);
rhs =y(13);
residual(25)= lhs-rhs;
if nargout >= 2,
  g1 = zeros(25, 34);

  %
  % Jacobian matrix
  %

  g1(1,1)=1-params(4);
  g1(1,7)=1;
  g1(2,4)=1;
  g1(2,7)=(-(T22*T174));
  g1(2,8)=(-(T19*T181));
  g1(3,4)=(-(1/y(7)));
  g1(3,7)=(-((-y(4))/(y(7)*y(7))));
  g1(3,11)=1;
  g1(4,4)=(-(1/y(8)));
  g1(4,8)=(-((-y(4))/(y(8)*y(8))));
  g1(4,9)=1;
  g1(5,4)=(-1);
  g1(5,1)=(-(1-params(4)));
  g1(5,5)=1;
  g1(6,5)=1;
  g1(6,6)=1;
  g1(7,15)=(-((1-params(5))*params(9)/(1-(1-params(4))*params(1))));
  g1(7,17)=1;
  g1(7,22)=(-1);
  g1(8,10)=(-((1-params(5))*params(9)));
  g1(8,22)=1;
  g1(8,30)=(-((1-params(4))*params(1)));
  g1(9,16)=(-((1-params(5))*params(9)/(1-(1-params(4))*params(1))));
  g1(9,18)=1;
  g1(9,23)=(-1);
  g1(10,10)=(-((1-params(5))*params(9)));
  g1(10,23)=1;
  g1(10,31)=(-((1-params(4))*params(1)));
  g1(11,9)=(-(y(14)*(y(18)-params(1)*y(29))));
  g1(11,14)=(-((-params(6))+y(9)*(y(18)-params(1)*y(29))));
  g1(11,18)=(-(y(9)*y(14)));
  g1(11,19)=1;
  g1(11,29)=(-(params(1)+y(9)*y(14)*(-params(1))));
  g1(12,9)=(-(y(17)-params(1)*y(29)));
  g1(12,17)=(-y(9));
  g1(12,29)=(-(y(9)*(-params(1))));
  g1(13,14)=1;
  g1(13,15)=T207*T213/2.506628274631;
  g1(14,15)=(-(T100*T220));
  g1(14,21)=1;
  g1(15,14)=(-((-(params(15)*y(21)))/(y(14)*y(14))));
  g1(15,16)=1;
  g1(15,21)=(-(params(15)/y(14)));
  g1(16,1)=(-(y(14)*(-(1-params(4)))));
  g1(16,8)=1;
  g1(16,14)=(-(params(16)-(1-params(4))*y(1)));
  g1(17,2)=(-params(11));
  g1(17,10)=1;
  g1(17,3)=(-x(it_, 1));
  g1(17,33)=(-y(3));
  g1(18,3)=(-params(12));
  g1(18,13)=1;
  g1(18,34)=(-params(14));
  g1(19,5)=(-(y(10)+y(16)));
  g1(19,10)=(-y(5));
  g1(19,12)=1;
  g1(19,16)=(-y(5));
  g1(20,20)=1;
  g1(20,25)=(-1);
  g1(21,10)=(-((1-params(5))*params(9)));
  g1(21,25)=1;
  g1(21,32)=(-((1-params(4))*params(1)));
  g1(22,29)=(-params(1));
  g1(22,24)=1;
  g1(23,5)=(-1);
  g1(23,26)=1;
  g1(24,10)=(-1);
  g1(24,27)=1;
  g1(25,13)=(-1);
  g1(25,28)=1;
end
if nargout >= 3,
  %
  % Hessian matrix
  %

  v2 = zeros(37,3);
  v2(1,1)=2;
  v2(1,2)=211;
  v2(1,3)=(-(T22*params(8)*getPowerDeriv(y(7),params(3),2)));
  v2(2,1)=2;
  v2(2,2)=245;
  v2(2,3)=(-(T174*T181));
  v2(3,1)=2;
  v2(3,2)=212;
  v2(3,3)=  v2(2,3);
  v2(4,1)=2;
  v2(4,2)=246;
  v2(4,3)=(-(T19*getPowerDeriv(y(8),1-params(3),2)));
  v2(5,1)=3;
  v2(5,2)=208;
  v2(5,3)=(-((-1)/(y(7)*y(7))));
  v2(6,1)=3;
  v2(6,2)=109;
  v2(6,3)=  v2(5,3);
  v2(7,1)=3;
  v2(7,2)=211;
  v2(7,3)=(-((-((-y(4))*(y(7)+y(7))))/(y(7)*y(7)*y(7)*y(7))));
  v2(8,1)=4;
  v2(8,2)=242;
  v2(8,3)=(-((-1)/(y(8)*y(8))));
  v2(9,1)=4;
  v2(9,2)=110;
  v2(9,3)=  v2(8,3);
  v2(10,1)=4;
  v2(10,2)=246;
  v2(10,3)=(-((-((-y(4))*(y(8)+y(8))))/(y(8)*y(8)*y(8)*y(8))));
  v2(11,1)=11;
  v2(11,2)=451;
  v2(11,3)=(-(y(18)-params(1)*y(29)));
  v2(12,1)=11;
  v2(12,2)=286;
  v2(12,3)=  v2(11,3);
  v2(13,1)=11;
  v2(13,2)=587;
  v2(13,3)=(-y(14));
  v2(14,1)=11;
  v2(14,2)=290;
  v2(14,3)=  v2(13,3);
  v2(15,1)=11;
  v2(15,2)=592;
  v2(15,3)=(-y(9));
  v2(16,1)=11;
  v2(16,2)=460;
  v2(16,3)=  v2(15,3);
  v2(17,1)=11;
  v2(17,2)=961;
  v2(17,3)=(-(y(14)*(-params(1))));
  v2(18,1)=11;
  v2(18,2)=301;
  v2(18,3)=  v2(17,3);
  v2(19,1)=11;
  v2(19,2)=966;
  v2(19,3)=(-(y(9)*(-params(1))));
  v2(20,1)=11;
  v2(20,2)=471;
  v2(20,3)=  v2(19,3);
  v2(21,1)=12;
  v2(21,2)=553;
  v2(21,3)=(-1);
  v2(22,1)=12;
  v2(22,2)=289;
  v2(22,3)=  v2(21,3);
  v2(23,1)=12;
  v2(23,2)=961;
  v2(23,3)=params(1);
  v2(24,1)=12;
  v2(24,2)=301;
  v2(24,3)=  v2(23,3);
  v2(25,1)=13;
  v2(25,2)=491;
  v2(25,3)=T207*T269/2.506628274631;
  v2(26,1)=14;
  v2(26,2)=491;
  v2(26,3)=(-(T100*((-0.5)*T207*2*T92*T220+T105*(-0.5)*T207*2*T207)));
  v2(27,1)=15;
  v2(27,2)=456;
  v2(27,3)=(-((-((-(params(15)*y(21)))*(y(14)+y(14))))/(y(14)*y(14)*y(14)*y(14))));
  v2(28,1)=15;
  v2(28,2)=694;
  v2(28,3)=(-((-params(15))/(y(14)*y(14))));
  v2(29,1)=15;
  v2(29,2)=463;
  v2(29,3)=  v2(28,3);
  v2(30,1)=16;
  v2(30,2)=443;
  v2(30,3)=1-params(4);
  v2(31,1)=16;
  v2(31,2)=14;
  v2(31,3)=  v2(30,3);
  v2(32,1)=17;
  v2(32,2)=1091;
  v2(32,3)=(-1);
  v2(33,1)=17;
  v2(33,2)=101;
  v2(33,3)=  v2(32,3);
  v2(34,1)=19;
  v2(34,2)=311;
  v2(34,3)=(-1);
  v2(35,1)=19;
  v2(35,2)=146;
  v2(35,3)=  v2(34,3);
  v2(36,1)=19;
  v2(36,2)=515;
  v2(36,3)=(-1);
  v2(37,1)=19;
  v2(37,2)=152;
  v2(37,3)=  v2(36,3);
  g2 = sparse(v2(:,1),v2(:,2),v2(:,3),25,1156);
end
if nargout >= 4,
  %
  % Third order derivatives
  %

  v3 = zeros(34,3);
  v3(1,1)=2;
  v3(1,2)=7147;
  v3(1,3)=(-(T22*params(8)*getPowerDeriv(y(7),params(3),3)));
  v3(2,1)=2;
  v3(2,2)=8303;
  v3(2,3)=(-(T181*params(8)*getPowerDeriv(y(7),params(3),2)));
  v3(3,1)=2;
  v3(3,2)=7148;
  v3(3,3)=  v3(2,3);
  v3(4,1)=2;
  v3(4,2)=7181;
  v3(4,3)=  v3(2,3);
  v3(5,1)=2;
  v3(5,2)=8337;
  v3(5,3)=(-(T174*getPowerDeriv(y(8),1-params(3),2)));
  v3(6,1)=2;
  v3(6,2)=7182;
  v3(6,3)=  v3(5,3);
  v3(7,1)=2;
  v3(7,2)=8304;
  v3(7,3)=  v3(5,3);
  v3(8,1)=2;
  v3(8,2)=8338;
  v3(8,3)=(-(T19*getPowerDeriv(y(8),1-params(3),3)));
  v3(9,1)=3;
  v3(9,2)=7144;
  v3(9,3)=(-((y(7)+y(7))/(y(7)*y(7)*y(7)*y(7))));
  v3(10,1)=3;
  v3(10,2)=3679;
  v3(10,3)=  v3(9,3);
  v3(11,1)=3;
  v3(11,2)=7045;
  v3(11,3)=  v3(9,3);
  v3(12,1)=3;
  v3(12,2)=7147;
  v3(12,3)=(-((y(7)*y(7)*y(7)*y(7)*(-(2*(-y(4))))-(-((-y(4))*(y(7)+y(7))))*(y(7)*y(7)*(y(7)+y(7))+y(7)*y(7)*(y(7)+y(7))))/(y(7)*y(7)*y(7)*y(7)*y(7)*y(7)*y(7)*y(7))));
  v3(13,1)=4;
  v3(13,2)=8334;
  v3(13,3)=(-((y(8)+y(8))/(y(8)*y(8)*y(8)*y(8))));
  v3(14,1)=4;
  v3(14,2)=3714;
  v3(14,3)=  v3(13,3);
  v3(15,1)=4;
  v3(15,2)=8202;
  v3(15,3)=  v3(13,3);
  v3(16,1)=4;
  v3(16,2)=8338;
  v3(16,3)=(-((y(8)*y(8)*y(8)*y(8)*(-(2*(-y(4))))-(-((-y(4))*(y(8)+y(8))))*(y(8)*y(8)*(y(8)+y(8))+y(8)*y(8)*(y(8)+y(8))))/(y(8)*y(8)*y(8)*y(8)*y(8)*y(8)*y(8)*y(8))));
  v3(17,1)=11;
  v3(17,2)=20103;
  v3(17,3)=(-1);
  v3(18,1)=11;
  v3(18,2)=9708;
  v3(18,3)=  v3(17,3);
  v3(19,1)=11;
  v3(19,2)=9840;
  v3(19,3)=  v3(17,3);
  v3(20,1)=11;
  v3(20,2)=15318;
  v3(20,3)=  v3(17,3);
  v3(21,1)=11;
  v3(21,2)=15615;
  v3(21,3)=  v3(17,3);
  v3(22,1)=11;
  v3(22,2)=19938;
  v3(22,3)=  v3(17,3);
  v3(23,1)=11;
  v3(23,2)=32819;
  v3(23,3)=params(1);
  v3(24,1)=11;
  v3(24,2)=9719;
  v3(24,3)=  v3(23,3);
  v3(25,1)=11;
  v3(25,2)=10214;
  v3(25,3)=  v3(23,3);
  v3(26,1)=11;
  v3(26,2)=15329;
  v3(26,3)=  v3(23,3);
  v3(27,1)=11;
  v3(27,2)=15989;
  v3(27,3)=  v3(23,3);
  v3(28,1)=11;
  v3(28,2)=32654;
  v3(28,3)=  v3(23,3);
  v3(29,1)=13;
  v3(29,2)=16675;
  v3(29,3)=T207*((-(T92*T207+T92*T207))/2*T269+T213*(-(T207*T207+T207*T207))/2)/2.506628274631;
  v3(30,1)=14;
  v3(30,2)=16675;
  v3(30,3)=(-(T100*(T220*(-0.5)*T207*2*T207+T220*(-0.5)*T207*2*T207+(-0.5)*T207*2*T92*((-0.5)*T207*2*T92*T220+T105*(-0.5)*T207*2*T207))));
  v3(31,1)=15;
  v3(31,2)=15484;
  v3(31,3)=(-((y(14)*y(14)*y(14)*y(14)*(-(2*(-(params(15)*y(21)))))-(-((-(params(15)*y(21)))*(y(14)+y(14))))*(y(14)*y(14)*(y(14)+y(14))+y(14)*y(14)*(y(14)+y(14))))/(y(14)*y(14)*y(14)*y(14)*y(14)*y(14)*y(14)*y(14))));
  v3(32,1)=15;
  v3(32,2)=23576;
  v3(32,3)=(-((-((y(14)+y(14))*(-params(15))))/(y(14)*y(14)*y(14)*y(14))));
  v3(33,1)=15;
  v3(33,2)=15491;
  v3(33,3)=  v3(32,3);
  v3(34,1)=15;
  v3(34,2)=15722;
  v3(34,3)=  v3(32,3);
  g3 = sparse(v3(:,1),v3(:,2),v3(:,3),25,39304);
end
end
