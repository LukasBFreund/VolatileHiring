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

residual = zeros(27, 1);
T19 = params(8)*y(8)^params(3);
T22 = y(9)^(1-params(3));
T183 = params(8)*getPowerDeriv(y(8),params(3),1);
T190 = getPowerDeriv(y(9),1-params(3),1);
lhs =y(8);
rhs =1-(1-params(4))*y(1);
residual(1)= lhs-rhs;
lhs =y(5);
rhs =T19*T22;
residual(2)= lhs-rhs;
lhs =y(12);
rhs =y(5)/y(8);
residual(3)= lhs-rhs;
lhs =y(10);
rhs =y(5)/y(9);
residual(4)= lhs-rhs;
lhs =y(6);
rhs =(1-params(4))*y(1)+y(5);
residual(5)= lhs-rhs;
lhs =y(7);
rhs =1-y(6);
residual(6)= lhs-rhs;
lhs =y(18);
rhs =(1-params(5))*params(9)*y(16)/(1-(1-params(4))*params(1))+y(22);
residual(7)= lhs-rhs;
lhs =y(22);
rhs =(1-params(5))*(params(9)*y(11)-params(7))+params(1)*((1-params(4))*y(33)+params(4)*y(36));
residual(8)= lhs-rhs;
lhs =y(19);
rhs =(1-params(5))*params(9)*y(17)/(1-(1-params(4))*params(1))+y(23);
residual(9)= lhs-rhs;
lhs =y(23);
rhs =(1-params(5))*(params(9)*y(11)-params(7))+params(1)*(params(4)*y(36)+(1-params(4))*y(34));
residual(10)= lhs-rhs;
lhs =y(20);
rhs =params(1)*y(32)-y(15)*params(6)+y(10)*y(15)*(y(19)-params(1)*y(32));
residual(11)= lhs-rhs;
lhs =params(6);
rhs =y(10)*(y(18)-params(1)*y(32));
residual(12)= lhs-rhs;
lhs =y(15);
rhs =1-(y(16)-params(19))/(params(18)-params(19));
residual(13)= lhs-rhs;
lhs =y(17);
rhs =0.5*(y(16)+params(18));
residual(14)= lhs-rhs;
lhs =y(9);
rhs =y(15)*y(29);
residual(15)= lhs-rhs;
lhs =y(29);
rhs =params(16)-(1-params(4))*y(1)-params(4)*params(20)*y(4);
residual(16)= lhs-rhs;
lhs =y(30);
rhs =y(6)+params(20)*y(4);
residual(17)= lhs-rhs;
lhs =y(31);
rhs =y(20)*(1-params(20))+y(36)*params(1)*params(20);
residual(18)= lhs-rhs;
lhs =y(11);
rhs =(1-params(11))*(steady_state(7))+params(11)*y(2)+y(3)*x(it_, 1);
residual(19)= lhs-rhs;
lhs =y(14);
rhs =(1-params(12))*params(13)+y(3)*params(12)+params(14)*x(it_, 2);
residual(20)= lhs-rhs;
lhs =y(13);
rhs =y(6)*(y(11)+y(17));
residual(21)= lhs-rhs;
lhs =y(21);
rhs =(1-params(5))*params(9)*(steady_state(12))/(1-(1-params(4))*params(1))+y(25);
residual(22)= lhs-rhs;
lhs =y(25);
rhs =(1-params(5))*(params(9)*y(11)-params(7))+params(1)*(params(4)*y(36)+(1-params(4))*y(35));
residual(23)= lhs-rhs;
lhs =y(24);
rhs =params(1)*y(32);
residual(24)= lhs-rhs;
lhs =y(26);
rhs =y(6);
residual(25)= lhs-rhs;
lhs =y(27);
rhs =y(11);
residual(26)= lhs-rhs;
lhs =y(28);
rhs =y(14);
residual(27)= lhs-rhs;
if nargout >= 2,
  g1 = zeros(27, 38);

  %
  % Jacobian matrix
  %

  g1(1,1)=1-params(4);
  g1(1,8)=1;
  g1(2,5)=1;
  g1(2,8)=(-(T22*T183));
  g1(2,9)=(-(T19*T190));
  g1(3,5)=(-(1/y(8)));
  g1(3,8)=(-((-y(5))/(y(8)*y(8))));
  g1(3,12)=1;
  g1(4,5)=(-(1/y(9)));
  g1(4,9)=(-((-y(5))/(y(9)*y(9))));
  g1(4,10)=1;
  g1(5,5)=(-1);
  g1(5,1)=(-(1-params(4)));
  g1(5,6)=1;
  g1(6,6)=1;
  g1(6,7)=1;
  g1(7,16)=(-((1-params(5))*params(9)/(1-(1-params(4))*params(1))));
  g1(7,18)=1;
  g1(7,22)=(-1);
  g1(8,11)=(-((1-params(5))*params(9)));
  g1(8,22)=1;
  g1(8,33)=(-((1-params(4))*params(1)));
  g1(8,36)=(-(params(4)*params(1)));
  g1(9,17)=(-((1-params(5))*params(9)/(1-(1-params(4))*params(1))));
  g1(9,19)=1;
  g1(9,23)=(-1);
  g1(10,11)=(-((1-params(5))*params(9)));
  g1(10,23)=1;
  g1(10,34)=(-((1-params(4))*params(1)));
  g1(10,36)=(-(params(4)*params(1)));
  g1(11,10)=(-(y(15)*(y(19)-params(1)*y(32))));
  g1(11,15)=(-((-params(6))+y(10)*(y(19)-params(1)*y(32))));
  g1(11,19)=(-(y(10)*y(15)));
  g1(11,20)=1;
  g1(11,32)=(-(params(1)+y(10)*y(15)*(-params(1))));
  g1(12,10)=(-(y(18)-params(1)*y(32)));
  g1(12,18)=(-y(10));
  g1(12,32)=(-(y(10)*(-params(1))));
  g1(13,15)=1;
  g1(13,16)=1/(params(18)-params(19));
  g1(14,16)=(-0.5);
  g1(14,17)=1;
  g1(15,9)=1;
  g1(15,15)=(-y(29));
  g1(15,29)=(-y(15));
  g1(16,1)=1-params(4);
  g1(16,29)=1;
  g1(16,4)=params(4)*params(20);
  g1(17,6)=(-1);
  g1(17,4)=(-params(20));
  g1(17,30)=1;
  g1(18,20)=(-(1-params(20)));
  g1(18,31)=1;
  g1(18,36)=(-(params(1)*params(20)));
  g1(19,2)=(-params(11));
  g1(19,11)=1;
  g1(19,3)=(-x(it_, 1));
  g1(19,37)=(-y(3));
  g1(20,3)=(-params(12));
  g1(20,14)=1;
  g1(20,38)=(-params(14));
  g1(21,6)=(-(y(11)+y(17)));
  g1(21,11)=(-y(6));
  g1(21,13)=1;
  g1(21,17)=(-y(6));
  g1(22,21)=1;
  g1(22,25)=(-1);
  g1(23,11)=(-((1-params(5))*params(9)));
  g1(23,25)=1;
  g1(23,35)=(-((1-params(4))*params(1)));
  g1(23,36)=(-(params(4)*params(1)));
  g1(24,32)=(-params(1));
  g1(24,24)=1;
  g1(25,6)=(-1);
  g1(25,26)=1;
  g1(26,11)=(-1);
  g1(26,27)=1;
  g1(27,14)=(-1);
  g1(27,28)=1;
end
if nargout >= 3,
  %
  % Hessian matrix
  %

  v2 = zeros(32,3);
  v2(1,1)=2;
  v2(1,2)=274;
  v2(1,3)=(-(T22*params(8)*getPowerDeriv(y(8),params(3),2)));
  v2(2,1)=2;
  v2(2,2)=312;
  v2(2,3)=(-(T183*T190));
  v2(3,1)=2;
  v2(3,2)=275;
  v2(3,3)=  v2(2,3);
  v2(4,1)=2;
  v2(4,2)=313;
  v2(4,3)=(-(T19*getPowerDeriv(y(9),1-params(3),2)));
  v2(5,1)=3;
  v2(5,2)=271;
  v2(5,3)=(-((-1)/(y(8)*y(8))));
  v2(6,1)=3;
  v2(6,2)=160;
  v2(6,3)=  v2(5,3);
  v2(7,1)=3;
  v2(7,2)=274;
  v2(7,3)=(-((-((-y(5))*(y(8)+y(8))))/(y(8)*y(8)*y(8)*y(8))));
  v2(8,1)=4;
  v2(8,2)=309;
  v2(8,3)=(-((-1)/(y(9)*y(9))));
  v2(9,1)=4;
  v2(9,2)=161;
  v2(9,3)=  v2(8,3);
  v2(10,1)=4;
  v2(10,2)=313;
  v2(10,3)=(-((-((-y(5))*(y(9)+y(9))))/(y(9)*y(9)*y(9)*y(9))));
  v2(11,1)=11;
  v2(11,2)=542;
  v2(11,3)=(-(y(19)-params(1)*y(32)));
  v2(12,1)=11;
  v2(12,2)=357;
  v2(12,3)=  v2(11,3);
  v2(13,1)=11;
  v2(13,2)=694;
  v2(13,3)=(-y(15));
  v2(14,1)=11;
  v2(14,2)=361;
  v2(14,3)=  v2(13,3);
  v2(15,1)=11;
  v2(15,2)=699;
  v2(15,3)=(-y(10));
  v2(16,1)=11;
  v2(16,2)=551;
  v2(16,3)=  v2(15,3);
  v2(17,1)=11;
  v2(17,2)=1188;
  v2(17,3)=(-(y(15)*(-params(1))));
  v2(18,1)=11;
  v2(18,2)=374;
  v2(18,3)=  v2(17,3);
  v2(19,1)=11;
  v2(19,2)=1193;
  v2(19,3)=(-(y(10)*(-params(1))));
  v2(20,1)=11;
  v2(20,2)=564;
  v2(20,3)=  v2(19,3);
  v2(21,1)=12;
  v2(21,2)=656;
  v2(21,3)=(-1);
  v2(22,1)=12;
  v2(22,2)=360;
  v2(22,3)=  v2(21,3);
  v2(23,1)=12;
  v2(23,2)=1188;
  v2(23,3)=params(1);
  v2(24,1)=12;
  v2(24,2)=374;
  v2(24,3)=  v2(23,3);
  v2(25,1)=15;
  v2(25,2)=1079;
  v2(25,3)=(-1);
  v2(26,1)=15;
  v2(26,2)=561;
  v2(26,3)=  v2(25,3);
  v2(27,1)=19;
  v2(27,2)=1371;
  v2(27,3)=(-1);
  v2(28,1)=19;
  v2(28,2)=113;
  v2(28,3)=  v2(27,3);
  v2(29,1)=21;
  v2(29,2)=386;
  v2(29,3)=(-1);
  v2(30,1)=21;
  v2(30,2)=201;
  v2(30,3)=  v2(29,3);
  v2(31,1)=21;
  v2(31,2)=614;
  v2(31,3)=(-1);
  v2(32,1)=21;
  v2(32,2)=207;
  v2(32,3)=  v2(31,3);
  g2 = sparse(v2(:,1),v2(:,2),v2(:,3),27,1444);
end
if nargout >= 4,
  %
  % Third order derivatives
  %

  v3 = zeros(28,3);
  v3(1,1)=2;
  v3(1,2)=10382;
  v3(1,3)=(-(T22*params(8)*getPowerDeriv(y(8),params(3),3)));
  v3(2,1)=2;
  v3(2,2)=11826;
  v3(2,3)=(-(T190*params(8)*getPowerDeriv(y(8),params(3),2)));
  v3(3,1)=2;
  v3(3,2)=10383;
  v3(3,3)=  v3(2,3);
  v3(4,1)=2;
  v3(4,2)=10420;
  v3(4,3)=  v3(2,3);
  v3(5,1)=2;
  v3(5,2)=11864;
  v3(5,3)=(-(T183*getPowerDeriv(y(9),1-params(3),2)));
  v3(6,1)=2;
  v3(6,2)=10421;
  v3(6,3)=  v3(5,3);
  v3(7,1)=2;
  v3(7,2)=11827;
  v3(7,3)=  v3(5,3);
  v3(8,1)=2;
  v3(8,2)=11865;
  v3(8,3)=(-(T19*getPowerDeriv(y(9),1-params(3),3)));
  v3(9,1)=3;
  v3(9,2)=10379;
  v3(9,3)=(-((y(8)+y(8))/(y(8)*y(8)*y(8)*y(8))));
  v3(10,1)=3;
  v3(10,2)=6050;
  v3(10,3)=  v3(9,3);
  v3(11,1)=3;
  v3(11,2)=10268;
  v3(11,3)=  v3(9,3);
  v3(12,1)=3;
  v3(12,2)=10382;
  v3(12,3)=(-((y(8)*y(8)*y(8)*y(8)*(-(2*(-y(5))))-(-((-y(5))*(y(8)+y(8))))*(y(8)*y(8)*(y(8)+y(8))+y(8)*y(8)*(y(8)+y(8))))/(y(8)*y(8)*y(8)*y(8)*y(8)*y(8)*y(8)*y(8))));
  v3(13,1)=4;
  v3(13,2)=11861;
  v3(13,3)=(-((y(9)+y(9))/(y(9)*y(9)*y(9)*y(9))));
  v3(14,1)=4;
  v3(14,2)=6089;
  v3(14,3)=  v3(13,3);
  v3(15,1)=4;
  v3(15,2)=11713;
  v3(15,3)=  v3(13,3);
  v3(16,1)=4;
  v3(16,2)=11865;
  v3(16,3)=(-((y(9)*y(9)*y(9)*y(9)*(-(2*(-y(5))))-(-((-y(5))*(y(9)+y(9))))*(y(9)*y(9)*(y(9)+y(9))+y(9)*y(9)*(y(9)+y(9))))/(y(9)*y(9)*y(9)*y(9)*y(9)*y(9)*y(9)*y(9))));
  v3(17,1)=11;
  v3(17,2)=26534;
  v3(17,3)=(-1);
  v3(18,1)=11;
  v3(18,2)=13547;
  v3(18,3)=  v3(17,3);
  v3(19,1)=11;
  v3(19,2)=13695;
  v3(19,3)=  v3(17,3);
  v3(20,1)=11;
  v3(20,2)=20577;
  v3(20,3)=  v3(17,3);
  v3(21,1)=11;
  v3(21,2)=20910;
  v3(21,3)=  v3(17,3);
  v3(22,1)=11;
  v3(22,2)=26349;
  v3(22,3)=  v3(17,3);
  v3(23,1)=11;
  v3(23,2)=45306;
  v3(23,3)=params(1);
  v3(24,1)=11;
  v3(24,2)=13560;
  v3(24,3)=  v3(23,3);
  v3(25,1)=11;
  v3(25,2)=14189;
  v3(25,3)=  v3(23,3);
  v3(26,1)=11;
  v3(26,2)=20590;
  v3(26,3)=  v3(23,3);
  v3(27,1)=11;
  v3(27,2)=21404;
  v3(27,3)=  v3(23,3);
  v3(28,1)=11;
  v3(28,2)=45121;
  v3(28,3)=  v3(23,3);
  g3 = sparse(v3(:,1),v3(:,2),v3(:,3),27,54872);
end
end
