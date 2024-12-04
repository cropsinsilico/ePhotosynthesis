function out = VARS_Dict()

names = [
"AVR", ...
"HPR", ...
"O2", ...
"CO2_cond", ... %"GLight", ...
"O2_cond", ...
"PS12ratio", ...
"ADP", ...
"Pi", ...
"TestATPCost", ...
"Air_CO2", ...
"CO2_in", ...
"Radiation_PAR", ...
"PPFD_in", ...
"TestLi", ...
"PS2BF_Pi", ...
"PS_PR_Param", ...
"Tp", ...
"alfa", ...
"fc", ...
"lightParam", ...
"alpha1", ...
"alpha2", ...
"input_LHCI", ...
"input_LHCII", ...
"input_PSI", ...
"input_PSIIcore", ...
"PSIIantennaSize", ...
"PSIantennaSize", ...
];

global AVR;
global HPR;
global CO2_cond;
%global GLight;
global O2_cond;
global PS12ratio;
global PS2BF_ADP;
global CO2_Env;
global LI;
global PS2BF_Pi;
global Tp;
global alfa;
global fc;
global input_LHCI;
global input_LHCII;
global input_PSI;
global input_PSIIcore;
global PSIIantennaSize;
global PSIantennaSize;

Air_CO2 = CO2_Env / 0.7;
PPFD_in = LI * 30;
Radiation_PAR = PPFD_in / (1E6/(2.35E5));

values = [
AVR, ...
HPR, ...
210.0, ... % O2 is CM_Con(12), 210 in C++, but unused
CO2_cond, ... %GLight, ...
O2_cond, ...
PS12ratio, ...
PS2BF_ADP, ...
PS2BF_Pi, ... % Pi is BF_con(16), 0 in C++, only used for C3, PS2BF_Pi for EPS
0.0, ... % TestATPCost not part of matlab code, only used w/o C3
Air_CO2, ...
CO2_Env, ...
Radiation_PAR, ...
PPFD_in, ...
LI, ...
PS2BF_Pi, ...
0.0, ... % PS_PR_Param unused by C++ or matlab
Tp, ...
alfa, ...
fc, ...
LI*30*alfa*(1.0-fc), ...
1.0, ... % alpha1 unused by matlab
1.0, ... % alpha2 unused by matlab
input_LHCI, ...
input_LHCII, ...
input_PSI, ...
input_PSIIcore, ...
PSIIantennaSize, ...
PSIantennaSize, ...
];

out_VARS = dictionary(names, values);
%out_VARS = remove(out_VARS, "GLight");

out = dictionary();
out('VARS') = out_VARS;

end
