%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%   Copyright   Xin-Guang Zhu, Yu Wang, Donald R. ORT and Stephen P. LONG
%CAS-MPG Partner Institute for Computational Biology, Shanghai Institutes for Biological Sciences, CAS, Shanghai,200031 
%China Institute of Genomic Biology and Department of Plant Biology, Shanghai Institutes for Biological Sciences, CAS, Shanghai,200031 
%University of Illinois at Urbana Champaign
%Global Change and Photosynthesis Research Unit, USDA/ARS, 1406 Institute of Genomic Biology, Urbana, IL 61801, USA.
 
%   This file is part of e-photosynthesis.
 
%    e-photosynthesis is free software; you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation; 
 
%    e-photosynthesis is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
 
%    You should have received a copy of the GNU General Public License (GPL)
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function Velocity = PRrate(t,PrS,PR_Param)

global NADHc;
global NADc;
global ADPc;
global ATPc;
global GLUc;
global KGc;


global PR_OLD_TIME;
global PR_TIME_N;
global PR_VEL;
global PR_CON;

Gcea = PrS(1);
Gca = PrS(2);
Pga = PrS(3);
Pgca = PrS(4);
Gcac = PrS(5);
Goac = PrS(6);
Serc = PrS(7);
Glyc = PrS(8);
Hprc = PrS(9);
Gceac = PrS(10);
Rubp = PrS(11);
C = PrS(12);
O = PrS(13);

% To set global information for different reactions


% Reaction: 111: RUBP+O2<-->PGlycolate + PGA
global V111;
global KO;
global KC;
global KR;

PrV111 = V111;
PrKO = KO;
PrKC = KC;
PrKR = KR;

% Reaction: 112: PGlycolate-->Pi+Glycolate;
global V112;
global KM112;       % Km112 for PGlycolate;
global KI1122;       % Inhibition constant for glycolate;
global KI1121;      % Inhibition constnat for Pi ( Competitive with PGlycolate)
PrV112 = V112;
PrKM112 = KM112;
PrKI1122 = KI1122;
PrKI1121 = KI1121;

% Reaction 113  : Gla+ATP<-->ADP + PGA
global V113;
global KM1131;  % Km for ATP;
global KM1132;  % Km for Gla;
global KI113;   % Ki for PGA;
global KE113;

PrV113 = V113;
PrKM1131 = KM1131;
PrKM1132 = KM1132;
PrKI113 = KI113;
PrKE113 = KE113;

% Reactoin 121; Glycolate +O2<-->H2O2+Glyoxylate
global V121;
global KM121;
PrV121 = V121;
PrKM121 = KM121;

% Reaction 122  : Glyoxylate + Serine<--> Hydoxypyruvate + Glycine;
global V122;
global KM1221;      % Michaelis constant for glyoxylate;
global KM1222;      % Michaelis constant for serinie;
global KI1221;      % Inhibition constant for Glycine;
global KE122;

PrV122 = V122;
PrKM1221 = KM1221;
PrKM1222 = KM1222;
PrKI1221 = KI1221;
PrKE122 = KE122;

% Reaction 123: HydroxylPyruvate + NAD <--> NADH + Glycerate
global V123;
global KM123;       % Michaelis constant for hydroxylpyruvate;
global KI123;
global KE123;

PrV123 = V123;
PrKM123 = KM123;
PrKI123 = KI123;
PrKE123 = KE123;


% Reaction 124: Glyoxylate + Glu  <--> KG + Glycine;
global V124;
global KM1241;  % Michaelis constant for glyoxylate
global KM1242;  % Michaelis constant for Glu
global KI124;   % This is a guessed value
global KE124;

PrV124 = V124;
PrKM1241 = KM1241;
PrKM1242 = KM1242;
PrKI124 = KI124;
PrKE124 = KE124;

% Reaction 131: LS2+Glycine <--> CO2+ AMDHL
global V131;
global KM1311;   % Michaelis constant for Glycine;
global KM1312;   % Michaelis constant for NAD;
global KI1311;   % Inibition constant for Serine;
global KI1312;   % Inhibition constant for NADH;

PrV131 = V131;
PrKM1311 = KM1311;
PrKM1312 = KM1312;
PrKI1311 = KI1311;
PrKI1312 = KI1312;

% Reaction 132; AMDHL + THF <--> MTHF+DHLA + NH3;
global V132;
global KM1321;  % Michaelis constant for AMDHL;
global KM1322;  % Michaelis constant for THF;
PrV132 = V132;
PrKM1321 = KM1321;
PrKM1322 = KM1322;

% Reaction 133: DHLA + MAD <--> NADH + LS2;
global V133;
global KM1331;  % Michaelis constant for NAD;
global KM1332;  % Michaelis constant for DHLA;
global KI1331;	% Inhibition constant for NADH;
global KI1332;	% Inhibition constant for LS2;

PrV133 = V133;	
PrKM1331 = KM1331;
PrKM1332 = KM1332;
PrKI1331 = KI1331;
PrKI1332 = KI1332;

% The consant for calculating the glycerate uptake.
global V1T;
global KM1011;
global KI1011;
PrKM1011 = KM1011;
PrKI1011 = KI1011;
PrV1T = V1T;

% The constant for calculating the glycolate uptake
global V2T;
global KM1012;
global KI1012;
PrKM1012 = KM1012;
PrKI1012 = KI1012;
PrV2T = V2T;

% Calculate the reactino inside chloroplast;
Pi = PR_Param(2);    % Value from spinach. Assume constant currently.

% Information from PS cycle is transfered back to the PR cycle for calcualtion.
global PR_PS_com;

global PS2PR_Pi;
global CO2_cond;
global O2_cond;

global PR_ADP;         % The chloroplast ADP concentration
global PR_ATP;         % The chloroplast ATP concentration
global PS2PR_ATP;
global PS2PR_ADP;
global PS2PR_Pi;

C = CO2_cond;
O = O2_cond;

if PR_PS_com ==1             % For the combined model. 
    ADP = PS2PR_ADP;
    ATP = PS2PR_ATP;
    Pi =  PS2PR_Pi;
    
    global StomCond_TrDynaPS_com;
    if StomCond_TrDynaPS_com ==1
        global PS2PRC;
        global PS2PRO;
        C = PS2PRC; 
        O = PS2PRO;
    end
else
    ADP = PR_ADP;
    ATP = PR_ATP;
end

global V1Reg;       % This is a parameter generated from PSRate routine.
global RUBISCOMETHOD;
global RUBISCOTOTAL;

global FIBF_PSPR_com;
global EPS_ATP_Rate;
global ATPActive;
global PsV1;   


if FIBF_PSPR_com ==1         
    if ATPActive == 0
        PrV111 = PrV111;
    end
end

if RUBISCOMETHOD ==2          
    if PR_PS_com ==1     
        PrV111t = PrV111*Rubp/(Rubp+PrKR*V1Reg);
    else                  
        PrV111t = PrV111*Rubp/(Rubp+PrKR);
    end
    v111 = PrV111t * O/(O+PrKO*(1+C/PrKC));
    
     if Rubp < PsV1/2.5;
       v111 = v111 * (2.5 * Rubp/PrV111t);
     end 
    
elseif RUBISCOMETHOD==1
    v111 = PrV111 * O/(O+PrKO*(1+C/PrKC));
    if Rubp < RUBISCOTOTAL
        v111 = v111 * Rubp/RUBISCOTOTAL;
    end
end

v112 = PrV112 * Pgca /(Pgca + PrKM112*(1+Gca/PrKI1122)*(1+Pi/PrKI1121));

 if PR_PS_com ==1            
      v113 = PrV113 * (ATP * Gcea - ADP * Pga/PrKE113)/((ATP + PrKM1131*(1 + Pga/PrKI113))*(Gcea + PrKM1132)); % This is the old version.
      
 else
     v113 = PrV113 * (ATP * Gcea - ADP * Pga/PrKE113)/((ATP + PrKM1131*(1 + 2.5/PrKI113))*(Gcea + PrKM1132)); % This is the old version.
 end

Ks = 0.4;       % From Ferjancic-Biagini et al 1998
Ki = 5;         % From Ferjancic-Biagini et al 1998
Kip = 14;       % From Ferjancic-Biagini et al 1998
v121 = PrV121 * Gcac/(Gcac + PrKM121);

% The reaction: Goac + Serine --> HPR + Gly

v122 = PrV122 * (Goac * Serc - Hprc * Glyc/PrKE122)/((Goac+PrKM1221)*(Serc + PrKM1222*(1+Glyc/PrKI1221)));

PrKM1232 = 0.5;    

v123 = PrV123 * (Hprc *NADHc-Gceac * NADc/PrKE123)/((Hprc+PrKM123*(1+Hprc/PrKI123))*(NADHc + PrKM1232));  

v124 = PrV124 * (Goac * GLUc-KGc * Glyc/PrKE124)/((Goac + PrKM1241)*(GLUc + PrKM1242*(1+Glyc/PrKI1221)));

v131 = PrV131 * Glyc/(Glyc + PrKM1311*(1+Serc/PrKI1311)); 

v2out = PrV2T * (Gca/(Gca + PrKM1012*(1+Gcea/PrKI1012)) - Gcac/(Gcac + PrKM1012*(1+Gceac/PrKI1012)));   % Competive inhibition

v1in = PrV1T *(Gceac/(Gceac + PrKM1011 *(1+ Gcac/PrKI1011))-Gcea /(Gcea + PrKM1011*(1+Gca/PrKI1011)));  % Competive inhibition


if (PR_TIME_N ==0)
    PR_TIME_N = 1;
end

if (t > PR_OLD_TIME)
    PR_TIME_N = PR_TIME_N + 1;
    PR_OLD_TIME = t;
end

PR_VEL(PR_TIME_N,1) = t;

PR_VEL(PR_TIME_N,2) = v111;
PR_VEL(PR_TIME_N,3) = v112;
PR_VEL(PR_TIME_N,4) = v113;
PR_VEL(PR_TIME_N,5) = v121;
PR_VEL(PR_TIME_N,6) = v122;
PR_VEL(PR_TIME_N,7) = v123;
PR_VEL(PR_TIME_N,8) = v124;
PR_VEL(PR_TIME_N,9) = v131;
PR_VEL(PR_TIME_N,10) = v1in;
PR_VEL(PR_TIME_N,11) = v2out;

% The following is used to take the information back to the PRmb routine.

Velocity = zeros(10,1);

Velocity(1) = v111;
Velocity(2) = v112;
Velocity(3) = v113;
Velocity(4) = v121;
Velocity(5) = v122;
Velocity(6) = v123;
Velocity(7) = v124;
Velocity(8) = v131;
Velocity(9) = v1in;
Velocity(10) =v2out;


% here is some parameters we need to output
global PR2OUT;
PR2OUT  = zeros(5,1);

PR2OUT(1)    =   Gcea;
PR2OUT(2)    =   Gca;
PR2OUT(3)    =   Pga;
PR2OUT(4)    =   Pgca;
PR2OUT(5)    =   Gcac;
PR2OUT(6)    =   Goac;
PR2OUT(7)    =   Serc;
PR2OUT(8)    =   Glyc;
PR2OUT(9)    =   Hprc;
PR2OUT(10)   =   Gceac;
PR2OUT(11)   =   Rubp;
PR2OUT(12)   =   v131;