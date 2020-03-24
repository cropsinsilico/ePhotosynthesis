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



function PrS = PRinitial(BE)
global PRRatio;
global NADHc;
global NADc;
global ADPc;
global ATPc;
global GLUc;
global KGc;
global PR_ADP;        
global PR_ATP;        

NADHc = 0.47*PRRatio(9);        
NADc = 0.4*PRRatio(10);         

ADPc = 0.64*PRRatio(11);
ATPc= 0.35*PRRatio(12);
GLUc=24*PRRatio(13);
KGc=0.4*PRRatio(14);

PR_ADP = 0.82*PRRatio(15);
PR_ATP = 0.68*PRRatio(16);

SERc= 7.5;                  % Serine in cytosol; 7.5 original value
GLYc = 1.8;                 % Glycine in cytosol; 1.8 original vlaue
PGA = 4.3;                  % PGA in chloroplast;4.3 is the original value;

GOAc = 0.028;              % Glyoxylate in cytosol; 0.028; EXPERIMENTAL DATA;

GCAc = 0.36;                   % See the note for GCA.
GCA = 0.36;                    % Derived from radioactive labelling experiment; assuem equal concenatration 
                               % inside and outshide chloroplast

PGCA= 0.0029;                % Phosphoglycolate in chloroplast derived based on the Km112; orignal value is : 0.0029; 
GCEA =0.1812;                  % Glycerate in chloroplast; derived based on V113
GCEAc = 0.1812;                 % Glycerate in cytosol; assume at equilibrium with GCEA initially.
HPRc = 0.0035;                % HydroxylPyruvate; derived from equation 123;
RUBP = 2;                   % RuBP concentration

CO2 = 0.012;                 % CO2 concentration(mM)
O2 = 0.264;                  % O2 concentration(mM)

PrS = zeros(10,1);

PrS(1) = GCEA;
PrS(2) = GCA;
PrS(3) = PGA;
PrS(4) = PGCA;

PrS(5) = GCAc;
PrS(6) = GOAc;
PrS(7) = SERc;
PrS(8) = GLYc;
PrS(9) = HPRc;
PrS(10) = GCEAc;
PrS(11) = RUBP;
PrS(12) = CO2;
PrS(13) = O2;

% To set global information for different reactions
% Reaction: 110: RuBP + CO2 <--> 2PGA

CE = 1; % This is the coefficient for calibrating the volume effect % Default is 4. 
CEV111 = 1;    % 1.72 was used to 
CE122 = 1;

% Reaction: 111: RUBP+O2<-->PGlycolate + PGA
global V111;
global KO;
global KC;
global KR;
global PR_PS_com;
global PS2PR_V1;


KO = 0.222*PRRatio(17);           % Michaelis constant for O2
KC = 0.0115*PRRatio(18);          % Michaelis constant for CO2  

if PR_PS_com ==1
    global	KM11	;
    global	KM12	;
    KC = KM11;
    KO = KM12; 
end

KR = 0.02*PRRatio(19);           % Michaelis constant for RUBP  

% Reaction: 112: PGlycolate-->Pi+Glycolate;
global V112;        
global KM112;       % Km112 for PGlycolate;
global KI1122;      % Inhibition constant for Glycolate;
global KI1121;      % The competitive Pi inhibition for PGlycolate

KM112 = 0.026*PRRatio(20);
KI1122 = 94*PRRatio(21);
KI1121 = 2.55*PRRatio(22);    


% Reaction 113  : Gcea+ATP<-->ADP + PGA
global V113;
global KM1131;  % Km for ATP;
global KM1132;  % Km for Gcea;
global KI113;   % Ki for ATP BY pga;
global  KE113;  % New
KM1131 = 0.21*PRRatio(23);
KM1132 = 0.25*PRRatio(24);
KI113 = 0.36*PRRatio(25);    %%%%%%%%%%%%%%%%%%%%%%%%% Competitive inhibition for ATP; in original paper it is 0.36;
KE113 = 300*PRRatio(26);     % New       Kleczkowski et al . 1985 Archives of Biochemistry and Biophysics  300, as default
                                       

% Reactoin 121; Glycolate +O2<-->H2O2+Glyoxylate
global V121;
global KM121;

KM121 = 0.1*PRRatio(27);

% Reaction 122  : Glyoxylate + Serine<--> Hydoxypyruvate + Glycine;
global V122;
global KM1221;      % Michaelis constant for glyoxylate;
global KM1222;      % Michaelis constant for serinie;  
global KI1221;      % Inhibition constant for Glycine; 
global KE122;       % New
KM1221 = 0.15*PRRatio(28);
KM1222 = 2.7*PRRatio(29);
KI1221 = 33*PRRatio(30);       
KE122 = 0.24*PRRatio(31);   %  New: Guynn, R.W.; Arch. Biochem. Biophys.; 218, 14 (1982).; 0.24. At 25 degree. 

% Reaction 123: HydroxylPyruvate + NAD <--> NADH + Glycerate
global V123;
global KM123;           %   Michaelis constant for hydroxylpyruvate;
global KI123;
global KE123;   % New

KM123 = 0.09*PRRatio(32);       
KI123 = 12*PRRatio(33);          % Inhibition constant for hydroxypyruvate;
KE123 = 1/(4*10^(-6))*PRRatio(34);  % Guynn, R.W.; Arch. Biochem. Biophys.; 218, 14 (1982).; 1/(4*10^(-6);

% Reaction 124: Glyoxylate + Glu  <--> KG + Glycine;
global V124;
global KM1241;  % Michaelis constant for glyoxylate
global KM1242;  % Michaelis constant for Glu
global KI124;   % This KI is one guessed.
global KE124;   % New       Cooper, A.J.L.; Meister, A.; Biochemistry; 11, 661 (1972).; K' 607. 

  
KM1241 = 0.15*PRRatio(35);
KM1242 = 1.7*PRRatio(36);
KI124 = 2*PRRatio(37);     % This is a guessed vlaue        ???????????????? To be calibrated.
KE124 = 607*PRRatio(38);

% Reaction 131: NAD+Glycine <--> CO2+ NADH + NH3
global V131;
global KM1311;   % Michaelis constant for Glycine;
global KI1311;   % Inhibition constant for Serine

KM1311 = 6*PRRatio(39);
KI1311 = 4*PRRatio(40);

global KI1312;   % Inhibition constant for NADH;    Since in the current program, we assume that P protein limit the 
                % rate of the overall glycin decarboxylase; the KI1312 and KM1312 were not used. 
global KM1312;   % Michaelis constant for NAD;
KM1312 = 0.075*PRRatio(41);
KI1312 = 0.015*PRRatio(42);

% The consant for calculating the glycerate uptake.
global V1T;
global KM1011;
global KI1011;

V1T = 0.25*CE *20*PRRatio(43);
KM1011 = 0.39*PRRatio(44);
KI1011 = 0.28*PRRatio(45);

% The constant for calculating the glycolate output
global V2T;
global KM1012;
global KI1012;
V2T = 0.32*CE * 10 *2*PRRatio(46) ;  
KM1012 = 0.2*PRRatio(47);
KI1012 = 0.22*PRRatio(48);


global gp2V111; 
V111 = gp2V111; 

global GP; 
if GP ==0
	if PR_PS_com ==1
        V111 = 0.24 * PS2PR_V1*PRRatio(1);
	
	else
        V111 = 3.7*0.24 * 1; 
	end
	V112 = 52.41992121*PRRatio(2);         
	V113 = 5.715787563*PRRatio(3) ;  
	V121 = 1.456108923*PRRatio(4);    
	V122 = 3.306190845 * 3*PRRatio(5) ;    
	V123 = 10.00978112*PRRatio(6);         
	V124 = 2.745819515*PRRatio(7);   
	V131 = 2.494745448 *PRRatio(8);   
end

