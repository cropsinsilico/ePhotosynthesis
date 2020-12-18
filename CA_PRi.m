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

function PrS = CA_PRi(BE)

global NADHc;
global NADc;
global ADPc;
global ATPc;
global GLUc;
global KGc;
global PR_ADP;          
global PR_ATP;        

NADHc = 0.47;        
NADc = 0.4;               

ADPc = 0.64;
ATPc= 0.35;
GLUc=24;
KGc=0.4;

PR_ADP = 0.82;
PR_ATP = 0.68;

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

CE = 1; 
CEV111 = 1;    
CE122 = 1;

    
	% Reaction: 111: RUBP+O2<-->PGlycolate + PGA
	global V111;
	global KO;
	global KC;
	global KR;
	
	global gp2V111;
	V111 = gp2V111;
	KO = 0.48;           % Michaelis constant for O2
	KC = 0.014;          % Michaelis constant for CO2  
	KR = 0.02;           % Michaelis constant for RUBP  
	
	% Reaction: 112: PGlycolate-->Pi+Glycolate;
	global V112;        
	global KM112;       % Km112 for PGlycolate;
	global KI1122;      % Inhibition constant for Glycolate;
	global KI1121;      % The competitive Pi inhibition for PGlycolate
	
	KM112 = 0.026;
	KI1122 = 94;
	KI1121 = 2.55;    
	
	% Reaction 113  : Gcea+ATP<-->ADP + PGA
	global V113;
	global KM1131;  % Km for ATP;
	global KM1132;  % Km for Gcea;
	global KI113;   % Ki for ATP BY pga;
	global  KE113;  % New
	
	KM1131 = 0.21;
	KM1132 = 0.25;
	KI113 = 0.36;                        %%%%%%%%%%%%%%%%%%%%%%%%% Competitive inhibition for ATP; in original paper it is 0.36;
	KE113 = 300;     % New       Kleczkowski et al . 1985 Archives of Biochemistry and Biophysics  
                                           
	
	% Reactoin 121; Glycolate +O2<-->H2O2+Glyoxylate
	global V121;
	global KM121;
	KM121 = 0.1;
	
	% Reaction 122  : Glyoxylate + Serine<--> Hydoxypyruvate + Glycine;
	global V122;
	global KM1221;      % Michaelis constant for glyoxylate;
	global KM1222;      % Michaelis constant for serinie;
	global KI1221;      % Inhibition constant for Glycine;
	global KE122;       % nEW
	
	KM1221 = 0.15;
	KM1222 = 2.7;
	KI1221 = 33;       
	KE122 = 0.24;  %  New: Guynn, R.W.; Arch. Biochem. Biophys.; 218, 14 (1982).; 0.24. At 25 degree. 
	
	% Reaction 123: HydroxylPyruvate + NAD <--> NADH + Glycerate
	global V123;
	global KM123;           %   Michaelis constant for hydroxylpyruvate;
	global KI123;
	global KE123;   % New
	
	KM123 = 0.09;       
	KI123 = 12;          % Inhibition constant for hydroxypyruvate;
	KE123 = 1/(4*10^(-6));  % Guynn, R.W.; Arch. Biochem. Biophys.; 218, 14 (1982).; 1/(4*10^(-6);
	
	
	
	% Reaction 124: Glyoxylate + Glu  <--> KG + Glycine;
	global V124;
	global KM1241;  % Michaelis constant for glyoxylate
	global KM1242;  % Michaelis constant for Glu
	global KI124;   % This KI is one guessed.
	global KE124;   % New       Cooper, A.J.L.; Meister, A.; Biochemistry; 11, 661 (1972).; K' 607. 
	
	
	KM1241 = 0.15;
	KM1242 = 1.7;
	KI124 = 2;        
	KE124 = 607;
	
	% Reaction 131: NAD+Glycine <--> CO2+ NADH + NH3
	global V131;
	global KM1311;   % Michaelis constant for Glycine;
	global KI1311;   % Inhibition constant for Serine
	KM1311 = 6;
	KI1311 = 4;
	
	global KI1312;   % Inhibition constant for NADH;    Since in the current program, we assume that P protein limit the 
                    % rate of the overall glycin decarboxylase; the KI1312 and KM1312 were not used. 
	global KM1312;   % Michaelis constant for NAD;
	KM1312 = 0.075;
	KI1312 = 0.015;
	
	% The consant for calculating the glycerate uptake.
	global V1T;
	global KM1011;
	global KI1011;
	
	V1T = 0.25*CE *10;
	KM1011 = 0.39;
	KI1011 = 0.28;
	
	% The constant for calculating the glycolate uptake
	global V2T;
	global KM1012;
	global KI1012;
	V2T = 0.32*CE * 10;     
	KM1012 = 0.2;
	KI1012 = 0.22;
	
   
global J;
J=1;

global VmaxXValue;

    V111 = VmaxXValue(21,J);   
    V112 = VmaxXValue(22,J);    
    V113 = VmaxXValue(23,J);    
    V121 = VmaxXValue(24,J);    
    V122 = VmaxXValue(25,J);   
    V123 = VmaxXValue(26,J);   
    V124 = VmaxXValue(27,J);       
    V131 = VmaxXValue(28,J);    
    V1T  = VmaxXValue(29,J);    
    V2T  = VmaxXValue(30,J);     
    
    
KO	=	VmaxXValue	(	121	,	J	)	;
KC	=	VmaxXValue	(	122	,	J	)	;
KR	=	VmaxXValue	(	123	,	J	)	;
KM112	=	VmaxXValue	(	124	,	J	)	;
KI1122	=	VmaxXValue	(	125	,	J	)	;
KI1121	=	VmaxXValue	(	126	,	J	)	;
KM1131	=	VmaxXValue	(	127	,	J	)	;
KM1132	=	VmaxXValue	(	128	,	J	)	;
KI113	=	VmaxXValue	(	129	,	J	)	;
KE113	=	VmaxXValue	(	130	,	J	)	;
KM121	=	VmaxXValue	(	131	,	J	)	;
KM1221	=	VmaxXValue	(	132	,	J	)	;
KM1222	=	VmaxXValue	(	133	,	J	)	;
KI1221	=	VmaxXValue	(	134	,	J	)	;
KE122	=	VmaxXValue	(	135	,	J	)	;
KM123	=	VmaxXValue	(	136	,	J	)	;
KI123	=	VmaxXValue	(	137	,	J	)	;
KE123	=	VmaxXValue	(	138	,	J	)	;
KM1241	=	VmaxXValue	(	139	,	J	)	;
KM1242	=	VmaxXValue	(	140	,	J	)	;
KI124	=	VmaxXValue	(	141	,	J	)	;
KE124	=	VmaxXValue	(	142	,	J	)	;
KM1311	=	VmaxXValue	(	143	,	J	)	;
KI1311	=	VmaxXValue	(	144	,	J	)	;
KI1312	=	VmaxXValue	(	145	,	J	)	;
KM1312	=	VmaxXValue	(	146	,	J	)	;
KM1011	=	VmaxXValue	(	147	,	J	)	;
KI1011	=	VmaxXValue	(	148	,	J	)	;
KM1012	=	VmaxXValue	(	149	,	J	)	;
KI1012	=	VmaxXValue	(	150	,	J	)	;

