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



global StepSize;
global VmaxMatrix;
global VmaxVary;
global MetConCoeff;;
global DriveSelector;
global NewCalculate;
NewCalculate =1;


StepSize = 5*10^-2;                       %   Distance used to calculate the tangent
                                        %    of the graph for some given Vmax value.
DriveSelector = 4;                      %   User determines which simulation to use:
                                        %    1 = CA_PSDrive.m
                                        %    2 = CA_PS_PRDrive.m
                                        %    3 = CA_trDynaPS_Drive.m
                                        %    4 = CA_CMDrive.m
global ControlAnalysis
ControlAnalysis = 3;

global InExcel;
global UseInputVmax;
UseInputVmax = 0;               % 1, used the V from InExcel; 0 used the default values. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The center values for each Vmax is chosen here, along with whether or not
% one wishes to vary the value of a particular Vmax.  If one wishes to vary
% Vn for some n, then put VnVary = 1, otherwise, enter 0.  Only one Vn may
% be varied at a time.  If more than one VnVary = 1, CA.m will not run.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

V1		= 3.38;                         %   Default is 7.46

V1Vary  = 1;

V2		= 73.4;                        %   Default is 23.52
V2Vary  = 1;

V3		=9.8;                        %   Default is 10.08
V3Vary  = 0;

V4      = 0;                            %   Default is 0
V4Vary  = 0;                            %   Leave at 0

V5		= 2.96;                        %   Default is 11.04       
V5Vary  = 0;

V6		= 1.76;                         %   Default is 2.31
V6Vary  = 0;

V7		= 7.605;                        %   Default is 2.016
V7Vary  = 0;

V8		= 2.96;                        %   Default is 11.04
V8Vary  = 0;


V9		= 0.79;                        %   Default is 0.336
V9Vary  = 0;

V10		= 7.605;                        %   Default is 2.016
V10Vary  = 0;

V11		= 0;                            %   Default is 0
V11Vary  = 0;

V12		= 0;                            %   Default is 0
V12Vary  = 0;

V13		= 26.39;                      %   Default is 16.0188
V13Vary  = 0;

V16		= 15;                            %   Default is 6
V16Vary  = 0; 

V21		= 0;                            %   Default is 0
V21Vary  = 0;

V22		= 0;                            %   Default is 0                          
V22Vary  = 0;

V23		= 0.65;                         %   Default is 3.36
V23Vary  = 0;

V31		= 3.73/3;                            %   Default is 0
V31Vary  = 0;

V32		= 3.73/3;                            %   Default is 0
V32Vary  = 0;

V33		= 3.73/3;                            %   Default is 0
V33Vary  = 0;

V111    = V1 * 0.2;                       %   Default is 0.0466    
V111Vary = 0;

V112    = 71;                         %   Default is 3.27                
V112Vary = 0;

V113    = 0.93;                        %   Default is 0.457    
V113Vary = 0;

V121    = 1.18;                          %   Default is 1.4  
V121Vary = 0;

V122    = 0.54;                          %   Default is 1.6  
V122Vary = 0;

V123    = 8.13;                          %   Default is 7.3  
V123Vary = 0;

V124    = 0.4465;                         %   Default is 0.55   
V124Vary = 0;

V131    = 0.225;                          %   Default is 0.7  
V131Vary = 0;


V1T     = 0.25;                         %   Default is 0.25   
V1TVary  = 0;

V2T     = 0.32;                         %   Default is 0.32   
V2TVary  = 0;

FIBF_PQT = 8;                           %   Default is 8
FIBF_PQTVary  = 0;

ChlT = 70;                              %   Default is 70
ChlTVary  = 0;

ChlT2 = 290;                            %   Default is 290
ChlT2Vary  = 0;

V2M = 27.8;                             %   Default is 27.8
V2MVary  = 0;

Vmax11 = 3;                             %   Default is 3
Vmax11Vary  = 0;

Tcyt = 1;                               %   Default is 1
TcytVary  = 0;

BF_An = 120;                            %   Default is 120
BF_AnVary  = 0;

BF_Cn = 80;                             %   Default is 80
BF_CnVary  = 0;



KM11	=	0.0115;		% 	CO2	1	RuBP+CO2->2PGA
KM12	=	0.222;		%	O2	1	RuBP+CO2->2PGA
KM13	=	0.02;		% 	RuBP	1	RuBP+CO2->2PGA
KI11    =  	0.84   ;    % PGA						
KI12    =	0.04   ;       % FBP						
KI13    = 	0.075 ;       % SBP						
KI14    = 	0.9   ;       % Pi						
KI15    = 	0.07  ;       % NADPH						
KM21	=	0.240;		%	PGA	2	PGA+ATP <-> ADP + DPGA
KM22	=	0.390;		% 	ATP	2	PGA+ATP <-> ADP + DPGA
KM23    =   	0.23  ;        %  ADP     						
KM31a	=	0.004;		%	BPGA	3	DPGA+NADPH <->GAP + OP+NADP 
KM32b	=	0.1	;	    % 	NADPH	3	DPGA+NADPH <->GAP + OP+NADP
KM41	=	2.5	;	    %	DHAP	4	DHAP <->GAP
KM42	=	0.68;		% 	GAP	4	DHAP <->GAP
KE4     =   	0.05;       %   Using the value from Patterson						
KM51	=	0.3	;	    %	GAP	5	GAP+DHAP <->FBP
KM52	=	0.4	;	    % 	DHAP	5	GAP+DHAP <->FBP
KM53	=	0.02;		%	FBP	5	GAP+DHAP <->FBP     % Original Value: 0.02
KE5     = 	7.100;          % Defult: 7.1						
KM61	=	0.033;		% 	FBP	6	FBP<->F6P+OP
KI61    = 	0.7   ;       %   F6P       						
KI62    = 	12    ;       %   Pi						
KE6     =   	6.66 * 10^5;    % The equilibrium constant for this reaction        % New    mM     Laisk or Bassham and Krause 1969 BBA						
KM71	=	0.100;		%	Xu5P	7	F6P+GAP<->E4P+Xu5P      % jn
KM72	=	0.100;		% 	E4P	7	F6P+GAP<->E4P+Xu5P
KM73    = 	0.1;         %   F6P This value was based on estimate						
KM74    = 	0.1000;         % Estimate for GAP ORIGINAL 0.1						
KE7     =   	10 ;       % The equilibrium constant for this reaction             % New           Laisk  Bassham and Krause 1969 BBA						
KM8	    =	0.02;		%	SBP	8	E4P+DHAP<->SBP
KM81    = 	0.4   ;       % DHAP						
KM82    = 	0.2   ;       % E4P estimate						
KE8     = 	1.017 ;     % The equilibrium constant for this reaction                  % New    mM-1         Laisk  Bassham and Krause 1969 BBA. Default: 1.107								
KM9	=	0.05;		% 	SBP	9	SBP<->S7P+OP    		
KI9     = 	12    ;       %   The inibintion constant for Pi; 								
KE9     =   	6.66 * 10^5 ; % The equilibrium constant of this reaction           % New   mM      Laisk  Bassham and Krause 1969 BBA								
KM10	=	1.5	;	    %	R5P	10	S7P+GAP<->Ri5P+Xu5P		
KM101   =   	0.1 ;       %   Xu5P								
KM102   = 	0.072 ;       %   Estimate for GAP								
KM103   = 	0.46 ;        %   Estimate for S7P                                    % New 								
KE10    = 	1/0.85 ;      %   The equilibrium constant for this reaction          % New From Laisk or Bassham and Krause 1969 BBA								
KE11	=	0.4	;	    %	Equilibrium Constant	11	Ri5P<-->Ru5P		
KE12	=	0.67;		% 	Equilibrium Constant	12	Xu5P<-->Ru5P		
KM131	=	0.05;		    %	Ru5P	13	Ru5P+ATP<->RuBP+ADP		
KM132	=	0.059;		    % 	ATP	13	Ru5P+ATP<->RuBP+ADP		
KI131	=	2	;			%	PGA	13	Ru5P+ATP<->RuBP+ADP
KI132	=	0.7	;			%	RuBP	13	Ru5P+ATP<->RuBP+ADP
KI133	=	4	;			%	Pi	13	Ru5P+ATP<->RuBP+ADP
KI134	=	2.5	;			%	ADP	13	Ru5P+ATP<->RuBP+ADP
KI135	=	0.4	;			%	ADP	13	Ru5P+ATP<->RuBP+ADP
KE13    =   	6.846 * 10^3;   %   The equilibrium constant for this reaction  % New From Laisk or Bassham and Krause 1969 BBA								
KM161	=	0.014;		%	ADP	16	ADP+Pi<->ATP 		
KM162	=	0.3;		% 	Pi	16	ADP+Pi<-> ATP		
KM163   =   	0.3;        %   ATP 16  ADP+Pi<-> ATP                           % New       Based on Laisk  								
KE16    =   	5.734;      %   The equilibrium constant for this reaction      % NEW, From Laisk or Bassham and Krause 1969 BBA								
KE21	=	2.3;		%	Equilibrium constant	21	F6P<->G6P		
KE22	=	0.058;		% 	Equilibrium constant	22	G6P<->G1P		
KM231	=	0.08;		%	G1P	23	G1P+ATP+Gn<->PPi+ADP+Gn+1		
KM232	=	0.08;		% 	ATP	23	G1P+ATP+Gn<->PPi+ADP+Gn+1		
KA231	=	0.1;		%	PGA	23	G1P+ATP+Gn<->PPi+ADP+Gn+1		
KA232	=	0.02;		% 	F6P	23	G1P+ATP+Gn<->PPi+ADP+Gn+1		
KA233	=	0.02;		%	FBP	23	G1P+ATP+Gn<->PPi+ADP+Gn+1		
KI23	=	10;		    % 	ADP	23	G1P+ATP+Gn<->PPi+ADP+Gn+1		
KM311	=	0.077;		%	DHAP	31	DHAPi<->DHAPo		
KM312	=	0.63;		% 	Pi	31	DHAPi<->DHAPo
KM313	=	0.74;		%	Pext	31	DHAPi<->DHAPo
KM32	=	0.25;		% 	PGA	32	PGAi<->PGAo
KM33	=	0.075;		%	GAP	33	GAPi<->GAPo

CE=1							
KO = 0.48;           % Michaelis constant for O2							
KC = 0.014;          % Michaelis constant for CO2  							
KR = 0.02;           % Michaelis constant for RUBP  							
KM112 = 0.026;							
KI1122 = 94;							
KI1121 = 2.55;    							
KM1131 = 0.21;							
KM1132 = 0.25;							
KI113 = 0.36;                        %%%%%%%%%%%%%%%%%%%%%%%%% Competitive inhibition for ATP; in original paper it is 0.36;							
KE113 = 300;     % New       Kleczkowski et al . 1985 Archives of Biochemistry and Biophysics  							
KM121 = 0.1;
KM1221 = 0.15;
KM1222 = 2.7;
KI1221 = 33;       
KE122 = 0.24;  %  New: Guynn, R.W.; Arch. Biochem. Biophys.; 218, 14 (1982).; 0.24. At 25 degree. 
KM123 = 0.09;       
KI123 = 12;          % Inhibition constant for hydroxypyruvate;
KE123 = 1/(4*10^(-6));  % Guynn, R.W.; Arch. Biochem. Biophys.; 218, 14 (1982).; 1/(4*10^(-6);
KM1241 = 0.15;
KM1242 = 1.7;
KI124 = 2;     % This is a guessed vlaue        ???????????????? To be calibrated.
KE124 = 607;
KM1311 = 6;
KI1311 = 4;
KM1312 = 0.075;
KI1312 = 0.015;
KM1011 = 0.39;
KI1011 = 0.28;
KM1012 = 0.2;
KI1012 = 0.22;


KE501	=	0.05	;	%	Equilibrium Constant		50		KE501		0.05		[Bassham, 1869 #832]
Km511	=	0.02	;	%	FBP	4.1.2.13	51		Km511	FBP	0.02	Pisum sativum	(Anderson, Heinrikson et al. 1975)
Km512	=	0.3	;	%	FBP	4.1.2.13	51		Km512	GAP	0.3	Spinacia oleracea	(Iwaki, Wadano et al. 1991)
Km513	=	0.4	;	%	FBP	4.1.2.13	51		Km513	DHAP	0.4	Spinacia oleracea	(Iwaki, Wadano et al. 1991)
KE51	=	12	;	%	Based on Thomas et al 1997 Biochem Journal. The fifth citation in the paper. 								
Km514	=	0.014	;	%	FBP	4.1.2.13	51		Km514	SBP	0.014	Spinacia oleracea	(Harris and Koniger 1997)
Km521	=	0.0025	;	%	FBPase[1]	3.1.3.11	52		Km521	FBP	0.0025	Pisum sativum	(Jang, Lee et al. 2003)
KI521	=	0.7	;	%	FBPase	3.1.3.11	52		KI521	F6P	0.7		[Heldt, 1983 #841]
KI522	=	12	;	%	FBPase	3.1.3.11	52		KI522	Pi	12	Pisum sativum	(Charles & Halliwell 1997)
KI523	=	7*10^(-5)	;	%	FBPase	3.1.3.11	52		KI523	F26BP	7*10^(-5)	Pisum sativum <Com>	{Jang, 2003 #2523}
KE52	=	6663	;	%	FBPase	3.1.3.11	52		KE52			6663	[Bassham, 1869 #832]
KE531	=	2.3	;	%	Equilibrium Constant	5.3.1.9	53		KE531		2.3[2]		[Bassham, 1869 #832]
KE541	=	0.0584	;	%	Equilibrium Constant	5.4.2.2	54	G1P G6P	KE541	G1P G6P	0.0584		[Bassham, 1869 #832]
Km551	=	0.14	;	%	UGPase	2.7.7.9	55		Km551	G1P	0.14	Solanum tuberosum	(Nakano, Omura et al. 1989)
Km552	=	0.1	;	%	UDPase	2.7.7.9	55		Km552	UTP	0.1	Solanum tuberosum	(Nakano, Omura et al. 1989)
Km553	=	0.11	;	%	UGPase	2.7.7.9	55		Km553	OPOP	0.11	Solanum tuberosum	(Nakano, Omura et al. 1989)
Km554	=	0.12	;	%	UGPase	2.7.7.9	55		Km554	UDPGlu	0.12	Solanum tuberosum	(Nakano, Omura et al. 1989)
KE55	=	0.31	;	%	UGPase	2.7.7.9	55		KE55	Equi	0.31		Lunn and Rees 1990
Km561	=	0.8	;	%	SPase	2.4.1.14	56		Km561	D-F6P	0.8	Pisum sativum	(Lunn and Ap Rees 1990)
Km562	=	2.4	;	%	Spase	2.4.1.14	56		Km562	UDP-glucose	2.4	Pisum sativum	(Lunn and Ap Rees 1990)
KI561	=	0.7	;	%				Inhibitor	KI561	UDP	0.7	Spinacia oleracea	(Harbron, Foyer et al. 1981)
KI562	=	0.8	;	%	Sucrose Synthesase			Inhibitor	KI562	FBP	0.8	Spinacia oleracea	(Harbron, Foyer et al. 1981)
KI563	=	0.4	;	%				Inhibitor	KI563	SUCP	0.4	Spinacia oleracea	(Harbron, Foyer et al. 1981)
KI564	=	11	;	%		2.4.1.14	56	Inhibitor	KI564	Pi	11	Spinacia oleracea	(Harbron, Foyer et al. 1981)
KI565	=	50	;	%		2.4.1.14	56	Inhibitor	KI565	Sucrose	50	Spinacia oleracea	{Salerno, 1978 #2525}
KE56	=	10	;	%					KE56		10	Pisum sativum	Lunn and Rees, 1990
Km571	=	0.35	;	%	SPP	3.1.3.24	57.1		Km571	SUCP	0.35	Pisum sativum	(Whitaker 1984)
Ki572	=	80	;	%	SPP	3.1.3.24	57.2		Ki572	SUC	80	Daucus carota	(Whitaker 1984)
KE57	=	780	;	%	SPP	3.1.3.24	57.2		KE57	Equili	780		Lunn and Rees 1990
Km581	=	0.032	;	%	F26BPa	3.1.3.46	58		Km581	F26BP	0.032	Spinacia oleracea	(Macdonald, Chou et al. 1989)
KI581	=	0.1	;	%	F26BPa	3.1.3.46	58		KI581	F6P	0.1	Arabidopsis thaliana	(Villadsen and Nielsen 2001)
KI582	=	0.5	;	%	F26BPa	3.1.3.46	58		KI582	OP	0.5	Arabidopsis thaliana	(Villadsen and Nielsen 2001)
Km591	=	0.5	;	%	6PF2K	2.7.1.105	59		Km591	ATP	0.5	Spinacia oleracea	(Walker and Huber 1987)
Km592	=	0.021	;	%	6PF2K	2.7.1.105	59		Km592	F26BP	0.021	Sparus aurate	(Garcia de Frutos and Baanante 1995)
Km593	=	0.5	;	%	6PF2K	2.7.1.105	59		Km593	F6P	0.5	Spinacia oleracea	(Walker and Huber 1987)
KI591	=	0.16	;	%			59		KI591	ADP	0.16	Rattus norvegicus	(Kretschmer and Hofmann 1984)
KI592	=	0.7	;	%	6PF2K	2.7.1.105	59		KI592	DHAP	0.7	Spinacia oleracea	{Markham, 2002 #2524}
KE59	=	590	;	%	6PF2K	2.7.1.105	59		KE59		590		Cornish-Bowden, 1997
Km601	=	0.042	;	%	Nucleoside Diphosphate Kinase	2.7.4.6	60	NI	Km601	ADP	0.042	Rat	Kamura and Shimada 1988
Km602	=	1.66	;	%	Nucleoside Diphosphate Kinase	2.7.4.6	60	NI	Km602	ATP	1.66	Rat	Kamura and Shimada 1988
Km603	=	0.28	;	%	Nucleoside Diphosphate Kinase	2.7.4.6	60	NI	Km603	UDP	0.28	Saccharomyces cerevisiae	{Jong, 1991 #2518}
Km604	=	16	;	%	Nucleoside Diphosphate Kinase	2.7.4.6	60	NI	Km604	UTP	16	Rattus norvegicus	{Fukuchi, 1994 #2519}
KE60	=	16	;	%	Nucleoside Diphosphate Kinase	2.7.4.6	60	NI	KE60		16	1.04	{Lynn, 1978 #2520}
KE61	=	1200	;	%	Pyrophosphate hydrolysis				KE61		1.2*107		{Flodgaard, 1974 #2521}
Km621	=	5	;	%	Vsink			Notice: pH dependent	Km621	Sucrose	5		{Weschke, 2000 #2522}

%V51	=	2.5	;	%	DHAP+GAP --FBP          % default 0.5								
%V52	=	2	;	%	FBP --F6P + Pi								
%V55	=	50	;	%	G1P+UTP --OPOP+UDPG 								
%V56	=	8.5	;	%	UDPG+F6P--SUCP + UDP
%V57	=	0.27	;	%	SUCP--Pi + SUC; This 0.01 is to enforce phosphate limitation
%V58	=	0.03	;	%	F26BP--F6P + Pi
V51	=	0.158116716	;
V52	=	0.09421173	;
V55	=	0.169935876	;
V56	=	0.081731064	;
V57	=	0.817310639	;
V58	=	0.024766989	;

V59	=	0.03	;	%	F6P + ATP --ADP + F26BP
V60	=	6.1	;	%	ATP+UDP --UTP + ADP
V61	=	10000	;	%	
V62	=	2	;	%	SUC Sink        0.2 works.
Vdhap_in	=	1.2	;	%	
Vgap_in 	=	1.2	;	%	
Vpga_in  	=	1.2	;	%	

% This section is copied directly from the gpmain

V1	=	2.913930914	;
V2	=	30.1408176	;
V3	=	4.039482839	;
V5	=	1.218890457	;
V6	=	0.726259575	;
V7	=	3.122215944	;
V8	=	1.218890364	;
V9	=	0.324190582	;
V10	=	3.122068837	;
V13	=	10.83475664	;
V23	=	0.26684349	;
V112	=	52.41992121	;
V113	=	5.715787563	;
V121	=	1.456108923	;
V122	=	3.306190845	;
V123	=	10.00978112	;
V124	=	2.745819515	;
V131	=	2.494745448	;
V51	=	0.107376831	;
V52	=	0.063979048	;
V55	=	0.115403205	;
V56	=	0.055503446	;
V57	=	0.55503446	;
V58	=	0.016819226	;


V16 = 15;
V31 = 3.73/3;
V32 = 3.73/3;
V33 = 3.73/3;
V111    = V1 * 0.24;                       %   Default is 0.0466    

KM11_vary	=	0	;
KM12_vary	=	0	;
KM13_vary	=	0	;
KI11_vary	=	0	;
KI12_vary	=	0	;
KI13_vary	=	0	;
KI14_vary	=	0	;
KI15_vary	=	0	;
KM21_vary	=	0	;
KM22_vary	=	0	;
KM23_vary	=	0	;
KM31a_vary	=	0	;
KM32b_vary	=	1	;
KM41_vary	=	0	;
KM42_vary	=	0	;
KE4_vary	=	0	;
KM51_vary	=	0	;
KM52_vary	=	0	;
KM53_vary	=	0	;
KE5_vary	=	0	;
KM61_vary	=	0	;
KI61_vary	=	0	;
KI62_vary	=	0	;
KE6_vary	=	0	;
KM71_vary	=	0	;
KM72_vary	=	0	;
KM73_vary	=	0	;
KM74_vary	=	0	;
KE7_vary	=	0	;
KM8_vary	=	0	;
KM81_vary	=	0	;
KM82_vary	=	0	;
KE8_vary	=	0	;
KM9_vary	=	0	;
KI9_vary	=	0	;
KE9_vary	=	0	;
KM10_vary	=	0	;
KM101_vary	=	0	;
KM102_vary	=	0	;
KM103_vary	=	0	;
KE10_vary	=	0	;
KE11_vary	=	0	;
KE12_vary	=	0	;
KM131_vary	=	0	;
KM132_vary	=	0	;
KI131_vary	=	0	;
KI132_vary	=	0	;
KI133_vary	=	0	;
KI134_vary	=	0	;
KI135_vary	=	0	;
KE13_vary	=	0	;
KM161_vary	=	0	;
KM162_vary	=	0	;
KM163_vary	=	0	;
KE16_vary	=	0	;
KE21_vary	=	0	;
KE22_vary	=	0	;
KM231_vary	=	0	;
KM232_vary	=	0	;
KA231_vary	=	0	;
KA232_vary	=	0	;
KA233_vary	=	0	;
KI23_vary	=	0	;
KM311_vary	=	0	;
KM312_vary	=	0	;
KM313_vary	=	0	;
KM32_vary	=	0	;
KM33_vary	=	0	;
			
			
KO_vary	=	0	;
KC_vary	=	0	;
KR_vary	=	0	;
KM112_vary	=	0	;
KI1122_vary	=	0	;
KI1121_vary	=	0	;
KM1131_vary	=	0	;
KM1132_vary	=	0	;
KI113_vary	=	0	;
KE113_vary	=	0	;
KM121_vary	=	0	;
KM1221_vary	=	0	;
KM1222_vary	=	0	;
KI1221_vary	=	0	;
KE122_vary	=	0	;
KM123_vary	=	0	;
KI123_vary	=	0	;
KE123_vary	=	0	;
KM1241_vary	=	0	;
KM1242_vary	=	0	;
KI124_vary	=	0	;
KE124_vary	=	0	;
KM1311_vary	=	0	;
KI1311_vary	=	0	;
KI1312_vary	=	0	;
KM1312_vary	=	0	;
KM1011_vary	=	0	;
KI1011_vary	=	0	;
KM1012_vary	=	0	;
KI1012_vary	=	0	;


KE501_vary	=	0	;
Km511_vary	=	0	;
Km512_vary	=	0	;
Km513_vary	=	0	;
KE51_vary	=	0	;
Km514_vary	=	0	;
Km521_vary	=	0	;
KI521_vary	=	0	;
KI522_vary	=	0	;
KI523_vary	=	0	;
KE52_vary	=	0	;
KE531_vary	=	0	;
KE541_vary	=	0	;
Km551_vary	=	0	;
Km552_vary	=	0	;
Km553_vary	=	0	;
Km554_vary	=	0	;
KE55_vary	=	0	;
Km561_vary	=	0	;
Km562_vary	=	0	;
KI561_vary	=	0	;
KI562_vary	=	0	;
KI563_vary	=	0	;
KI564_vary	=	0	;
KI565_vary	=	0	;
KE56_vary	=	0	;
Km571_vary	=	0	;
Ki572_vary	=	0	;
KE57_vary	=	0	;
Km581_vary	=	0	;
KI581_vary	=	0	;
KI582_vary	=	0	;
Km591_vary	=	0	;
Km592_vary	=	0	;
Km593_vary	=	0	;
KI591_vary	=	0	;
KI592_vary	=	0	;
KE59_vary	=	0	;
Km601_vary	=	0	;
Km602_vary	=	0	;
Km603_vary	=	0	;
Km604_vary	=	0	;
KE60_vary	=	0	;
KE61_vary	=	0	;
Km621_vary	=	0	;
V51_vary	=	0	;
V52_vary	=	0	;
V55_vary	=	0	;
V56_vary	=	0	;
V57_vary	=	0	;
V58_vary	=	0	;
V59_vary	=	0	;
V60_vary	=	0	;
V61_vary	=	0	;
V62_vary	=	0	;
Vdhap_in_vary	=	0	;
Vgap_in_vary	=	0	;
Vpga_in_vary	=	0	;


if UseInputVmax ==1
    V1 = InExcel(3);
    V2 = InExcel(4);
    V3 = InExcel(5);
    V5 = InExcel(6);
    V6 = InExcel(7);
    V7 = InExcel(11);
    V8 = InExcel(9);
    V9 = InExcel(10);
    V10 = InExcel(11);
    V13 = InExcel(12);
    V16 = 15;               % Assumed in this current study. In the future, this line need to be changed into : V16 = InExcel(13); 
    V23 = InExcel(14);
    V112 = InExcel(15);
    V113 = InExcel(16);
    V121 = InExcel(17);
    V122 = InExcel(18);
    V123 = InExcel(19);
    V124 = InExcel(20);
    V131 = InExcel(21);
    V51 = InExcel(22);
    V52 = InExcel(23);
    V55 = InExcel(24);
    V56 = InExcel(25);
    V57 = InExcel(26);
    V58 = InExcel(27);
end




VmaxMatrix  = zeros(30,1);    
VmaxVary    = zeros(30,1);              

VmaxMatrix(1,1)  = V1;
VmaxVary(1,1)    = V1Vary;
VmaxMatrix(2,1)  = V2;
VmaxVary(2,1)    = V2Vary;
VmaxMatrix(3,1)  = V3;
VmaxVary(3,1)    = V3Vary;
VmaxMatrix(4,1)  = V4;
VmaxVary(4,1)    = V4Vary;
VmaxMatrix(5,1)  = V5;
VmaxVary(5,1)    = V5Vary;
VmaxMatrix(6,1)  = V6;
VmaxVary(6,1)    = V6Vary;
VmaxMatrix(7,1)  = V7;
VmaxVary(7,1)    = V7Vary;
VmaxMatrix(8,1)  = V8;
VmaxVary(8,1)    = V8Vary;
VmaxMatrix(9,1)  = V9;
VmaxVary(9,1)    = V9Vary;
VmaxMatrix(10,1) = V10;
VmaxVary(10,1)   = V10Vary;
VmaxMatrix(11,1) = V11;
VmaxVary(11,1)   = V11Vary;
VmaxMatrix(12,1) = V12;
VmaxVary(12,1)   = V12Vary;
VmaxMatrix(13,1) = V13;
VmaxVary(13,1)   = V13Vary;
VmaxMatrix(14,1) = V16;
VmaxVary(14,1)   = V16Vary;
VmaxMatrix(15,1) = V21;
VmaxVary(15,1)   = V21Vary;
VmaxMatrix(16,1) = V22;
VmaxVary(16,1)   = V22Vary;
VmaxMatrix(17,1) = V23;
VmaxVary(17,1)   = V23Vary;
VmaxMatrix(18,1) = V31;
VmaxVary(18,1)   = V31Vary;
VmaxMatrix(19,1) = V32;
VmaxVary(19,1)   = V32Vary;
VmaxMatrix(20,1) = V33;
VmaxVary(20,1)   = V33Vary;
VmaxMatrix(21,1) = V111;
VmaxVary(21,1)   = V111Vary;
VmaxMatrix(22,1) = V112;
VmaxVary(22,1)   = V112Vary;
VmaxMatrix(23,1) = V113;
VmaxVary(23,1)   = V113Vary;
VmaxMatrix(24,1) = V121;
VmaxVary(24,1)   = V121Vary;
VmaxMatrix(25,1) = V122;
VmaxVary(25,1)   = V122Vary;
VmaxMatrix(26,1) = V123;
VmaxVary(26,1)   = V123Vary;
VmaxMatrix(27,1) = V124;
VmaxVary(27,1)   = V124Vary;
VmaxMatrix(28,1) = V131;
VmaxVary(28,1)   = V131Vary;
VmaxMatrix(29,1) = V1T;
VmaxVary(29,1)   = V1TVary;
VmaxMatrix(30,1) = V2T;
VmaxVary(30,1)   = V2TVary;
VmaxMatrix(31,1) = FIBF_PQT;
VmaxVary(31,1)   = FIBF_PQTVary;
VmaxMatrix(32,1) = ChlT;
VmaxVary(32,1)   = ChlTVary;
VmaxMatrix(33,1) = ChlT2;
VmaxVary(33,1)   = ChlT2Vary;
VmaxMatrix(34,1) = V2M;
VmaxVary(34,1)   = V2MVary; 
VmaxMatrix(35,1) = Vmax11;
VmaxVary(35,1)   = Vmax11Vary;
VmaxMatrix(36,1) = Tcyt;
VmaxVary(36,1)   = TcytVary;
VmaxMatrix(37,1) = BF_An;
VmaxVary(37,1)   = BF_AnVary;
VmaxMatrix(38,1) = BF_Cn;
VmaxVary(38,1)   = BF_CnVary;

VmaxMatrix	(	41	,	1	)	=	KM11	;
VmaxMatrix	(	42	,	1	)	=	KM12	;
VmaxMatrix	(	43	,	1	)	=	KM13	;
VmaxMatrix	(	44	,	1	)	=	KI11	;
VmaxMatrix	(	45	,	1	)	=	KI12	;
VmaxMatrix	(	46	,	1	)	=	KI13	;
VmaxMatrix	(	47	,	1	)	=	KI14	;
VmaxMatrix	(	48	,	1	)	=	KI15	;
VmaxMatrix	(	49	,	1	)	=	KM21	;
VmaxMatrix	(	50	,	1	)	=	KM22	;
VmaxMatrix	(	51	,	1	)	=	KM23	;
VmaxMatrix	(	52	,	1	)	=	KM31a	;
VmaxMatrix	(	53	,	1	)	=	KM32b	;
VmaxMatrix	(	54	,	1	)	=	KM41	;
VmaxMatrix	(	55	,	1	)	=	KM42	;
VmaxMatrix	(	56	,	1	)	=	KE4	;
VmaxMatrix	(	57	,	1	)	=	KM51	;
VmaxMatrix	(	58	,	1	)	=	KM52	;
VmaxMatrix	(	59	,	1	)	=	KM53	;
VmaxMatrix	(	60	,	1	)	=	KE5	;
VmaxMatrix	(	61	,	1	)	=	KM61	;
VmaxMatrix	(	62	,	1	)	=	KI61	;
VmaxMatrix	(	63	,	1	)	=	KI62	;
VmaxMatrix	(	64	,	1	)	=	KE6	;
VmaxMatrix	(	65	,	1	)	=	KM71	;
VmaxMatrix	(	66	,	1	)	=	KM72	;
VmaxMatrix	(	67	,	1	)	=	KM73	;
VmaxMatrix	(	68	,	1	)	=	KM74	;
VmaxMatrix	(	69	,	1	)	=	KE7	;
VmaxMatrix	(	70	,	1	)	=	KM8	;
VmaxMatrix	(	71	,	1	)	=	KM81	;
VmaxMatrix	(	72	,	1	)	=	KM82	;
VmaxMatrix	(	73	,	1	)	=	KE8	;
VmaxMatrix	(	74	,	1	)	=	KM9	;
VmaxMatrix	(	75	,	1	)	=	KI9	;
VmaxMatrix	(	76	,	1	)	=	KE9	;
VmaxMatrix	(	77	,	1	)	=	KM10	;
VmaxMatrix	(	78	,	1	)	=	KM101	;
VmaxMatrix	(	79	,	1	)	=	KM102	;
VmaxMatrix	(	80	,	1	)	=	KM103	;
VmaxMatrix	(	81	,	1	)	=	KE10	;
VmaxMatrix	(	82	,	1	)	=	KE11	;
VmaxMatrix	(	83	,	1	)	=	KE12	;
VmaxMatrix	(	84	,	1	)	=	KM131	;
VmaxMatrix	(	85	,	1	)	=	KM132	;
VmaxMatrix	(	86	,	1	)	=	KI131	;
VmaxMatrix	(	87	,	1	)	=	KI132	;
VmaxMatrix	(	88	,	1	)	=	KI133	;
VmaxMatrix	(	89	,	1	)	=	KI134	;
VmaxMatrix	(	90	,	1	)	=	KI135	;
VmaxMatrix	(	91	,	1	)	=	KE13	;
VmaxMatrix	(	92	,	1	)	=	KM161	;
VmaxMatrix	(	93	,	1	)	=	KM162	;
VmaxMatrix	(	94	,	1	)	=	KM163	;
VmaxMatrix	(	95	,	1	)	=	KE16	;
VmaxMatrix	(	96	,	1	)	=	KE21	;
VmaxMatrix	(	97	,	1	)	=	KE22	;
VmaxMatrix	(	98	,	1	)	=	KM231	;
VmaxMatrix	(	99	,	1	)	=	KM232	;
VmaxMatrix	(	100	,	1	)	=	KA231	;
VmaxMatrix	(	101	,	1	)	=	KA232	;
VmaxMatrix	(	102	,	1	)	=	KA233	;
VmaxMatrix	(	103	,	1	)	=	KI23	;
VmaxMatrix	(	104	,	1	)	=	KM311	;
VmaxMatrix	(	105	,	1	)	=	KM312	;
VmaxMatrix	(	106	,	1	)	=	KM313	;
VmaxMatrix	(	107	,	1	)	=	KM32	;
VmaxMatrix	(	108	,	1	)	=	KM33	;

								
VmaxMatrix	(	121	,	1	)	=	KO	;
VmaxMatrix	(	122	,	1	)	=	KC	;
VmaxMatrix	(	123	,	1	)	=	KR	;
VmaxMatrix	(	124	,	1	)	=	KM112	;
VmaxMatrix	(	125	,	1	)	=	KI1122	;
VmaxMatrix	(	126	,	1	)	=	KI1121	;
VmaxMatrix	(	127	,	1	)	=	KM1131	;
VmaxMatrix	(	128	,	1	)	=	KM1132	;
VmaxMatrix	(	129	,	1	)	=	KI113	;
VmaxMatrix	(	130	,	1	)	=	KE113	;
VmaxMatrix	(	131	,	1	)	=	KM121	;
VmaxMatrix	(	132	,	1	)	=	KM1221	;
VmaxMatrix	(	133	,	1	)	=	KM1222	;
VmaxMatrix	(	134	,	1	)	=	KI1221	;
VmaxMatrix	(	135	,	1   )	=	KE122	;
VmaxMatrix	(	136	,	1	)	=	KM123	;
VmaxMatrix	(	137	,	1	)	=	KI123	;
VmaxMatrix	(	138	,	1	)	=	KE123	;
VmaxMatrix	(	139	,	1	)	=	KM1241	;
VmaxMatrix	(	140	,	1	)	=	KM1242	;
VmaxMatrix	(	141	,	1	)	=	KI124	;
VmaxMatrix	(	142	,	1	)	=	KE124	;
VmaxMatrix	(	143	,	1	)	=	KM1311	;
VmaxMatrix	(	144	,	1	)	=	KI1311	;
VmaxMatrix	(	145	,	1	)	=	KI1312	;
VmaxMatrix	(	146	,	1	)	=	KM1312	;
VmaxMatrix	(	147	,	1	)	=	KM1011	;
VmaxMatrix	(	148	,	1	)	=	KI1011	;
VmaxMatrix	(	149	,	1	)	=	KM1012	;
VmaxMatrix	(	150	,	1	)	=	KI1012	;


VmaxMatrix	(	161	,1)=	KE501	;	%
VmaxMatrix	(	162	,1)=	Km511	;	%
VmaxMatrix	(	163	,1)=	Km512	;	%
VmaxMatrix	(	164	,1)=	Km513	;	%
VmaxMatrix	(	165	,1)=	KE51	;	%
VmaxMatrix	(	166	,1)=	Km514	;	%
VmaxMatrix	(	167	,1)=	Km521	;	%
VmaxMatrix	(	168	,1)=	KI521	;	%
VmaxMatrix	(	169	,1)=	KI522	;	%
VmaxMatrix	(	170	,1)=	KI523	;	%
VmaxMatrix	(	171	,1)=	KE52	;	%
VmaxMatrix	(	172	,1)=	KE531	;	%
VmaxMatrix	(	173	,1)=	KE541	;	%
VmaxMatrix	(	174	,1)=	Km551	;	%
VmaxMatrix	(	175	,1)=	Km552	;	%
VmaxMatrix	(	176	,1)=	Km553	;	%
VmaxMatrix	(	177	,1)=	Km554	;	%
VmaxMatrix	(	178	,1)=	KE55	;	%
VmaxMatrix	(	179	,1)=	Km561	;	%
VmaxMatrix	(	180	,1)=	Km562	;	%
VmaxMatrix	(	181	,1)=	KI561	;	%
VmaxMatrix	(	182	,1)=	KI562	;	%
VmaxMatrix	(	183	,1)=	KI563	;	%
VmaxMatrix	(	184	,1)=	KI564	;	%
VmaxMatrix	(	185	,1)=	KI565	;	%
VmaxMatrix	(	186	,1)=	KE56	;	%
VmaxMatrix	(	187	,1)=	Km571	;	%
VmaxMatrix	(	188	,1)=	Ki572	;	%
VmaxMatrix	(	189	,1)=	KE57	;	%
VmaxMatrix	(	190	,1)=	Km581	;	%
VmaxMatrix	(	191	,1)=	KI581	;	%
VmaxMatrix	(	192	,1)=	KI582	;	%
VmaxMatrix	(	193	,1)=	Km591	;	%
VmaxMatrix	(	194	,1)=	Km592	;	%
VmaxMatrix	(	195	,1)=	Km593	;	%
VmaxMatrix	(	196	,1)=	KI591	;	%
VmaxMatrix	(	197	,1)=	KI592	;	%
VmaxMatrix	(	198	,1)=	KE59	;	%
VmaxMatrix	(	199	,1)=	Km601	;	%
VmaxMatrix	(	200	,1)=	Km602	;	%
VmaxMatrix	(	201	,1)=	Km603	;	%
VmaxMatrix	(	202	,1)=	Km604	;	%
VmaxMatrix	(	203	,1)=	KE60	;	%
VmaxMatrix	(	204	,1)=	KE61	;	%
VmaxMatrix	(	205	,1)=	Km621	;	%
VmaxMatrix	(	206	,1)=	V51	;	%
VmaxMatrix	(	207	,1)=	V52	;	%
VmaxMatrix	(	208	,1)=	V55	;	%
VmaxMatrix	(	209	,1)=	V56	;	%
VmaxMatrix	(	210	,1)=	V57	;	%
VmaxMatrix	(	211	,1)=	V58	;	%
VmaxMatrix	(	212	,1)=	V59	;	%
VmaxMatrix	(	213	,1)=	V60	;	%
VmaxMatrix	(	214	,1)=	V61	;	%
VmaxMatrix	(	215	,1)=	V62	;	%
VmaxMatrix	(	216	,1)=	Vdhap_in	;	%
VmaxMatrix	(	217	,1)=	Vgap_in 	;	%
VmaxMatrix	(	218	,1)=	Vpga_in  	;	%





VmaxVary	(	41	,	1	)	=	KM11_vary	;
VmaxVary	(	42	,	1	)	=	KM12_vary	;
VmaxVary	(	43	,	1	)	=	KM13_vary	;
VmaxVary	(	44	,	1	)	=	KI11_vary	;
VmaxVary	(	45	,	1	)	=	KI12_vary	;
VmaxVary	(	46	,	1	)	=	KI13_vary	;
VmaxVary	(	47	,	1	)	=	KI14_vary	;
VmaxVary	(	48	,	1	)	=	KI15_vary	;
VmaxVary	(	49	,	1	)	=	KM21_vary	;
VmaxVary	(	50	,	1	)	=	KM22_vary	;
VmaxVary	(	51	,	1	)	=	KM23_vary	;
VmaxVary	(	52	,	1	)	=	KM31a_vary	;
VmaxVary	(	53	,	1	)	=	KM32b_vary	;
VmaxVary	(	54	,	1	)	=	KM41_vary	;
VmaxVary	(	55	,	1	)	=	KM42_vary	;
VmaxVary	(	56	,	1	)	=	KE4_vary	;
VmaxVary	(	57	,	1	)	=	KM51_vary	;
VmaxVary	(	58	,	1	)	=	KM52_vary	;
VmaxVary	(	59	,	1	)	=	KM53_vary	;
VmaxVary	(	60	,	1	)	=	KE5_vary	;
VmaxVary	(	61	,	1	)	=	KM61_vary	;
VmaxVary	(	62	,	1	)	=	KI61_vary	;
VmaxVary	(	63	,	1	)	=	KI62_vary	;
VmaxVary	(	64	,	1	)	=	KE6_vary	;
VmaxVary	(	65	,	1	)	=	KM71_vary	;
VmaxVary	(	66	,	1	)	=	KM72_vary	;
VmaxVary	(	67	,	1	)	=	KM73_vary	;
VmaxVary	(	68	,	1	)	=	KM74_vary	;
VmaxVary	(	69	,	1	)	=	KE7_vary	;
VmaxVary	(	70	,	1	)	=	KM8_vary	;
VmaxVary	(	71	,	1	)	=	KM81_vary	;
VmaxVary	(	72	,	1	)	=	KM82_vary	;
VmaxVary	(	73	,	1	)	=	KE8_vary	;
VmaxVary	(	74	,	1	)	=	KM9_vary	;
VmaxVary	(	75	,	1	)	=	KI9_vary	;
VmaxVary	(	76	,	1	)	=	KE9_vary	;
VmaxVary	(	77	,	1	)	=	KM10_vary	;
VmaxVary	(	78	,	1	)	=	KM101_vary	;
VmaxVary	(	79	,	1	)	=	KM102_vary	;
VmaxVary	(	80	,	1	)	=	KM103_vary	;
VmaxVary	(	81	,	1	)	=	KE10_vary	;
VmaxVary	(	82	,	1	)	=	KE11_vary	;
VmaxVary	(	83	,	1	)	=	KE12_vary	;
VmaxVary	(	84	,	1	)	=	KM131_vary	;
VmaxVary	(	85	,	1	)	=	KM132_vary	;
VmaxVary	(	86	,	1	)	=	KI131_vary	;
VmaxVary	(	87	,	1	)	=	KI132_vary	;
VmaxVary	(	88	,	1	)	=	KI133_vary	;
VmaxVary	(	89	,	1	)	=	KI134_vary	;
VmaxVary	(	90	,	1	)	=	KI135_vary	;
VmaxVary	(	91	,	1	)	=	KE13_vary	;
VmaxVary	(	92	,	1	)	=	KM161_vary	;
VmaxVary	(	93	,	1	)	=	KM162_vary	;
VmaxVary	(	94	,	1	)	=	KM163_vary	;
VmaxVary	(	95	,	1	)	=	KE16_vary	;
VmaxVary	(	96	,	1	)	=	KE21_vary	;
VmaxVary	(	97	,	1	)	=	KE22_vary	;
VmaxVary	(	98	,	1	)	=	KM231_vary	;
VmaxVary	(	99	,	1	)	=	KM232_vary	;
VmaxVary	(	100	,	1	)	=	KA231_vary	;
VmaxVary	(	101	,	1	)	=	KA232_vary	;
VmaxVary	(	102	,	1	)	=	KA233_vary	;
VmaxVary	(	103	,	1	)	=	KI23_vary	;
VmaxVary	(	104	,	1	)	=	KM311_vary	;
VmaxVary	(	105	,	1	)	=	KM312_vary	;
VmaxVary	(	106	,	1	)	=	KM313_vary	;
VmaxVary	(	107	,	1	)	=	KM32_vary	;
VmaxVary	(	108	,	1	)	=	KM33_vary	;

VmaxVary	(	121	,	1	)	=	KO_vary	;
VmaxVary	(	122	,	1	)	=	KC_vary	;
VmaxVary	(	123	,	1	)	=	KR_vary	;
VmaxVary	(	124	,	1	)	=	KM112_vary	;
VmaxVary	(	125	,	1	)	=	KI1122_vary	;
VmaxVary	(	126	,	1	)	=	KI1121_vary	;
VmaxVary	(	127	,	1	)	=	KM1131_vary	;
VmaxVary	(	128	,	1	)	=	KM1132_vary	;
VmaxVary	(	129	,	1	)	=	KI113_vary	;
VmaxVary	(	130	,	1	)	=	KE113_vary	;
VmaxVary	(	131	,	1	)	=	KM121_vary	;
VmaxVary	(	132	,	1	)	=	KM1221_vary	;
VmaxVary	(	133	,	1	)	=	KM1222_vary	;
VmaxVary	(	134	,	1	)	=	KI1221_vary	;
VmaxVary	(	135	,	1	)	=	KE122_vary	;
VmaxVary	(	136	,	1	)	=	KM123_vary	;
VmaxVary	(	137	,	1	)	=	KI123_vary	;
VmaxVary	(	138	,	1	)	=	KE123_vary	;
VmaxVary	(	139	,	1	)	=	KM1241_vary	;
VmaxVary	(	140	,	1	)	=	KM1242_vary	;
VmaxVary	(	141	,	1	)	=	KI124_vary	;
VmaxVary	(	142	,	1	)	=	KE124_vary	;
VmaxVary	(	143	,	1	)	=	KM1311_vary	;
VmaxVary	(	144	,	1	)	=	KI1311_vary	;
VmaxVary	(	145	,	1	)	=	KI1312_vary	;
VmaxVary	(	146	,	1	)	=	KM1312_vary	;
VmaxVary	(	147	,	1	)	=	KM1011_vary	;
VmaxVary	(	148	,	1	)	=	KI1011_vary	;
VmaxVary	(	149	,	1	)	=	KM1012_vary	;
VmaxVary	(	150	,	1	)	=	KI1012_vary	;

VmaxVary	(161	,1)=	KE501_vary	;
VmaxVary	(162	,1)=	Km511_vary	;
VmaxVary	(163	,1)=	Km512_vary	;
VmaxVary	(164	,1)=	Km513_vary	;
VmaxVary	(165	,1)=	KE51_vary	;
VmaxVary	(166	,1)=	Km514_vary	;
VmaxVary	(167	,1)=	Km521_vary	;
VmaxVary	(168	,1)=	KI521_vary	;
VmaxVary	(169	,1)=	KI522_vary	;
VmaxVary	(170	,1)=	KI523_vary	;
VmaxVary	(171	,1)=	KE52_vary	;
VmaxVary	(172	,1)=	KE531_vary	;
VmaxVary	(173	,1)=	KE541_vary	;
VmaxVary	(174	,1)=	Km551_vary	;
VmaxVary	(175	,1)=	Km552_vary	;
VmaxVary	(176	,1)=	Km553_vary	;
VmaxVary	(177	,1)=	Km554_vary	;
VmaxVary	(178	,1)=	KE55_vary	;
VmaxVary	(179	,1)=	Km561_vary	;
VmaxVary	(180	,1)=	Km562_vary	;
VmaxVary	(181	,1)=	KI561_vary	;
VmaxVary	(182	,1)=	KI562_vary	;
VmaxVary	(183	,1)=	KI563_vary	;
VmaxVary	(184	,1)=	KI564_vary	;
VmaxVary	(185	,1)=	KI565_vary	;
VmaxVary	(186	,1)=	KE56_vary	;
VmaxVary	(187	,1)=	Km571_vary	;
VmaxVary	(188	,1)=	Ki572_vary	;
VmaxVary	(189	,1)=	KE57_vary	;
VmaxVary	(190	,1)=	Km581_vary	;
VmaxVary	(191	,1)=	KI581_vary	;
VmaxVary	(192	,1)=	KI582_vary	;
VmaxVary	(193	,1)=	Km591_vary	;
VmaxVary	(194	,1)=	Km592_vary	;
VmaxVary	(195	,1)=	Km593_vary	;
VmaxVary	(196	,1)=	KI591_vary	;
VmaxVary	(197	,1)=	KI592_vary	;
VmaxVary	(198	,1)=	KE59_vary	;
VmaxVary	(199	,1)=	Km601_vary	;
VmaxVary	(200	,1)=	Km602_vary	;
VmaxVary	(201	,1)=	Km603_vary	;
VmaxVary	(202	,1)=	Km604_vary	;
VmaxVary	(203	,1)=	KE60_vary	;
VmaxVary	(204	,1)=	KE61_vary	;
VmaxVary	(205	,1)=	Km621_vary	;
VmaxVary	(206	,1)=	V51_vary	;
VmaxVary	(207	,1)=	V52_vary	;
VmaxVary	(208	,1)=	V55_vary	;
VmaxVary	(209	,1)=	V56_vary	;
VmaxVary	(210	,1)=	V57_vary	;
VmaxVary	(211	,1)=	V58_vary	;
VmaxVary	(212	,1)=	V59_vary	;
VmaxVary	(213	,1)=	V60_vary	;
VmaxVary	(214	,1)=	V61_vary	;
VmaxVary	(215	,1)=	V62_vary	;
VmaxVary	(216	,1)=	Vdhap_in_vary	;
VmaxVary	(217	,1)=	Vgap_in_vary	;
VmaxVary	(218	,1)=	Vpga_in_vary	;

VmaxVary = zeros(218,1);
VmaxVary	(1:218	,	1	)	=	1	;
%VmaxVary(1,1) = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is where the user chooses which metabolite to calculate elasticity 
% coefficients against.  Also, one can choose 'the rate of carbon
% consumption' as the dependent variable to which one calculates flux control
% coefficients. If you would like to calculate for one metabolites, assign
% the value for that metabolites to be 1. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MetConCoeff2(1)  = 0;         %   RuBP
MetConCoeff2(2)  = 0;         %   PGA
MetConCoeff2(3)  = 0;         %   DPGA
MetConCoeff2(4)  = 0;         %   T3P
MetConCoeff2(5)  = 0;         %   FBP
MetConCoeff2(6)  = 0;         %   E4P
MetConCoeff2(7)  = 0;         %   S7P
MetConCoeff2(8)  = 0;         %   SBP
MetConCoeff2(9)  = 0;         %   ATP
MetConCoeff2(10) = 0;         %   NADPH
MetConCoeff2(11) = 0;         %   CO2
MetConCoeff2(12) = 0;         %   O2
MetConCoeff2(13) = 0;         %   HexP
MetConCoeff2(14) = 0;         %   PenP
MetConCoeff2(15) = 0;         %   GCEA  
MetConCoeff2(16) = 0;         %   GCA
MetConCoeff2(17) = 0;         %   PGCA
MetConCoeff2(18) = 0;         %   GCAc
MetConCoeff2(19) = 0;         %   GOAc
MetConCoeff2(20) = 0;         %   SERc
MetConCoeff2(21) = 0;         %   GLYc
MetConCoeff2(22) = 0;         %   HPRc
MetConCoeff2(23) = 0;         %   GCEAc
MetConCoeff2(24) = 1;         %   Carbon Rate

MetConCoeff = MetConCoeff2';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global RedoxReg_RA_com;
RedoxReg_RA_com =0;


global RuACT_EPS_com;
RuACT_EPS_com =1;

global GP ;
GP = 0;
