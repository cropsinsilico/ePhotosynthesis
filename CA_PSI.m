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



function PSs = CA_PSI(Begin);

global PS_C_CA;             %   Global constant for the total adenylates
global PS_C_CP;             %   Global constant for the total phosphate
global PS_C_CN;             %   Global constant for the total NADP+NADPH
global PS_PEXT;             %   Global constant for the cytosolic Phosphate concentration;

PS_C_CP	= 15;   % 15
PS_C_CA	=1.5;           % 1.5 is the default
PS_C_CN	=1;
PS_PEXT	=0.5;    

global PSPR_RA_CA;
PSPR_RA_CA = PS_C_CA;


RuBP	= 2.000;
PGA	    =2.400;
DPGA	=0.0011;
T3P	    =0.5;
NONE	=0;
FBP	    =0.670;
E4P	    =0.050;
S7P	    =2.000;
SBP	    =0.300;
ATP	    =0.68;         % 0.68 is the default for the dark adapted leaves when the total ATP and ADP is 1.5. 
NADPH	=0.21;         % Set this value to be 0.5
HexP    = 2.2;
PenP    = 0.25;

CO2	    = 0.012;        
O2	    = 0.264;        

PSs(1)	=RuBP;
PSs(2)	=PGA;
PSs(3)	=DPGA;
PSs(4)	=T3P;
PSs(5)	=NONE;
PSs(6)	=FBP;
PSs(7)	=E4P;
PSs(8)	=S7P;
PSs(9)	=SBP;
PSs(10)	=ATP;
PSs(11)	=NADPH;
PSs(12)	= CO2;
PSs(13)	= O2;
PSs(14) = HexP;
PSs(15) = PenP;



		% Initialize the constants for the different reactions
		
		global	KM11	;
		global	KM12	;
		global	KM13	;
		global  KI11    ;
		global  KI12    ;
		global  KI13    ;
		global  KI14    ;
		global  KI15    ;
		
		global	KM21	;
		global	KM22	;
		global  KM23    ;
		
		global	KM31a	;
		global	KM32b	;
		
		global	KM41	;
		global	KM42	;
		global  KE4     ;
		
		global	KM51	;
		global	KM52	;
		global	KM53	;
		global  KE5     ;
		
		global	KM61	;
		global  KI61    ;
		global  KI62    ;
		
		global	KM71	;
		global  KM72    ;
		global	KM73	;
		global  KM74    ;
		
		global	KM8	    ;
		global  KM81    ;
		global  KM82    ;
		
		global	KM9	    ;
		global  KI9     ;
		
		global	KM10	;
		global  KM101   ;
		global  KM102   ;
		global  KM103   ;
		
		global	KE11	;
		
		global	KE12	;
		
		global	KM131	;
		global	KM132	;
		global	KI131	;
		global	KI132	;
		global	KI133	;
		global	KI134	;
		global	KI135	;
		
		global	KM161	;
		global	KM162	;
		
		global	KE21	;
		
		global	KE22	;
		
		global	KM231	;
		global	KM232	;
		global	KA231	;
		global	KA232	;
		global	KA233	;
		global	KI23	;
		
		global	KM311	;
		global	KM312	;
		global	KM313	;
		
		global	KM32	;
		global	KM33	;
		
		
		global KE6;
		global KE7;
		global KE8;
		global KE9;
		global KE10;
		global KE13;
		global KE16;
		
		global KM103;
		global KM163;
		
		
		KM11	=	0.014;		% 	CO2	1	RuBP+CO2->2PGA
		KM12	=	0.28;		%	O2	1	RuBP+CO2->2PGA
		KM13	=	0.02;		% 	RuBP	1	RuBP+CO2->2PGA
		
		KI11    =   0.84   ;    % PGA
		KI12    =0.04   ;       % FBP
		KI13    = 0.075 ;       % SBP
		KI14    = 0.9   ;       % Pi
		KI15    = 0.07  ;       % NADPH
		
		
		KM21	=	0.240;		%	PGA	2	PGA+ATP <-> ADP + DPGA
		KM22	=	0.390;		% 	ATP	2	PGA+ATP <-> ADP + DPGA
		KM23    =   0.23  ;        %  ADP     
		
		KM31a	=	0.004;		%	BPGA	3	DPGA+NADPH <->GAP + OP+NADP 
		KM32b	=	0.1	;	    % 	NADPH	3	DPGA+NADPH <->GAP + OP+NADP
		
		KM41	=	2.5	;	    %	DHAP	4	DHAP <->GAP
		KM42	=	0.68;		% 	GAP	4	DHAP <->GAP
		KE4     =   0.05;       %   Using the value from Patterson
		
		KM51	=	0.3	;	    %	GAP	5	GAP+DHAP <->FBP
		KM52	=	0.4	;	    % 	DHAP	5	GAP+DHAP <->FBP
		KM53	=	0.02;		%	FBP	5	GAP+DHAP <->FBP     % Original Value: 0.02
		KE5     = 7.100;          % Defult: 7.1
		
		KM61	=	0.033;		% 	FBP	6	FBP<->F6P+OP
		KI61    = 0.7   ;       %   F6P       
		KI62    = 12    ;       %   Pi
		KE6     =   6.66 * 10^5;    % The equilibrium constant for this reaction        % New    mM     Laisk or Bassham and Krause 1969 BBA
		
		KM71	=	0.100;		%	Xu5P	7	F6P+GAP<->E4P+Xu5P      % jn
		KM72	=	0.100;		% 	E4P	7	F6P+GAP<->E4P+Xu5P
		KM73    = 0.1;         %   F6P This value was based on estimate
		KM74    = 0.1000;         % Estimate for GAP ORIGINAL 0.1
		KE7     =   10 ;       % The equilibrium constant for this reaction             % New           Laisk  Bassham and Krause 1969 BBA
		
		KM8	    =	0.02;		%	SBP	8	E4P+DHAP<->SBP
		KM81    = 0.4   ;       % DHAP
		KM82    = 0.2   ;       % E4P estimate
		KE8     = 1.017 ;     % The equilibrium constant for this reaction                  % New    mM-1         Laisk  Bassham and Krause 1969 BBA. Default: 1.107
		
		KM9	    =	0.05;		% 	SBP	9	SBP<->S7P+OP    
		KI9     = 12    ;       %   The inibintion constant for Pi; 
		KE9     =   6.66 * 10^5 ; % The equilibrium constant of this reaction           % New   mM      Laisk  Bassham and Krause 1969 BBA
		
		
		KM10	=	1.5	;	    %	R5P	10	S7P+GAP<->Ri5P+Xu5P
		KM101   =   0.1 ;       %   Xu5P
		KM102   = 0.072 ;       %   Estimate for GAP
		KM103   = 0.46 ;        %   Estimate for S7P                                    % New 
		KE10    = 1/0.85 ;      %   The equilibrium constant for this reaction          % New From Laisk or Bassham and Krause 1969 BBA
		
		KE11	=	0.4	;	    %	Equilibrium Constant	11	Ri5P<-->Ru5P
		KE12	=	0.67;		% 	Equilibrium Constant	12	Xu5P<-->Ru5P
		
		KM131	=	0.05;		    %	Ru5P	13	Ru5P+ATP<->RuBP+ADP
		KM132	=	0.059;		    % 	ATP	13	Ru5P+ATP<->RuBP+ADP
		KI131	=	2	;			%	PGA	13	Ru5P+ATP<->RuBP+ADP
		KI132	=	0.7	;			%	RuBP	13	Ru5P+ATP<->RuBP+ADP
		KI133	=	4	;			%	Pi	13	Ru5P+ATP<->RuBP+ADP
		KI134	=	2.5	;			%	ADP	13	Ru5P+ATP<->RuBP+ADP
		KI135	=	0.4	;			%	ADP	13	Ru5P+ATP<->RuBP+ADP
		KE13    =   6.846 * 10^3;   %   The equilibrium constant for this reaction  % New From Laisk or Bassham and Krause 1969 BBA
		
		KM161	=	0.014;		%	ADP	16	ADP+Pi<->ATP 
		KM162	=	0.3;		% 	Pi	16	ADP+Pi<-> ATP
		KM163   =   0.3;        %   ATP 16  ADP+Pi<-> ATP                           % New       Based on Laisk  
		KE16    =   16;      %   The equilibrium constant for this reaction      % NEW, From Laisk or Bassham and Krause 1969 BBA
		
		
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
		


% Initialize the Vmax for different reactions

global	V1	;
global	V2	;
global	V3	;
% global	V4	;
global	V5	;
global	V6	;
global	V7	;
global	V8	;
global	V9	;
global	V10	;
global	V11	;
global	V12	;
global	V13	;
global	V16	;
global	V21	;
global	V22	;
global	V23	;
global	V31	;
global	V32	;
global	V33	;

% FC is a fussl factor here. 
FC = 1;       % Defulat is 2.5. 
fc16 = 1;     % default 1.5.
% Initialize the values of the global variables
SC = 1;        % Scalling coefficient for the stroma volume per mg chl. defualt 2

global J; 
J =1;
global VmaxXValue;

V1		=	3.38 *SC	;	%	(Harris & Koniger, 1997)	1	Rubisco	RuBP+CO2<->2PGA         3.73
V2		=	73.4*SC	;	%	(Harris & Koniger, 1997)	2	PGA Kinase	PGA+ATP <-> ADP + DPGA; Initial as 11.76
V3		=	9.8*SC	    ;	%	(Harris & Koniger, 1997)	3	GAP dehydragenase	DPGA+NADPH <->GAP + OP+NADP 
V4		=	0	        ;	%	(Harris & Koniger, 1997)	4	Triose phosphate isomerase	DHAP <->GAP Original 40.53
V5		=	2.96*SC	    ;	%	(Harris & Koniger, 1997)	5	Aldolase	GAP+DHAP <->FBP 5.52
V6		=	1.76*SC	;	%	(Harris & Koniger, 1997)	6	FBPase	FBP<->F6P+OP    1.155
V7		=	7.6*SC	;	%	(Harris & Koniger, 1997)	7	Transketolase	F6P+GAP<->E4P+Xu5P  Original vlaue 1.008
V8		=	2.96*SC	    ;	%	(Harris & Koniger, 1997)	8	Aldolase	E4P+DHAP<->SBP 5.52
V9		=	0.79*SC *FC	;	%	(Harris & Koniger, 1997)	9	SBPase	SBP<->S7P+OP    0.168 as original value; 0.4168 was its value.
V10		=	7.605*SC	    ;	%	(Harris & Koniger, 1997)	10	Transketolase	S7P+GAP<->Ri5P+Xu5P 1.008
V11		=	0	        ;	%	(Harris & Koniger, 1997)	11	Pentosephosphate isomerase	Ri5P<-->Ru5P
V12		=	0	        ;	%	(Harris & Koniger, 1997)	12	Pentosephosphate epimerase	Xu5P<-->Ru5P
V13		=	26.39*SC	    ;	%	(Harris & Koniger, 1997)	13	Ribulosebiphosphate kinase	Ru5P+ATP<->RuBP+ADP
V16		=	15   ;  %;5.47 * SC	* fc16;	%	(Aflalo & Shavit, 1983, Davenport & McLeod, 1986)	16	ATP synthase	ADP+Pi<->ATP    1.47
V21		=	0	        ;	%	                            21	Hexose phosphate isomerase	F6P<->G6P
V22		=	0	        ;	%		                        22	Phosphoglucomutase	G6P<->G1P
V23		=	0.65;   % 0.5 * SC	* FC * 1.1;	%	(Latzko, Steup & Schachtele, 1981)	23	ADP-glucose pyrophosphorylase and	ADPG+Gn<->G(n+1)+ADP 0.18
V31 = 1;
V32=1;
V33=1;


SC = 1;

ControlAnalysis=1;
    
if ControlAnalysis ==1

V1      = VmaxXValue(1,J);  
V2		= VmaxXValue(2,J);  
V3		= VmaxXValue(3,J);  
V4		= VmaxXValue(4,J);  
V5		= VmaxXValue(5,J);  
V6		= VmaxXValue(6,J);    
V7		= VmaxXValue(7,J); 

V8		= VmaxXValue(8,J);  
V9		= VmaxXValue(9,J);  
V10		= VmaxXValue(10,J); 
V11		= VmaxXValue(11,J); 

V12		= VmaxXValue(12,J); 
V13		= VmaxXValue(13,J); 
V16		= VmaxXValue(14,J);                          	
V21		= VmaxXValue(15,J); 
V22		= VmaxXValue(16,J); 
V23		= VmaxXValue(17,J); 
V31		= VmaxXValue(18,J); 
V32		= VmaxXValue(19,J); 
V33		= VmaxXValue(20,J); 


KM11	=	VmaxXValue	(	41	,	J	)	;
KM12	=	VmaxXValue	(	42	,	J	)	;
KM13	=	VmaxXValue	(	43	,	J	)	;
KI11	=	VmaxXValue	(	44	,	J	)	;
KI12	=	VmaxXValue	(	45	,	J	)	;
KI13	=	VmaxXValue	(	46	,	J	)	;
KI14	=	VmaxXValue	(	47	,	J	)	;
KI15	=	VmaxXValue	(	48	,	J	)	;
KM21	=	VmaxXValue	(	49	,	J	)	;
KM22	=	VmaxXValue	(	50	,	J	)	;
KM23	=	VmaxXValue	(	51	,	J	)	;
KM31a	=	VmaxXValue	(	52	,	J	)	;
KM32b	=	VmaxXValue	(	53	,	J	)	;
KM41	=	VmaxXValue	(	54	,	J	)	;
KM42	=	VmaxXValue	(	55	,	J	)	;
KE4	=	VmaxXValue	(	56	,	J	)	;
KM51	=	VmaxXValue	(	57	,	J	)	;
KM52	=	VmaxXValue	(	58	,	J	)	;
KM53	=	VmaxXValue	(	59	,	J	)	;
KE5	=	VmaxXValue	(	60	,	J	)	;
KM61	=	VmaxXValue	(	61	,	J	)	;
KI61	=	VmaxXValue	(	62	,	J	)	;
KI62	=	VmaxXValue	(	63	,	J	)	;
KE6	=	VmaxXValue	(	64	,	J	)	;
KM71	=	VmaxXValue	(	65	,	J	)	;
KM72	=	VmaxXValue	(	66	,	J	)	;
KM73	=	VmaxXValue	(	67	,	J	)	;
KM74	=	VmaxXValue	(	68	,	J	)	;
KE7	=	VmaxXValue	(	69	,	J	)	;
KM8	=	VmaxXValue	(	70	,	J	)	;
KM81	=	VmaxXValue	(	71	,	J	)	;
KM82	=	VmaxXValue	(	72	,	J	)	;
KE8	=	VmaxXValue	(	73	,	J	)	;
KM9	=	VmaxXValue	(	74	,	J	)	;
KI9	=	VmaxXValue	(	75	,	J	)	;
KE9	=	VmaxXValue	(	76	,	J	)	;
KM10	=	VmaxXValue	(	77	,	J	)	;
KM101	=	VmaxXValue	(	78	,	J	)	;
KM102	=	VmaxXValue	(	79	,	J	)	;
KM103	=	VmaxXValue	(	80	,	J	)	;
KE10	=	VmaxXValue	(	81	,	J	)	;
KE11	=	VmaxXValue	(	82	,	J	)	;
KE12	=	VmaxXValue	(	83	,	J	)	;
KM131	=	VmaxXValue	(	84	,	J	)	;
KM132	=	VmaxXValue	(	85	,	J	)	;
KI131	=	VmaxXValue	(	86	,	J	)	;
KI132	=	VmaxXValue	(	87	,	J	)	;
KI133	=	VmaxXValue	(	88	,	J	)	;
KI134	=	VmaxXValue	(	89	,	J	)	;
KI135	=	VmaxXValue	(	90	,	J	)	;
KE13	=	VmaxXValue	(	91	,	J	)	;
KM161	=	VmaxXValue	(	92	,	J	)	;
KM162	=	VmaxXValue	(	93	,	J	)	;
KM163	=	VmaxXValue	(	94	,	J	)	;
KE16	=	VmaxXValue	(	95	,	J	)	;
KE21	=	VmaxXValue	(	96	,	J	)	;
KE22	=	VmaxXValue	(	97	,	J	)	;
KM231	=	VmaxXValue	(	98	,	J	)	;
KM232	=	VmaxXValue	(	99	,	J	)	;
KA231	=	VmaxXValue	(	100	,	J	)	;
KA232	=	VmaxXValue	(	101	,	J	)	;
KA233	=	VmaxXValue	(	102	,	J	)	;
KI23	=	VmaxXValue	(	103	,	J	)	;
KM311	=	VmaxXValue	(	104	,	J	)	;
KM312	=	VmaxXValue	(	105	,	J	)	;
KM313	=	VmaxXValue	(	106	,	J	)	;
KM32	=	VmaxXValue	(	107	,	J	)	;
KM33	=	VmaxXValue	(	108	,	J	)	;


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%H%%%%%%%%%
%   Here is the location for transfering variables  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global PS2RA_RuBP_ini;
PS2RA_RuBP_ini = RuBP;

global BF_FI_com;
global PS2BF_ATP;
global PS2BF_ADP;
global PS2BF_Pi;

PS2BF_ATP = ATP;
PS2BF_ADP = PS_C_CA - ATP;