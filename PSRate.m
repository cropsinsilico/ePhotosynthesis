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



function PSr = PSRate(t,PSs,Param)


global PS_C_CA;             %   Global constant for the total adenylates
global PS_C_CP;             %   Global constant for the total phosphate
global PS_C_CN;             %   Global constant for the total NADP+NADPH
global PS_PEXT;             %   Global constant for the cytosolic Phosphate concentration;

% For output of the flux

global PS_VEL;
global PS_TIME_N;
global PS_OLD_TIME;

PsCA = PS_C_CA;
PsCP = PS_C_CP;
PsCN = PS_C_CN;
PsPEXT = PS_PEXT;

% First the physical and chemical constant for all the reactions

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
global	KM72	;
global  KM73    ;
global  KM74    ;
global	KM8	;
global  KM81    ;
global  KM82    ;
global	KM9	;
global  KI9;
global	KM10	;
global  KM101   ;
global  KM102   ;
global  KM103   ;
global	KM11	;
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
%global	KM231	;
%global	KM232	;
%global	KA231	;
%global	KA232	;
%global	KA233	;
%global	KI23	;
global	KM311	;
global	KM312	;
global	KM313	;
global	KM32	;
global	KM33	;

global	KM231	;
global	KM232	;
global	KM233	;
global	KM234	;
global	KA231	;
global	KI231	;
global  KVmo    ;
global  KE23    ;

global KM241;
global KM242;
global KE24;

global KE25;    

global KE6;
global KE7;
global KE8;
global KE9;
global KE10;
global KE13;
global KE16;

global KM103;
global KM163;


% Set the value to the local variables: for example: PrKM11

PsKM11	=	KM11	;	% 	CO2	1	RuBP+CO2->2PGA
PsKM12	=	KM12	;	%	O2	1	RuBP+CO2->2PGA
PsKM13	=	KM13	;	% 	RuBP	1	RuBP+CO2->2PGA
PsKI11  =   KI11    ;   %   PGA
PsKI12  =   KI12    ;   %   FBP
PsKI13  =   KI13    ;   %   SBP
PsKI14  =   KI14    ;   %   Pi
PsKI15  =   KI15    ;   %   NADPH
PsKM21	=	KM21	;	%	PGA	2	PGA+ATP <-> ADP + DPGA
PsKM22	=	KM22	;	% 	ATP	2	PGA+ATP <-> ADP + DPGA
PsKM23  =   KM23    ;   %   ADP
PsKM31a	=	KM31a	;	%	BPGA	3	DPGA+NADPH <->GAP + OP+NADP 
PsKM32b	=	KM32b	;	% 	NADPH	3	DPGA+NADPH <->GAP + OP+NADP
PsKM41	=	KM41	;	%	DHAP	4	DHAP <->GAP
PsKM42	=	KM42	;	% 	GAP	4	DHAP <->GAP
PsKE4   =   KE4     ;
PsKM51	=	KM51	;	%	GAP	5	GAP+DHAP <->FBP
PsKM52	=	KM52	;	% 	DHAP	5	GAP+DHAP <->FBP
PsKM53	=	KM53	;	%	FBP	5	GAP+DHAP <->FBP
PsKE5   =   KE5     ;
PsKM61	=	KM61	;	% 	FBP	6	FBP<->F6P+OP
PsKI61  =   KI61    ;
PsKI62  =   KI62    ;
PsKM71	=	KM71	;	%	Xu5P	7	F6P+GAP<->E4P+Xu5P
PsKM72	=	KM72	;	% 	E4P	7	F6P+GAP<->E4P+Xu5P
PsKM73  =   KM73    ;   %   Estimate for F6P
PsKM74  =   KM74    ;   %   Estimate for GAP

PsKM8	=	KM8	    ;	%	SBP	8	E4P+DHAP<->SBP
PsKM81   =  KM81    ;   % DHAP
PsKM82  =   KM82    ;   % E4P
PsKM9	=	KM9	    ;	% 	SBP	9	SBP<->S7P+OP
PsKI9   =   KI9     ;
PsKM10	=	KM10	;	%	R5P	10	S7P+GAP<->Ri5P+Xu5P
PsKM101 =   KM101   ;   %   Xu5P
PsKM102 =   KM102   ;   %   GAP estimate
PsKM103 =   KM103   ;   %   S7P estimate
PsKE11	=	KE11	;	%	Equilibrium Constant	11	Ri5P<-->Ru5P
PsKE12	=	KE12	;	% 	Equilibrium Constant	12	Xu5P<-->Ru5P
PsKM131	=	KM131	;	%	Ru5P	13	Ru5P+ATP<->RuBP+ADP
PsKM132	=	KM132	;	% 	ATP	13	Ru5P+ATP<->RuBP+ADP
PsKI131	=	KI131	;	%	PGA
PsKI132	=	KI132	;	%	RuBP
PsKI133	=	KI133	;	%	Pi
PsKI134	=	KI134	;	%	ADP
PsKI135	=	KI135	;	%	ADP
PsKM161	=	KM161	;	%	ADP	16	ADP+Pi<->ATP
PsKM162	=	KM162	;	% 	Pi	16	ADP+Pi<-> ATP
PsKE21	=	KE21	;	%	Equilibrium constant	21	F6P<->G6P
PsKE22	=	KE22	;	% 	Equilibrium constant	22	G6P<->G1P
PsKM311	=	KM311	;	%	DHAP	31	DHAPi<->DHAPo
PsKM312	=	KM312	;	% 	Pi	31	DHAPi<->DHAPo
PsKM313	=	KM313	;	%	Pext	31	DHAPi<->DHAPo
PsKM32	=	KM32	;	% 	PGA	32	PGAi<->PGAo
PsKM33	=	KM33	;	%	GAP	33	GAPi<->GAPo


PsKM231 = KM231	;
PsKM232	=   KM232;
PsKM233 =	KM233	;
PsKM234 =   KM234	;
PsKA231	=   KA231   ;
PsKI231  =	KI231	;
PsKVmo  =   KVmo    ;
PsKE23  =   KE23    ;

PsKM241 = KM241;
PsKM242 = KM242;
PsKE24  = KE24;
PsKE25  = KE25;


PsKE6   =   KE6;
PsKE7   =   KE7;
PsKE8   =   KE8;
PsKE9   =   KE9;
PsKE10   =   KE10;
PsKE13   =   KE13;
PsKE16   =   KE16;
PsKM103 =   KM103;
PsKM163 =   KM163;

% Initialize the PrVmax of the different reactions based on the global variables Vmax

global	V1;	%	(Harris & Koniger, 1997)	1	Rubisco	RuBP+CO2<->2PGA
global	V2;	%	(Harris & Koniger, 1997)	2	PGA Kinase	PGA+ATP <-> ADP + DPGA
global	V3;	%	(Harris & Koniger, 1997)	3	GAP dehydragenase	DPGA+NADPH <->GAP + OP+NADP 
global	V4;	%	(Harris & Koniger, 1997)	4	Triose phosphate isomerase	DHAP <->GAP
global	V5;	%	(Harris & Koniger, 1997)	5	Aldolase	GAP+DHAP <->FBP
global	V6;	%	(Harris & Koniger, 1997)	6	FBPase	FBP<->F6P+OP
global	V7;	%	(Harris & Koniger, 1997)	7	Transketolase	F6P+GAP<->E4P+Xu5P
global	V8;	%	(Harris & Koniger, 1997)	8	Aldolase	E4P+DHAP<->SBP
global	V9;	%	(Harris & Koniger, 1997)	9	SBPase	SBP<->S7P+OP
global	V10;	%	(Harris & Koniger, 1997)	10	Transketolase	S7P+GAP<->Ri5P+Xu5P
global	V11;	%	(Harris & Koniger, 1997)	11	Pentosephosphate isomerase	Ri5P<-->Ru5P
global	V12;	%	(Harris & Koniger, 1997)	12	Pentosephosphate epimerase	Xu5P<-->Ru5P
global	V13;	%	(Harris & Koniger, 1997)	13	Ribulosebiphosphate kinase	Ru5P+ATP<->RuBP+ADP
global	V16;	%	(Aflalo & Shavit, 1983, Davenport & McLeod, 1986)	16	ATP synthase	ADP+Pi<->ATP
global	V21;	%		                        21	Hexose phosphate isomerase	F6P<->G6P
global	V22;	%		                        22	Phosphoglucomutase	G6P<->G1P
global	V23;	%	(Latzko, Steup & Schachtele, 1981)	23	ATP + G-1P  --> ADPG + PPi
global	V31;	%	(Lilley, Chon, Mosbach & Heldt, 1977b)	31	Phosphate translocator	DHAPi<->DHAPo
global	V32;	%	(Lilley et al., 1977b)	32	Phosphate translocator	PGAi<->PGAo
global	V33;	%	(Lilley et al., 1977b)	33	Phosphate translocator	GAPi<->GAPo
global  V24; %    %   ADPG --> ADP + Gn

% Get the values of the global variables of Vmax for different reactions
global FIBF_PSPR_com;
global DPH;
global ATPActive;

RegFactor = 1;

PsV1	=	V1	;	%	1	Rubisco	RuBP+CO2<->2PGA
PsV2	=	V2  ;	%	2	PGA Kinase	PGA+ATP <-> ADP + DPGA
PsV3	=	V3	;	%	3	GAP dehydragenase	DPGA+NADPH <->GAP + OP+NADP 
PsV4	=	V4	;	%	4	Triose phosphate isomerase	DHAP <->GAP
PsV5	=	V5	;	%	5	Aldolase	GAP+DHAP <->FBP
PsV6	=	V6	;	%	6	FBPase	FBP<->F6P+OP
PsV7	=	V7	;	%	7	Transketolase	F6P+GAP<->E4P+Xu5P
PsV8	=	V8	;	%	8	Aldolase	E4P+DHAP<->SBP
PsV9	=	V9	;	%	9	SBPase	SBP<->S7P+OP
PsV10	=	V10	;	%	10	Transketolase	S7P+GAP<->Ri5P+Xu5P
PsV11	=	V11	;	%	11	Pentosephosphate isomerase	Ri5P<-->Ru5P
PsV12	=	V12	;	%	12	Pentosephosphate epimerase	Xu5P<-->Ru5P
PsV13	=	V13;	%	13	Ribulosebiphosphate kinase	Ru5P+ATP<->RuBP+ADP
PsV16	=	V16	;	%	16	ATP synthase	ADP+Pi<->ATP
PsV21	=	V21	;	%	21	Hexose phosphate isomerase	F6P<->G6P
PsV22	=	V22	;	%	22	Phosphoglucomutase	G6P<->G1P
PsV23	=	V23 ;%	23	ATP + G-1P -> ADPG + PPi
PsV31	=	V31	* RegFactor;	%	31	Phosphate translocator	DHAPi<->DHAPo
PsV32	=	V32	* RegFactor;	%	32	Phosphate translocator	PGAi<->PGAo
PsV33	=	V33	* RegFactor;	%	33	Phosphate translocator	GAPi<->GAPo
PsV24   =   V24;        % 24    ADPG --> ADP + Gn

global SUCS2PS_Pic;

global PSPR_SUCS_com;
if PSPR_SUCS_com ==1
    PsPEXT = SUCS2PS_Pic;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First here is one way of the redox regulation, assuming the regulation is instataneous.
global RedoxReg_RA_com;     % This part is essentially not used now. This part is left here only 
                            % in case that there are more work using the equilibrium of Thio with enzyme
                            % as a way to regulate enzyme activities. 

global Redox2PS_V6;
global Redox2PS_V9;
global Redox2PS_V13;
global Redox2PS_V16;


if RedoxReg_RA_com == 2
    PsV6 =  Redox2PS_V6;
    PsV9 =  Redox2PS_V9;
    PsV13 =  Redox2PS_V13;
    PsV16 =  Redox2PS_V16;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The second method of the redox regulation. 
% First set the global poool from the RROEA. 

global RROEA_Pool;
global RROEA_EPS_com;
global RROEA2PS_GAPDH;
global RROEA2PS_FBPase;
global RROEA2PS_SBPase;
global RROEA2PS_PRK;
global RROEA2PS_ATPase;
global RROEA2PS_ATPGPP;

if RROEA_EPS_com ==0 
    % ATPreg = ATP/3;
    % ATPreg = PGA/3;         % If there is no regulation of enzyme activity, some forcing needed to be added. 
else
    ATPreg = 1;

    GAPDHT = RROEA_Pool	    (	1	);
    FBPaseT = RROEA_Pool	(	2	);
    SBPaseT = RROEA_Pool	(	3	);
    PRKT = RROEA_Pool	    (	4	);
    ATPaseT = RROEA_Pool	(	5   );
    ATPGPPT = RROEA_Pool	(	6	);

    PsV3 = V3 * RROEA2PS_GAPDH/GAPDHT;
    PsV9 = V9 * RROEA2PS_SBPase/SBPaseT;
    PsV13 = V13 * RROEA2PS_PRK/PRKT;
    PsV16 = V16 * RROEA2PS_ATPase/ATPaseT;
    PsV23 = V23 * RROEA2PS_ATPGPP/ATPGPPT;
end

% Setting the concentration 

RuBP	=	PSs(1)	;
PGA	    =	PSs(2)	;
DPGA	=	PSs(3)	;
T3P	    =	PSs(4)	;   
ADPG	=	PSs(5)	;
FBP	    =	PSs(6)	;
E4P	    =	PSs(7)	;
S7P	    =	PSs(8)	;
SBP	    =	PSs(9)	;
ATP	    =	PSs(10)	;
NADPH	=	PSs(11)	;
CO2	    =	PSs(12)	;
O2	    =	PSs(13)	;
HexP    =   PSs(14);
PenP    =   PSs(15);    

% Assuming that the regulation exists no matter there is enzyme regulation or not. ATPReg is 
% used to regulate the TP export and starch synthesis. 


% Now Calculate the concentration of the auxiliary variables.

global RuACT_EPS_com;

global PS2EPS_NADPH;
PS2EPS_NADPH = NADPH; 

global PR_PS_com;    % This is a variable indicating whether the PR model is actually need to be combined with PS or not. If 1 then means combined; 0 means not. 

global StomCond_TrDynaPS_com;       % Notice here only if there is no stomata conductance we need to use the 
                                    % external CO2 directly.                                   

if StomCond_TrDynaPS_com ==0
    
    global O2_cond;
    global CO2_cond;
    
    CO2 = CO2_cond;
    O2 = O2_cond; 
end

DHAP=  T3P/(1+KE4);
GAP =  KE4*T3P/(1+KE4);

ATPreg = DHAP/3;         % If there is no regulation of enzyme activity, some forcing needed to be added. 
if ATPreg >1
    ATPreg =1;
end

NADP = PS_C_CN - NADPH;
ADP =  PS_C_CA - ATP ;

F6P = (HexP /PsKE21)/(1+1/PsKE21+PsKE22);
G6P = HexP /(1+1/PsKE21+PsKE22);
G1P = (HexP *PsKE22)/(1+1/PsKE21+PsKE22);

Ru5P = PenP /(1+1/PsKE11+1/PsKE12);
Ri5P = (PenP/PsKE11) /(1+1/PsKE11+1/PsKE12);
Xu5P = (PenP/PsKE12)  /(1+1/PsKE11+1/PsKE12);

PR2PS_Pgca = Param(2);
Pit = PS_C_CP - PGA - 2*DPGA-GAP-DHAP-2*FBP-F6P-E4P-2*SBP-S7P-Xu5P-Ri5P-Ru5P-2*RuBP-G6P-G1P- ATP -PR2PS_Pgca;
Pi = 0.5 * (-PsKE25 + (PsKE25 * PsKE25 + 4 * Pit * PsKE25)^0.5);
OPOP = Pit - Pi;

global V1Reg;
global RUBISCOMETHOD;
global RUBISCOTOTAL;

global ATPActive;

global RedoxReg_RA_com;

if RedoxReg_RA_com ==0
    ATPreg = PGA/3;   
else
    ATPreg = 1;
end

V1Reg = 1+PGA/PsKI11+FBP/PsKI12+SBP/PsKI13 + Pi/PsKI14+NADPH/PsKI15;

if RUBISCOMETHOD ==2
    tmp = PsV1 * RuBP/(RuBP+PsKM13*V1Reg);
    v1 = tmp*CO2/(CO2+PsKM11*(1+O2/PsKM12));
    
   if RuBP<PsV1/2.5 
        v1 = v1 * RuBP/(PsV1/2.5);
   end

elseif RUBISCOMETHOD==1
    v1 = PsV1*CO2/(CO2+PsKM11*(1+O2/PsKM12));   
    if RuBP<PsV1/2.5
        v1 = v1 * RuBP/(PsV1/2.5);
    end
end      
 
v2 = PsV2 * PGA * ATP /((PGA + PsKM21)*(ATP+PsKM22*(1+ADP/PsKM23)));
v3 = PsV3 * DPGA * NADPH/((DPGA+PsKM31a)*(NADPH+PsKM32b));
v5 = PsV5 * (GAP * DHAP-FBP/PsKE5)/((PsKM51*PsKM52)*(1+GAP/PsKM51+DHAP/PsKM52+FBP/PsKM53+GAP*DHAP/(PsKM51*PsKM52)));
v8 = PsV8 * (DHAP *E4P-SBP/PsKE8)/((E4P+PsKM82)*(DHAP + PsKM81));

KE57 = 1.005 * 0.1;
Km8p5p = 0.118;
Km5p5p = 0.616;
KE810 = 0.8446;
Km5gap = 0.2727;
Km8f6p = 0.5443;
Km8s7p = 0.01576;
Km8gap = 0.09;
Den = 1 + (1+GAP/Km5gap)*(F6P/Km8f6p+S7P/Km8s7p)+GAP/Km8gap + 1/Km8p5p*(Xu5P*(1+E4P*Ri5P/Km5p5p)+E4P+Ri5P);

v7 = PsV7 * (F6P * GAP *KE57 - E4P * Xu5P)/(Km8p5p*Km5p5p*Den);
v10 = PsV7 * (S7P * GAP * KE810 - Xu5P * Ri5P)/(Km8p5p*Km5p5p*Den);
v6 = PsV6 * (FBP-F6P * Pi/PsKE6)/(FBP + PsKM61*(1+F6P/PsKI61+Pi/PsKI62));
v9 = PsV9 * (SBP-Pi * S7P/PsKE9) /(SBP + PsKM9*(1+Pi/PsKI9));
v13 = PsV13 * (ATP * Ru5P-ADP * RuBP/PsKE13)/((ATP*(1+ADP/PsKI134) + PsKM132*(1+ADP/PsKI135))*(Ru5P+PsKM131*(1+PGA/PsKI131+RuBP/PsKI132+Pi/PsKI133)));
v16 = PsV16 * (ADP * Pi-ATP/PsKE16)/(PsKM161*PsKM162 * (1+ADP/PsKM161 + Pi/PsKM162 + ATP/PsKM163 + ADP * Pi /(PsKM161 * PsKM162)));

Va = PsKVmo + PsV23 * (PGA/(PsKA231*(1+PGA/PsKA231)));
v23num = Va * (ATP * G1P - ADPG * OPOP/PsKE23);         % The reason we set this here is to assume that we can obtain a reverse reaction here. However, a more realistic
                                                        % way to achieve the homeostasis might be to allow starch breakdown and allow regulation of SBPase and FBPase. 

% WY 201803
%v23den = PsKM231 * PsKM232 * (1 + ATP/PsKM232 + G1P/PsKM231 + ATP * G1P/(PsKM231 * PsKM232) + ADPG/PsKM233 + OPOP/PsKM234 + ADPG * OPOP/(PsKM233 * PsKM234) + Pi/PsKI231);
v23den2 = (1+ Pi/PsKI231)*PsKM231 * PsKM232 * (1 + ATP/PsKM232 + G1P/PsKM231 + ATP * G1P/(PsKM231 * PsKM232) + ADPG/PsKM233 + OPOP/PsKM234 + ADPG * OPOP/(PsKM233 * PsKM234) );
v23 = v23num/v23den2;
%v23 = v23num/v23den;
v24num = PsV24 * (ADPG);                % Similar to the argument for reaction 23. The control of homeostasis might be better enforced at the point of reaction 23. 

v24dem = PsKM241 * (1+ ADPG/PsKM241);
v24 = v24num/v24dem;

MaxCoeff = 5; 
V25max =  0.5; 

v25 = V25max * (1- RuBP/MaxCoeff) * ATP/(ATP + 1); 
N = 1 + (1+ PsKM313/PsPEXT)*(Pi/PsKM312+PGA/PsKM32+GAP/PsKM33+DHAP/PsKM311);
v31 = PsV31 * DHAP/(N*PsKM311)  ;
v32 = PsV32 * PGA/(N*PsKM32);
v33 = PsV33 * GAP/(N * PsKM33);

v31 = v31 * ATPreg ;
v32 = v32 * ATPreg ;
v33 = v33 * ATPreg ;


global EPS_ADP;         % This variable is used in the BF_Rate when EPS is used. 
EPS_ADP = ADP;

global EPS_Pi;          % This variable is used in the BF_Rate when EPS is used. 
EPS_Pi = Pi;


global PS2EPS_V16;     
PS2EPS_V16 = v16;

global EPS_ADP;         % This variable is used in the BF_Rate when EPS is used. 
EPS_ADP = ADP;

global EPS_Pi;          % This variable is used in the BF_Rate when EPS is used. 
EPS_Pi = Pi;

global PS2EPS_v3;          
PS2EPS_v3 = v3;

global PSPR_RA_O2;          % RA is the combined EPS and Rubisco activase. 
global PSPR_RA_CO2;
PSPR_RA_O2 = O2;
PSPR_RA_CO2 = CO2;

global PS2RA_ATP;
PS2RA_ATP = ATP;

% information is sent back to PR by PS2PR_Pi global variable.

global PS2PR_Pi;
PS2PR_Pi = Pi;

global PS2PR_ATP;
PS2PR_ATP = ATP;

global PS2PR_ADP;
PS2PR_ADP = ADP;

global PS2BF_ATP;
PS2BF_ATP = ATP;

global PS2BF_ADP;
PS2BF_ADP = ADP;

global PS2BF_Pi;
PS2BF_Pi = Pi;


global PS2SUCS_PGA; 
PS2SUCS_PGA = PGA; 

% Notice the series PS2CM is used both in the CM model and the FPSReg model and thereafter. 

global PS2CM_vdhap;
PS2CM_vdhap = v31;

global PS2CM_vpga;
PS2CM_vpga = v32;

global PS2CM_vgap;
PS2CM_vgap = v33;

global PS2PRC; 
global PS2PRO;
PS2PRC = CO2; 
PS2PRO = O2; 

global PS2RubACC;
global PS2RubACO;
PS2RubACC = CO2; 
PS2RubACO = O2; 
        
global PS2Stom_CO2_consum;
PS2Stom_CO2_consum = v1; 

PSr = zeros(16,1);

PSr(1)	=	v1	;
PSr(2)	=	v2	;
PSr(3)	=	v3	;
PSr(4)	=	0	;
PSr(5)	=	v5	;
PSr(6)	=	v6	;
PSr(7)	=	v7	;
PSr(8)	=	v8	;
PSr(9)	=	v9	;
PSr(10)	=	v10	;
PSr(11)	=	v13	;
PSr(12)	=	v16	;
PSr(13)	=	v23	;
PSr(14)	=	v31	;
PSr(15) =   v32;
PSr(16) =   v33;
PSr(17) =   v24;
PSr(18) =   v25; 
PSr;
% Getting the information for output as figures.

if (PS_TIME_N ==0)
    PS_TIME_N = 1;
end

if (t > PS_OLD_TIME)
    PS_TIME_N = PS_TIME_N + 1;
    PS_OLD_TIME = t;
end

PS_VEL(1,PS_TIME_N) = t;

PS_VEL(2,PS_TIME_N) = v1;
PS_VEL(3,PS_TIME_N) = v2;
PS_VEL(4,PS_TIME_N) = v3;
PS_VEL(5,PS_TIME_N) = 0;
PS_VEL(6,PS_TIME_N) = v5;
PS_VEL(7,PS_TIME_N) = v6;
PS_VEL(8,PS_TIME_N) = v7;
PS_VEL(9,PS_TIME_N) = v8;
PS_VEL(10,PS_TIME_N) = v9;
PS_VEL(11,PS_TIME_N) = v10;
PS_VEL(12,PS_TIME_N) = v13;
PS_VEL(13,PS_TIME_N) = v16;
PS_VEL(14,PS_TIME_N) = v23;
PS_VEL(15,PS_TIME_N) = v31;
PS_VEL(16,PS_TIME_N) = v32;
PS_VEL(17,PS_TIME_N) = v33;
PS_VEL(18,PS_TIME_N) = Pi;
PS_VEL(19,PS_TIME_N) = v24;
PS_VEL(20,PS_TIME_N) = v25;



% Transfer the variables for output

global PS2OUT;
PS2OUT = zeros(5,1);

PS2OUT(1)	=   RuBP;
PS2OUT(2)	=   PGA;
PS2OUT(3)	=   DPGA;
PS2OUT(4)	=   T3P;   
PS2OUT(5)	=   ADPG;
PS2OUT(6)	=   FBP;
PS2OUT(7)	=   E4P;
PS2OUT(8)	=   S7P;
PS2OUT(9)	=   SBP;
PS2OUT(10)	=   ATP;
PS2OUT(11)	=   NADPH;
PS2OUT(12)	=   CO2;    
PS2OUT(13)	=   O2;     
PS2OUT(14)  =   HexP;   
PS2OUT(15)  =   PenP;        
PS2OUT(16)  =   Pi;    
PS2OUT(17)  =   ADP;    
PS2OUT(18) =    v1;   