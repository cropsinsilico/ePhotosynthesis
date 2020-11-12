function CMr = CM_Rate(t,CM_Con, SUCS_Param)
global PS2EPS_V16; 
global PS2EPS_v3;
global CO2_cond;
global VfactorC;
global VfactorT;
global GRNC;
global GRNT;
global Tp;
global pcfactor;
Vfactor1=1;
Vfactor2=1;
Vfactor3=1;
Vfactor4=1;
Vfactor5=1;
Vfactor6=1;
Vfactor7=1;
Vfactor9=1;
Vfactor11=1;
Vfactor12=1;
Vfactor124=1;
Vfactor13=1;
Vfactor21=1;
Vfactor22=1;
Vfactor23=1;
Vfactor25=1;
Vfactor52=1;
Vfactor53=1;
Vfactor54=1;
Vfactor56=1;
Vfactor57=1;
Vfactor59=1;
Vfactor112=1;
Vfactor113=1;
Vfactor121=1;
Vfactor122=1;
Vfactor123=1;
Vfactor131=1;
Vfactor51=1;
if GRNC==1
if CO2_cond>0;%0.0107
VfactorCp=VfactorC;
    
%     if pcfactor==1.059;
%     Vfactor1=VfactorCp(1)/1.059;
%    else
%     Vfactor1=VfactorCp(1);
%    end
Vfactor1=VfactorCp(1);
Vfactor2=VfactorCp(3);
Vfactor3=VfactorCp(4);
Vfactor4=VfactorCp(5);
Vfactor5=VfactorCp(6);
Vfactor6=VfactorCp(7);
Vfactor7=VfactorCp(8);
Vfactor9=VfactorCp(9);
Vfactor11=VfactorCp(10);
Vfactor12=VfactorCp(11);
Vfactor124=VfactorCp(12);
Vfactor13=VfactorCp(13);
Vfactor21=VfactorCp(14);
Vfactor22=VfactorCp(15);
Vfactor23=VfactorCp(17);
Vfactor25=VfactorCp(19);
Vfactor52=VfactorCp(20);
Vfactor53=VfactorCp(21);
Vfactor54=VfactorCp(22);
Vfactor56=VfactorCp(23);
Vfactor57=VfactorCp(24);
Vfactor59=VfactorCp(25);
Vfactor112=VfactorCp(26);
Vfactor113=VfactorCp(27);
Vfactor121=VfactorCp(29);
Vfactor122=VfactorCp(30);
Vfactor123=VfactorCp(31);
Vfactor131=VfactorCp(33);
end
end

Vf_T52=1;
Vf_T5=1;
Vf_T23=1;
Vf_T13=1;
Vf_T131=1;
Vf_T25=1;
Vf_T113=1;
Vf_T123=1;
Vf_T121=1;
Vf_T3=1;
Vf_T122=1;
Vf_T59=1;
Vf_T57=1;
Vf_T12=1;
Vf_T4=1;
Vf_T6=1;
Vf_T51=1;
Vf_T53=1;
Vf_T54=1;
Vf_T2=1;
Vf_T112=1;
Vf_T11=1;
Vf_T1=1;
Vf_T9=1;
%Vf_T23=1;
Vf_T56=1;
if GRNT==1
if Tp>25
Vf_T52=VfactorT(1);
Vf_T5=VfactorT(2);
Vf_T23=VfactorT(3);
Vf_T13=VfactorT(4);
Vf_T131=VfactorT(5);
Vf_T25=VfactorT(6);
Vf_T113=VfactorT(7);
Vf_T123=VfactorT(8);
Vf_T121=VfactorT(9);
Vf_T3=VfactorT(10);
Vf_T122=VfactorT(12);
Vf_T59=VfactorT(13);
Vf_T57=VfactorT(14);
Vf_T12=VfactorT(15);
Vf_T4=VfactorT(16);
Vf_T6=VfactorT(17);
Vf_T51=VfactorT(18);
Vf_T53=VfactorT(20);
Vf_T54=VfactorT(21);
Vf_T2=VfactorT(22);
Vf_T112=VfactorT(23);
Vf_T11=VfactorT(24);
Vf_T1=VfactorT(25);
Vf_T9=VfactorT(26);
%Vf_T23=VfactorT(27);
Vf_T56=VfactorT(28);
end
end
global LI;
global Jmax;
global alfa;
global fc;
global Theta;
global beta;

% The global constant used for conservation law

global PS_C_CA;             %   Global constant for the total adenylates
global PS_C_CP;             %   Global constant for the total phosphate
global PS_C_CN;             %   Global constant for the total NADP+NADPH
global NADHc;
global NADc;
global GLUc;
global KGc;

% Variables for flux output

global PS_VEL;
global PS_TIME_N;
global PS_OLD_TIME;
% The following defines the global variables used for output velocity and concentration information.
global PR_OLD_TIME;
global PR_TIME_N;
global PR_VEL;
global SUCS_OLD_TIME;
global SUCS_TIME_N;
global SUCS_VEL;
global SUCS_CON;

R=8.314;
c_c=38.28;
dHa_c=80.99;
c_o=14.68;
dHa_o=23.72;
% 
% First the physical and chemical constant for all the reactions

%PsKM11_0	=	0.0115;		% 	CO2	1	RuBP+CO2->2PGA
%%%%%%%%%%soy
PsKM11_0	=	0.0097;
PsKM12_0	=	0.244;		%	O2	1	RuBP+CO2->2PGA  0.28 DEFAUL. 
PsKM11=PsKM11_0*exp(c_c-dHa_c*1000/(R*(Tp+273.15)))/272.38;
PsKM12=PsKM12_0*exp(c_o-dHa_o*1000/(R*(Tp+273.15)))/165.82;
PsKM13	=	0.02;		% 	RuBP	1	RuBP+CO2->2PGA
PsKI11    =   0.84   ;    % PGA  
PsKI12    =0.04   ;       % FBP
PsKI13    = 0.075 ;       % SBP
PsKI14    = 0.9   ;       % Pi
PsKI15    = 0.07  ;       % NADPH
PsKM21	=	0.240;		%	PGA	2	PGA+ATP <-> ADP + DPGA
PsKM22	=	0.390;		% 	ATP	2	PGA+ATP <-> ADP + DPGA
PsKM23    =   0.23  ;        %  ADP     
PsKM31a	=	0.004;		%	BPGA	3	DPGA+NADPH <->GAP + OP+NADP 
PsKM32b	=	0.1	;	    % 	NADPH	3	DPGA+NADPH <->GAP + OP+NADP
KM41	=	2.5	;	    %	DHAP	4	DHAP <->GAP
KM42	=	0.68;		% 	GAP	4	DHAP <->GAP
KE4     =   0.05;       %   Using the value from Patterson
PsKM51	=	0.3	;	    %	GAP	5	GAP+DHAP <->FBP
PsKM52	=	0.4	;	    % 	DHAP	5	GAP+DHAP <->FBP
PsKM53	=	0.02;		%	FBP	5	GAP+DHAP <->FBP     % Original Value: 0.02
PsKE5     = 7.100;          % Defult: 7.1
PsKM61	=	0.033;		% 	FBP	6	FBP<->F6P+OP
PsKI61    = 0.7   ;       %   F6P       
PsKI62    = 12    ;       %   Pi
PsKE6     =   6.66 * 10^5;    % The equilibrium constant for this reaction        % New    mM     Laisk or Bassham and Krause 1969 BBA
PsKM71	=	0.100;		%	Xu5P	7	F6P+GAP<->E4P+Xu5P      % jn
PsKM72	=	0.100;		% 	E4P	7	F6P+GAP<->E4P+Xu5P
PsKM73    = 0.1;         %   F6P This value was based on estimate
PsKM74    = 0.1000;         % Estimate for GAP ORIGINAL 0.1
PsKE7     =   10 ;       % The equilibrium constant for this reaction             % New           Laisk  Bassham and Krause 1969 BBA
PsKM8	    =	0.02;		%	SBP	8	E4P+DHAP<->SBP
PsKM81    = 0.4   ;       % DHAP
PsKM82    = 0.2   ;       % E4P estimate
PsKE8     = 1.017 ;     % The equilibrium constant for this reaction                  % New    mM-1         Laisk  Bassham and Krause 1969 BBA. Default: 1.107
PsKM9	    =	0.05;		% 	SBP	9	SBP<->S7P+OP    
PsKI9     = 12    ;       %   The inibintion constant for Pi; 
PsKE9     =   6.66 * 10^5 ; % The equilibrium constant of this reaction           % New   mM      Laisk  Bassham and Krause 1969 BBA
PsKM10	=	1.5	;	    %	R5P	10	S7P+GAP<->Ri5P+Xu5P
PsKM101   =   0.1 ;       %   Xu5P
PsKM102   = 0.072 ;       %   Estimate for GAP
PsKM103   = 0.46 ;        %   Estimate for S7P                                    % New 
PsKE10    = 1/0.85 ;      %   The equilibrium constant for this reaction          % New From Laisk or Bassham and Krause 1969 BBA
PsKE11	=	0.4	;	    %	Equilibrium Constant	11	Ri5P<-->Ru5P
PsKE12	=	0.67;		% 	Equilibrium Constant	12	Xu5P<-->Ru5P

PsKM131	=	0.05;		    %	Ru5P	13	Ru5P+ATP<->RuBP+ADP
PsKM132	=	0.059;		    % 	ATP	13	Ru5P+ATP<->RuBP+ADP
PsKI131	=	2	;			%	PGA	13	Ru5P+ATP<->RuBP+ADP
PsKI132	=	0.7	;			%	RuBP	13	Ru5P+ATP<->RuBP+ADP
PsKI133	=	4	;			%	Pi	13	Ru5P+ATP<->RuBP+ADP
PsKI134	=	2.5	;			%	ADP	13	Ru5P+ATP<->RuBP+ADP
PsKI135	=	0.4	;			%	ADP	13	Ru5P+ATP<->RuBP+ADP
PsKE13    =   6.846 * 10^3;   %   The equilibrium constant for this reaction  % New From Laisk or Bassham and Krause 1969 BBA
PsKM161	=	0.014;		%	ADP	16	ADP+Pi<->ATP
PsKM162	=	0.3;		% 	Pi	16	ADP+Pi<-> ATP
PsKM163   =   0.3;        %   ATP 16  ADP+Pi<-> ATP                           % New       Based on Laisk  
PsKE16    =   5.734;      %   The equilibrium constant for this reaction      % NEW, From Laisk or Bassham and Krause 1969 BBA

PsKE21	=	2.3;		%	Equilibrium constant	21	F6P<->G6P
PsKE22	=	0.058;		% 	Equilibrium constant	22	G6P<->G1P
PsKM231	=	0.08;		%	G1P	23	G1P+ATP+Gn<->PPi+ADP+Gn+1
PsKM232	=	0.08;		% 	ATP	23	G1P+ATP+Gn<->PPi+ADP+Gn+1
PsKA231	=	0.1;		%	PGA	23	G1P+ATP+Gn<->PPi+ADP+Gn+1
PsKA232	=	0.02;		% 	F6P	23	G1P+ATP+Gn<->PPi+ADP+Gn+1
PsKA233	=	0.02;		%	FBP	23	G1P+ATP+Gn<->PPi+ADP+Gn+1
PsKI23	=	10;		    % 	ADP	23	G1P+ATP+Gn<->PPi+ADP+Gn+1

PsKM311	=	0.077;		%	DHAP	31	DHAPi<->DHAPo
PsKM312	=	0.63;		% 	Pi	31	DHAPi<->DHAPo
PsKM313	=	0.74;		%	Pext	31	DHAPi<->DHAPo
PsKM32	=	0.25;		% 	PGA	32	PGAi<->PGAo
PsKM33	=	0.075;		%	GAP	33	GAPi<->GAPo

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
global	V23;	%	(Latzko, Steup & Schachtele, 1981)	23	ADP-glucose pyrophosphorylase and	ADPG+Gn<->G(n+1)+ADP + 2Pi
global	V31;	%	(Lilley, Chon, Mosbach & Heldt, 1977b)	31	Phosphate translocator	DHAPi<->DHAPo
global	V32;	%	(Lilley et al., 1977b)	32	Phosphate translocator	PGAi<->PGAo
global	V33;	%	(Lilley et al., 1977b)	33	Phosphate translocator	GAPi<->GAPo

% Get the values of the global variables of Vmax for different reactions
global FIBF_PSPR_com;
%global DPH;
global ATPActive;

PsV1_0	=	V1*Vfactor1*Vf_T1	;	%	1	Rubisco	RuBP+CO2<->2PGA
PsV2_0	=	V2*Vfactor2*Vf_T2 ;	%	2	PGA Kinase	PGA+ATP <-> ADP + DPGA
PsV3_0	=	V3*Vfactor3*Vf_T3	;	%	3	GAP dehydragenase	DPGA+NADPH <->GAP + OP+NADP 
PsV4	=	V4*Vfactor4*Vf_T4	;	%	4	Triose phosphate isomerase	DHAP <->GAP
PsV5_0	=	V5*Vfactor5*Vf_T5;	%	5	Aldolase	GAP+DHAP <->FBP
PsV6_0	=	V6*Vf_T6	;	%	6	FBPase	FBP<->F6P+OP
PsV7_0	=	V7*Vfactor7	;	%	7	Transketolase	F6P+GAP<->E4P+Xu5P
PsV8_0	=	V8*Vfactor5*Vf_T5	;	%	8	Aldolase	E4P+DHAP<->SBP
PsV9_0	=	V9*Vf_T9	;	%	9	SBPase	SBP<->S7P+OP
PsV10_0	=	V10*Vfactor7	;	%	10	Transketolase	S7P+GAP<->Ri5P+Xu5P
PsV11	=	V11*Vfactor11*Vf_T11	;	%	11	Pentosephosphate isomerase	Ri5P<-->Ru5P
PsV12	=	V12*Vfactor12*Vf_T12	;	%	12	Pentosephosphate epimerase	Xu5P<-->Ru5P
PsV13_0	=	V13*Vfactor13*Vf_T13;	%	13	Ribulosebiphosphate kinase	Ru5P+ATP<->RuBP+ADP
PsV16	=	V16	;	%	16	ATP synthase	ADP+Pi<->ATP
PsV21	=	V21*Vfactor21	;	%	21	Hexose phosphate isomerase	F6P<->G6P
PsV22	=	V22	;	%	22	Phosphoglucomutase	G6P<->G1P
PsV23_0	=	V23*Vfactor23*Vf_T23;%	23	ADP-glucose pyrophosphorylase and	ADPG+Gn<->G(n+1)+ADP
PsV31	=	V31 ;	%	31	Phosphate translocator	DHAPi<->DHAPo
PsV32	=	V32	;	%	32	Phosphate translocator	PGAi<->PGAo
PsV33	=	V33	;	%	33	Phosphate translocator	GAPi<->GAPo


global V111;
global V112;
global V113;
global V121;
global V122;
global V123;
global V124;
global V131;
global V1T;
global V2T;

PrV111_0 = V111;
PrV112_0 = V112*Vfactor112*Vf_T112;
PrV113_0 = V113*Vfactor113*Vf_T113;
PrV121_0 = V121*Vfactor121*Vf_T121;
PrV122_0 = V122*Vfactor122*Vf_T122;
PrV123_0 = V123*Vfactor123*Vf_T123;
PrV124_0 = V124*Vfactor124;
PrV131_0 = V131*Vfactor131*Vf_T131;
PrV1T = V1T;
PrV2T = V2T;

global 	V51	;
global 	V52	;
global 	V55	;
global 	V56	;
global 	V57	;
global 	V58	;
global 	V59	;
global  V60;
global  V61;
global	V62;	

SUCSV51_0	=	V51*Vfactor51*Vf_T51	;	%	;		DHAP+GAP --FBP
SUCSV52_0	=	V52*Vfactor52*Vf_T52	;	%	;		FBP --F6P + Pi
SUCSV55_0	=	V55	;	%	;		G1P+UTP --OPOP+UDPG 
SUCSV56_0	=	V56*Vfactor56*Vf_T56	;	%	;		UDPG+F6P--SUCP + UDP
SUCSV57_0	=	V57*Vfactor57*Vf_T57	;	%	;		SUCP--Pi + SUC
SUCSV58_0	=	V58	;	%	;		F26BP--F6P + Pi
SUCSV59	=	V59*Vfactor59*Vf_T59	;	%	;		F6P + ATP --ADP + F26BP
SUCSV60	=	V60	;	%	;		ATP+UDP --UTP + ADP
SUCSV61	=	V61	;	%	;		POPO --2PO
SUCSV62	=	V62	;	%	;		SUC Sink

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%WY Temp response

Q10_1=1.93;
Q10_2=2;
Q10_3=2;
Q10_5=2;
Q10_6=2;
Q10_7=2;
Q10_8=2;
Q10_9=2;
Q10_10=2;
Q10_13=2;
Q10_23=2;
Q10_112=1.81;
Q10_113=2;
Q10_121=2;
Q10_122=2.01;
Q10_123=2;
Q10_124=2;
Q10_131=2;
Q10_51=2;
Q10_52=1.60;
Q10_55=2;
Q10_56=2;
Q10_57=2;
Q10_58=2;

Ru_Act=-3E-05*Tp^3 + 0.0013*Tp^2 - 0.0106*Tp + 0.8839;%Rubisco activition state
PsV1 =PsV1_0*Ru_Act*Q10_1^((Tp-25)/10);
PsV2 =PsV2_0*Q10_2^((Tp-25)/10);
PsV3 =PsV3_0*Q10_3^((Tp-25)/10);
PsV5 =PsV5_0*Q10_5^((Tp-25)/10);
PsV6 =PsV6_0*Q10_6^((Tp-25)/10);
PsV7 =PsV7_0*Q10_7^((Tp-25)/10);
PsV8 =PsV8_0*Q10_8^((Tp-25)/10);
PsV9 =PsV9_0*Q10_9^((Tp-25)/10);
PsV10=PsV10_0*Q10_10^((Tp-25)/10);
PsV13=PsV13_0*Q10_13^((Tp-25)/10);
PsV23=PsV23_0*Q10_23^((Tp-25)/10);
PrV112=PrV112_0*Q10_112^((Tp-25)/10);
PrV113=PrV113_0*Q10_113^((Tp-25)/10);
PrV121=PrV121_0*Q10_121^((Tp-25)/10);
PrV122=PrV122_0*Q10_122^((Tp-25)/10);
PrV123=PrV123_0*Q10_123^((Tp-25)/10);
PrV124=PrV124_0*Q10_124^((Tp-25)/10);
PrV131=PrV131_0*Q10_131^((Tp-25)/10);
SUCSV51=SUCSV51_0*Q10_51^((Tp-25)/10);
SUCSV52=SUCSV52_0*Q10_52^((Tp-25)/10);
SUCSV55=SUCSV55_0*Q10_55^((Tp-25)/10);
SUCSV56=SUCSV56_0*Q10_56^((Tp-25)/10);
SUCSV57=SUCSV57_0*Q10_57^((Tp-25)/10);
SUCSV58=SUCSV58_0*Q10_58^((Tp-25)/10);

PrV111= PsV1* 0.24;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Here is the region to transfer parameters


% Two different way of redox regulation of enzyme activities. These are not used in the current simulations. They are used in the full photosynthesis model. 
% First here is one way of the redox regulation, assuming the regulation is instataneous.

global Redox2PS_V6;
global Redox2PS_V9;
global Redox2PS_V13;
global Redox2PS_V16;

global RedoxReg_RA_com;

if RedoxReg_RA_com == 2
    PsV6 =  Redox2PS_V6;
    PsV9 =  Redox2PS_V9;
    PsV13 =  Redox2PS_V13;
    PsV16 =  Redox2PS_V16;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 1 Get the rate constant and the initial concentrations % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%global SUCS_RC;

light = SUCS_Param(1);

KE501	=	0.05	;	%	Equilibrium Constant		50		KE501		0.05		[Bassham, 1869 #832]
Km511	=	0.02	;	%	FBP	4.1.2.13	51		Km511	FBP	0.02	Pisum sativum	(Anderson, Heinrikson et al. 1975)
Km512	=	0.3	;	%	FBP	4.1.2.13	51		Km512	GAP	0.3	Spinacia oleracea	(Iwaki, Wadano et al. 1991)
Km513	=	0.4	;	%	FBP	4.1.2.13	51		Km513	DHAP	0.4	Spinacia oleracea	(Iwaki, Wadano et al. 1991)
KE51	=	12;     %  Based on Thomas et al 1997 Biochem Journal. The fifth citation in the paper. 
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
KE61	=	1.2*107	;	%	Pyrophosphate hydrolysis				KE61		1.2*107		{Flodgaard, 1974 #2521}
Km621	=	5	;	%	Vsink			Notice: pH dependent	Km621	Sucrose	5		{Weschke, 2000 #2522}



% SUCSVdhap_in  = Vdhap_in;    %   DHAP export from chloroplast
% SUCSVgap_in   = Vgap_in;    %   GAP export from chloroplast
% SUCSVpga_in   = Vpga_in;    %   PGA export from chloropalst


% The rate constant used in the model	

global SUCS_Pool;		

ATc = SUCS_Pool	(	1	)	;
UTc = SUCS_Pool	(	2	)		;
PTc = SUCS_Pool	(	3	)		;
% Setting the initial concentrations 

RuBP	=	CM_Con(1)	;
PGA	    =	CM_Con(2)	;
DPGA	=	CM_Con(3)	;
T3P	    =	CM_Con(4)	;   
FBP	    =	CM_Con(5)	;
E4P	    =	CM_Con(6)	;
S7P	    =	CM_Con(7)	;
SBP	    =	CM_Con(8)	;
ATP	    =	CM_Con(9)	;
NADPH	=	CM_Con(10)	;
CO2	    =	CM_Con(11)	;
O2	    =	CM_Con(12)	;
HexP    =   CM_Con(13);
PenP    =   CM_Con(14);  

Gcea = CM_Con(15);
Gca = CM_Con(16);
Pga = PGA;
Pgca = CM_Con(17);
Gcac = CM_Con(18);
Goac = CM_Con(19);
Serc = CM_Con(20);
Glyc = CM_Con(21);
Hprc = CM_Con(22);
Gceac = CM_Con(23);
Rubp = RuBP;
C = CO2;
O = O2;

T3Pc	=	CM_Con	(24);	
FBPc	=	CM_Con	(25);	
HexPc	=	CM_Con	(26);
F26BPc	=	CM_Con	(27);	
ATPc	=	CM_Con	(28);	
ADPc	=	CM_Con	(29);	
OPOPc	=	CM_Con	(30);
UDPGc	=	CM_Con	(31);	
UTPc	=	CM_Con	(32);	
SUCP	=	CM_Con	(33);
SUC	    =	CM_Con	(34);
PGAc	=	CM_Con	(35);	

global PS2BF_ADP;
global PS2BF_Pi;
global EPS_NADPH;


global PR_PS_com;    % This is a variable indicating whether the PR model is actually need to be combined with PS or not. If 1 then means combined; 0 means not. 

global CO2_cond;
global O2_cond;

CO2 = CO2_cond;
O2 = O2_cond;

DHAP= T3P/(1+KE4);
GAP = KE4*T3P/(1+KE4);

NADP = PS_C_CN - NADPH;
ADP = PS_C_CA - ATP ;

F6P = (HexP /PsKE21)/(1+1/PsKE21+PsKE22);
G6P = HexP /(1+1/PsKE21+PsKE22);
G1P = (HexP *PsKE22)/(1+1/PsKE21+PsKE22);

Ru5P = PenP /(1+1/PsKE11+1/PsKE12);
Ri5P = (PenP/PsKE11) /(1+1/PsKE11+1/PsKE12);
Xu5P = (PenP/PsKE12)  /(1+1/PsKE11+1/PsKE12);

% Notice here the concentration of PGCA is transfered from the photorespiratory pathway. 
PR2PS_Pgca = CM_Con(17);
Pi = PS_C_CP - PGA - 2*DPGA-GAP-DHAP-2*FBP-F6P-E4P-2*SBP-S7P-Xu5P-Ri5P-Ru5P-2*RuBP-G6P-G1P- ATP-PR2PS_Pgca;
% HexP
TEMP = 1 + KE541 + 1/KE531;

G6Pc = HexPc /TEMP;
F6Pc = G6Pc/KE531;
G1Pc = HexPc * KE541 /TEMP;

% T3P
GAPc = T3Pc /(1 + KE501);
DHAPc = T3Pc * KE501/(1+KE501);

% UDP
UDPc = UTc - UTPc - UDPGc;
ADPc = ATc - ATPc; 
% OP
PiTc = PTc - 2 * ( FBPc + F26BPc) - (PGAc + T3Pc + HexPc + SUCP + UTPc + ATPc);        %   ???? Where the PGA coming from?
Pic   = ((KE61^2 + 4 * KE61 * PiTc)^0.5 - KE61)/2;
OPOPc = PiTc - Pic;
PsPEXT = Pic;

global V1Reg;
global RUBISCOMETHOD;
global RUBISCOTOTAL;


if RedoxReg_RA_com ==0
    ATPreg = ATP/3;
    ATPreg = PGA/3;   
else
    ATPreg = 1;
end
V1Reg = 1+PGA/PsKI11+FBP/PsKI12+SBP/PsKI13 + Pi/PsKI14+NADPH/PsKI15;

if RUBISCOMETHOD ==2
    tmp = PsV1 * RuBP/(RuBP+PsKM13*V1Reg);
    v1 = tmp*CO2/(CO2+PsKM11*(1+O2/PsKM12));
   if RuBP<PsV1/2
        v1 = v1 * RuBP/(PsV1/2);
    end

elseif RUBISCOMETHOD==1
    v1 = PsV1*CO2/(CO2+PsKM11*(1+O2/PsKM12));   
    if RuBP<PsV1/2
        v1 = v1 * RuBP/(PsV1/2);
    end
end      
 
v2 = PsV2 * PGA * ATP /((PGA + PsKM21)*(ATP+PsKM22*(1+ADP/PsKM23)));
v3 = PsV3 * DPGA * NADPH/((DPGA+PsKM31a)*(NADPH+PsKM32b));
v5 = PsV5 * (GAP * DHAP-FBP/PsKE5)/((PsKM51*PsKM52)*(1+GAP/PsKM51+DHAP/PsKM52+FBP/PsKM53+GAP*DHAP/(PsKM51*PsKM52)));
v8 = PsV8 * (DHAP *E4P-SBP/PsKE8)/((E4P+PsKM82)*(DHAP + PsKM81));
v6 = PsV6 * (FBP-F6P * Pi/PsKE6)/(FBP + PsKM61*(1+F6P/PsKI61+Pi/PsKI62));
v7 = PsV7 * (F6P*GAP-Xu5P * E4P/PsKE7)/((F6P + PsKM73*(1+Xu5P/PsKM71+E4P/PsKM72))*(GAP+PsKM74));
v9 = PsV9 * (SBP-Pi * S7P/PsKE9) /(SBP + PsKM9*(1+Pi/PsKI9));
v10 = PsV10 * (GAP *S7P - Ri5P * Xu5P/PsKE10)/((GAP+PsKM102*(1+Xu5P/PsKM101+Ri5P/PsKM10))*(S7P+PsKM103));
v13 = PsV13 * (ATP * Ru5P-ADP * RuBP/PsKE13)/((ATP*(1+ADP/PsKI134) + PsKM132*(1+ADP/PsKI135))*(Ru5P+PsKM131*(1+PGA/PsKI131+RuBP/PsKI132+Pi/PsKI133)));


I2=LI*alfa*(1-fc)/2;
J=(I2+Jmax-sqrt((I2+Jmax)^2-4*Theta*I2*Jmax))/(2*Theta);
v16 = min(beta*J,PsV16* (ADP * Pi-ATP/PsKE16)/(PsKM161*PsKM162 * (1+ADP/PsKM161 + Pi/PsKM162 + ATP/PsKM163 + ADP * Pi /(PsKM161 * PsKM162))));
%v16o = PsV16* (ADP * Pi-ATP/PsKE16)/(PsKM161*PsKM162 * (1+ADP/PsKM161 + Pi/PsKM162 + ATP/PsKM163 + ADP * Pi /(PsKM161 * PsKM162)));
PS2EPS_V16=v16;
v23 = PsV23 * G1P *ATP /((G1P+PsKM231)*((1+ADP/PsKI23)*(ATP+PsKM232)+(PsKM232*Pi/(PsKA231*PGA+PsKA232*F6P+PsKA233*FBP))));
N = 1 + (1+ PsKM313/PsPEXT)*(Pi/PsKM312+PGA/PsKM32+GAP/PsKM33+DHAP/PsKM311);

% The ATP regualtion really is implicit in the light regulation of sucrose synthesis. 
v31 = PsV31 * DHAP/(N*PsKM311)  ;
v32 = PsV32 * PGA/(N*PsKM32);
v33 = PsV33 * GAP/(N * PsKM33);

v23 = v23 * ATPreg;
v31 = v31 * ATPreg;
v32 = v32 * ATPreg;
v33 = v33 * ATPreg;

global EPS_ATP_Rate;

if FIBF_PSPR_com ==0          % ModelMethod = 0 means that there is no connection between FIBF and PSPR. 
    if v16 == 0             % This assmed that light reguate the export of triose phosphate export. This function should use 
                            % ATP as a signal. 
         v23 = 0;
    end
else
    if EPS_ATP_Rate ==0
         v23 = 0;
    end
end

%global PS2BF_Pi;
PS2BF_Pi = Pi;
EPS_NADPH = NADPH; 
PS2BF_ADP =    ADP;

% Notice the series PS2CM is used both in the CM model and the FPSReg model and thereafter. 

CMr = zeros(16,1);

CMr(1)	=	v1	;
CMr(2)	=	v2	;
CMr(3)	=	v3	;
CMr(4)	=	0	;
CMr(5)	=	v5	;
CMr(6)	=	v6	;
CMr(7)	=	v7	;
CMr(8)	=	v8	;
CMr(9)	=	v9	;
CMr(10)	=	v10	;
CMr(11)	=	v13	;
CMr(12)	=	v16	;
CMr(13)	=	v23	;
CMr(14)	=	v31	;
CMr(15) =   v32;
CMr(16) =   v33;

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
PS2EPS_v3=v3;

% Reaction: 111: RUBP+O2<-->PGlycolate + PGA
KO = PsKM12;           % Michaelis constant for O2
KC = PsKM11;          % Michaelis constant for CO2  
KR = 0.02;           % Michaelis constant for RUBP  

% Reaction: 112: PGlycolate-->Pi+Glycolate;
KM112 = 0.026;   % Km112 for PGlycolate;
KI1122 = 94;     % Inhibition constant for Glycolate;
KI1121 = 2.55;   % The competitive Pi inhibition for PGlycolate   

% Reaction 113  : Gcea+ATP<-->ADP + PGA
KM1131 = 0.21;  % Km for ATP;
KM1132 = 0.25;  % Km for Gcea;
KI113 = 0.36;   % Ki for ATP BY pga;  %%%%%%%%%%%%%%%%%%%%%%%%% Competitive inhibition for ATP; in original paper it is 0.36;
KE113 = 300;     % New       Kleczkowski et al . 1985 Archives of Biochemistry and Biophysics  

% To set global information for different reactions
KM121 = 0.1;%Glycolate +O2<-->H2O2+Glyoxylate

KM1221 = 0.15; % Michaelis constant for glyoxylate;
KM1222 = 2.7;  % Michaelis constant for serinie;
KI1221 = 33;   % Inhibition constant for Glycine;  
KE122 = 0.24;  %  New: Guynn, R.W.; Arch. Biochem. Biophys.; 218, 14 (1982).; 0.24. At 25 degree. 

KM123 = 0.09;       %   Michaelis constant for hydroxylpyruvate;
KI123 = 12;          % Inhibition constant for hydroxypyruvate;
KE123 = 1/(4*10^(-6));  % Guynn, R.W.; Arch. Biochem. Biophys.; 218, 14 (1982).; 1/(4*10^(-6);

KM1241 = 0.15; % Michaelis constant for glyoxylate
KM1242 = 1.7;  % Michaelis constant for Glu
KI124 = 2;     % This KI is one guessed.
KE124 = 607;   % New       Cooper, A.J.L.; Meister, A.; Biochemistry; 11, 661 (1972).; K' 607. 

KM1311 = 6;% Michaelis constant for Glycine;
KI1311 = 4; % Inhibition constant for Serine

KM1312 = 0.075;% Michaelis constant for NAD;
KI1312 = 0.015;% Inhibition constant for NADH;    Since in the current program, we assume that P protein limit the 
                % rate of the overall glycin decarboxylase; the KI1312 and
                % KM1312 were not used.


KM1011 = 0.39;
KI1011 = 0.28;

KM1012 = 0.2;
KI1012 = 0.22;



% Reaction: 111: RUBP+O2<-->PGlycolate + PGA

PrKO = KO;
PrKC = KC;
PrKR = KR;

% Reaction: 112: PGlycolate-->Pi+Glycolate;
PrKM112 = KM112;
PrKI1122 = KI1122;
PrKI1121 = KI1121;

% Reaction 113  : Gla+ATP<-->ADP + PGA

PrKM1131 = KM1131;
PrKM1132 = KM1132;
PrKI113 = KI113;
PrKE113 = KE113;

% Reactoin 121; Glycolate +O2<-->H2O2+Glyoxylate

PrKM121 = KM121;

% Reaction 122  : Glyoxylate + Serine<--> Hydoxypyruvate + Glycine;

PrKM1221 = KM1221;
PrKM1222 = KM1222;
PrKI1221 = KI1221;
PrKE122 = KE122;

% Reaction 123: HydroxylPyruvate + NAD <--> NADH + Glycerate

PrKM123 = KM123;
PrKI123 = KI123;
PrKE123 = KE123;


% Reaction 124: Glyoxylate + Glu  <--> KG + Glycine;

PrKM1241 = KM1241;
PrKM1242 = KM1242;
PrKI124 = KI124;
PrKE124 = KE124;

% Reaction 131: LS2+Glycine <--> CO2+ AMDHL

PrKM1311 = KM1311;
PrKM1312 = KM1312;
PrKI1311 = KI1311;
PrKI1312 = KI1312;

% The consant for calculating the glycerate uptake.

PrKM1011 = KM1011;
PrKI1011 = KI1011;


% The constant for calculating the glycolate uptake

PrKM1012 = KM1012;
PrKI1012 = KI1012;


C = CO2_cond;
O = O2_cond;


if FIBF_PSPR_com ==1          % This term indicates whether the FIBF model should be connected with the PSPR model. 1 means connected. 
    if ATPActive == 0
        PrV111 = PrV111;
    end
end

if PR_PS_com ==1
%     global KM11;
%     global KM12;

    PrKC= PsKM11;
    PrKO= PsKM12;
end

if RUBISCOMETHOD ==2            % Using michelies and enzyme information
    if PR_PS_com ==1       % FOr the combined PS-PR model
        PrV111t = PrV111*Rubp/(Rubp+PrKR*V1Reg);
    else                    % For the PR model
        PrV111t = PrV111*Rubp/(Rubp+PrKR);
    end
    v111 = PrV111t * O/(O+PrKO*(1+C/PrKC));
    
elseif RUBISCOMETHOD==1
    v111 = PrV111 * O/(O+PrKO*(1+C/PrKC));
    if Rubp < RUBISCOTOTAL
        v111 = v111 * Rubp/RUBISCOTOTAL;
    end
end

v112 = PrV112 * Pgca /(Pgca + PrKM112*(1+Gca/PrKI1122)*(1+Pi/PrKI1121));

 if PR_PS_com ==1  ;           % For the combined PS-PR MODEL
     v113 = PrV113 * (ATP * Gcea-ADP * Pga/PrKE113)/((ATP + PrKM1131*(1 + Pga/PrKI113))*(Gcea + PrKM1132));
 else
     v113 = PrV113 * (ATP * Gcea-ADP * Pga/PrKE113)/((ATP + PrKM1131*(1 + 2.5/PrKI113))*(Gcea + PrKM1132));
 end

v121 = PrV121 * Gcac/(Gcac + PrKM121);

v122 = PrV122 * (Goac * Serc - Hprc * Glyc/PrKE122)/((Goac+PrKM1221)*(Serc + PrKM1222*(1+Glyc/PrKI1221)));
v123 = PrV123 * (Hprc *NADHc-Gcea * NADc/PrKE123)/(Hprc+PrKM123*(1+Hprc/PrKI123));  
v124 = PrV124 * (Goac * GLUc-KGc * Glyc/PrKE124)/((Goac + PrKM1241)*(GLUc + PrKM1242*(1+Glyc/PrKI124)));

v131 = PrV131 * Glyc/(Glyc + PrKM1311*(1+Serc/PrKI1311));
v2out = PrV2T * (Gca/(Gca + PrKM1012*(1+Gcea/PrKI1012)) - Gcac/(Gcac + PrKM1012*(1+Gceac/PrKI1012)));   % Competive inhibition
v1in = PrV1T *(Gceac/(Gceac + PrKM1011 *(1+ Gcac/PrKI1011))-Gcea /(Gcea + PrKM1011*(1+Gca/PrKI1011)));  % Competive inhibition

% Getting the information for output
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


CMr(17) = v111;
CMr(18)= v112;
CMr(19)= v113;
CMr(20)= v121;
CMr(21)= v122;
CMr(22)= v123;
CMr(23)= v124;
CMr(24)= v131;
CMr(25)= v1in;
CMr(26)= v2out;

%%% Calculate the rate equations

temp51 =Km512 * Km513 * ( 1 + GAPc/Km512 + DHAPc/Km513 + FBPc/Km511 + GAPc * DHAPc/(Km512 * Km513));
v51 = SUCSV51 * (GAPc * DHAPc - FBPc/KE51)/temp51;
Km521AP = Km521 * (1 + F26BPc/KI523);
temp52 = Km521AP * (1 + FBPc/Km521AP + Pic/KI522 + F6Pc/KI521 + Pic * F6Pc /(KI521 * KI522));
v52 = SUCSV52 * (FBPc - F6Pc  * Pic/KE52)/temp52;

temp55 = Km551 * Km552 * ( 1 + UTPc/Km551 + G1Pc/Km552 + UDPGc/Km553 + OPOPc/Km554 + UTPc * G1Pc/(Km551 * Km552) + UDPGc * OPOPc/(Km553 * Km554));
v55 = SUCSV55 * (UTPc * G1Pc - UDPGc * OPOPc/KE55)/temp55;

temp56 = (F6Pc + Km561 * (1 + FBPc/KI562))*(UDPGc + Km562 * ( 1 + UDPc/KI561)*(1+SUCP/KI563)*(1+Pic/KI564)*(1+SUC/KI565));
v56 = SUCSV56 * (F6Pc * UDPGc - SUCP * UDPc/KE56)/temp56;

temp57 = SUCP + Km571 * (1 + SUC/Ki572);
v57 = SUCSV57 * ( SUCP - SUC * Pic/KE57)/temp57;

temp58 = Km581 * (1 + F26BPc/Km581) * (1 + Pic/KI582) * (1+F6Pc/KI581);
v58 = SUCSV58 * F26BPc/temp58;

temp59 = (F6Pc + Km593 * ( 1 + F26BPc/Km592)*(1+DHAPc/KI592)) * (ATPc + Km591 * (1 + ADPc/KI591));
v59 = SUCSV59 * (ATPc * F6Pc - ADPc * F26BPc/KE59)/temp59;

temp60 = Km602 * Km603 * (1 + ATPc/Km602 + UDPc/Km603 + ATPc * UDPc/(Km602 * Km603) + ADPc/Km601 + UTPc/Km604 + ADPc * UTPc/(Km601 * Km604));
v60 = SUCSV60 * (ATPc * UDPc - ADPc * UTPc/KE60 )/temp60;

v62 = SUCSV62 * SUC/(SUC + Km621);

Vmatpf  = 0.25;
vatpf = Vmatpf * ADPc * Pic /((ADPc + 0.014)*(Pic + 0.3));

vdhap_in = 2 * Pic/(Pic + 2) ;
vgap_in  = 2 * Pic/(Pic + 2);

vpga_in = 0;
vpga_use = PGAc * 0.5/(PGAc + 1);



if (SUCS_TIME_N ==0)
    SUCS_TIME_N = 1;
end

if (t > SUCS_OLD_TIME)
    SUCS_TIME_N = SUCS_TIME_N + 1;
    SUCS_OLD_TIME = t;
end

SUCS_VEL	(	SUCS_TIME_N	,	1	)	=	t	;%	
SUCS_VEL	(	SUCS_TIME_N	,	2	)	=	v51	;%	DHAP+GAP --FBP

SUCS_VEL	(	SUCS_TIME_N	,	3	)	=	v52	;%	FBP --F6P + Pi
SUCS_VEL	(	SUCS_TIME_N	,	4	)	=	v55	;%	G1P+UTP --OPOP+UDPG 
SUCS_VEL	(	SUCS_TIME_N	,	5	)	=	v56	;%	UDPG+F6P--SUCP + UDP
SUCS_VEL	(	SUCS_TIME_N	,	6	)	=	v57	;%	SUCP--Pi + SUC
SUCS_VEL	(	SUCS_TIME_N	,	7	)	=	v58	;%	F26BP--F6P + Pi
SUCS_VEL	(	SUCS_TIME_N	,	8	)	=	v59	;%	F6P + ATP --ADP + F26BP
SUCS_VEL	(	SUCS_TIME_N	,	9	)	=	v60	;%	ATP+UDP --UTP + ADP
SUCS_VEL	(	SUCS_TIME_N	,	10	)	=	0	;%	POPO --2PO
SUCS_VEL	(	SUCS_TIME_N	,	11	)	=	v62	;%	SUC SINK 
SUCS_VEL	(	SUCS_TIME_N	,	12	)	=	vdhap_in	;%	DHAP IN
SUCS_VEL	(	SUCS_TIME_N	,	13	)	=	vgap_in	;%	GAP Export from chloroplast
SUCS_VEL	(	SUCS_TIME_N	,	14	)	=	vpga_in	;%	PGA export from chloroplast
SUCS_VEL	(	SUCS_TIME_N	,	15	)	=	vpga_use	;%	PGA utilisation in cytosol
SUCS_VEL	(	SUCS_TIME_N	,	16	)	=	vatpf	;%	ATP synthesis rate


SUCS_CON(SUCS_TIME_N,1) = t;
SUCS_CON(SUCS_TIME_N,2) = Pic;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assign table
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

CMr(27)	=	v51	;%	DHAP+GAP --FBP
CMr(28)=	v52	;%	FBP --F6P + Pi
CMr(29)=	v55	;%	G1P+UTP --OPOP+UDPG 
CMr(30)=	v56	;%	UDPG+F6P--SUCP + UDP
CMr(31)=	v57	;%	SUCP--Pi + SUC
CMr(32)=	v58	;%	F26BP--F6P + Pi
CMr(33)=	v59	;%	F6P + ATP --ADP + F26BP
CMr(34)=	v60	;%	ATP+UDP --UTP + ADP
CMr(35)=	0	;%	POPO --2PO
CMr(36)=	v62	;%	SUC SINK 
CMr(37)=	vdhap_in	;%	DHAP IN
CMr(38)=	vgap_in	;%	GAP Export from chloroplast
CMr(39)=	vpga_in	;%	PGA export from chloroplast
CMr(40)=	vpga_use	;%	PGA utilisation in cytosol
CMr(41)=	vatpf	;%	ATP synthesis rate	