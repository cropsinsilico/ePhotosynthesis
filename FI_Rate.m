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




function FI_Vel = FI_Rate(t,FI_Con, FI_Param)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 1 Get the rate constant and the initial concentrations % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global FI_RC;

% The rate constant used in the model												
kA_d	=	FI_RC	(	1	)				;	%	The rate constant of heat dissipation from peripheral antenna	Lazar (1999), 0.25~1 *10^(9)			
kA_f	=	FI_RC	(	2	)				;	%	The rate constant of fluorescence emission from peripheral antenna	Lazar 1999, with a lifetime of 5 ns at closed reaction center			
kA_U	=	FI_RC	(	3	)				;	%	The rate constant of exciton transfer from periphral antenna to core antenna	Reference needed, a guess			
kU_A	=	FI_RC	(	4	)				;	%	The rate constant of exciton transfer from core antenna to peripheral antenna	Reference needed, a guess			
kU_d	=	FI_RC	(	5	)				;	%	The rate constant of  heat emission from core antenna				
kU_f	=	FI_RC	(	6	)				;	%	The rate constant of fluorescence emission from core antenna				
k1	    =	FI_RC	(	7	)		        ;	%	The rate constant of primary charge separation for open reaction center				
k_r1	=	FI_RC	(	8	)				;	%	The rate constant of charge recombination for open reactoin center				
kz	=	FI_RC	(	9	)				    ;	%	The rate constant of the Tyrosine oxidation	Lazar (1999); 3.8~50 * 10^6 			
k12	=	FI_RC	(	10	)				    ;	%	The rate constant of the S1 to S2 transition	Lazar (1999); 0.667~33.3 * 10^3			
k23	=	FI_RC	(	11	)				    ;	%	The rate constant of the S2 to S3 transition	Lazar (1999); 0.667~33.3 * 10^3			
k30	=	FI_RC	(	12	)				    ;	%	The rate constant of the S3 to S0 transition	Lazar (1999); 0.667~33.3 * 10^3			
k01	=	FI_RC	(	13	)				    ;	%	The rate constant of the S0 to S1 transition	Lazar (1999); 0.667~33.3 * 10^3			
k2	=	FI_RC	(	14	)				    ;	%	The rate constant of the QA reduction by Pheo-	Lazar (1999); 2~2.3 * 10^9			
kAB1	=	FI_RC	(	15	)				;	%	The rate constant of QAQB-->QAQB-	Lazar (1999); 2.5~5 * 10^3			
kBA1	=	FI_RC	(	16	)				;	%	The rate constant of the QAQB- -->QAQB	Lazar (1999)			
kAB2	=	FI_RC	(	17	)				;	%	The rate constant of the QAQB- --> QAQB2-	Lazar (1999); 1.25~3.33 * 10^3			
kBA2	=	FI_RC	(	18	)				;	%	The rate constant of the QAQB2- --> QAQB- 	Lazar (1999), or same as kAB2 depend on the equilibium constant			
k3	=	FI_RC	(	19	)				    ;	%	The rate constant of the exchange of PQ and QBH2	Lazar (1999),0.12~1 for the fast PQ pool,  or 3~8 for the slow recycling PQ pool			
k_r3	=	FI_RC	(	20	)				;	%	The rate constant of the exchange of QB and PQH2	Lazar (1999), since the equilibrium constant is 1 (205 in Lazar, 1999) 			
k_pq_oxy	=	FI_RC	(	21	)				;	%	The rate constant of the PQH2 oxidation	Lazar (1999),50~500			
															
               
A	=	FI_Con	(	1	)			;	% 	The concentration of excitons in the peripheral antenna
U	=	FI_Con	(	2	)			;	%	The concentration fo excitons in the core antenna
P680ePheo   =   FI_Con  (   3   );          %QF add
P680pPheon	=	FI_Con	(	4	)			;	%	The concentration for the P680+ Pheo-
P680pPheo	=	FI_Con	(	5	)			;	% 	The concentration of P680+ Pheo
P680Pheon	=	FI_Con	(	6	)			;	%	The concentration of P680Pheo-
Yz	=	FI_Con	(	7	)			;	%	The concentration of reduced tyrosine
S1T	=	FI_Con	(	8	)			;	% 	The concentration of S1 in combination with reduced tyrosine 
S2T	=	FI_Con	(	9	)			;	%	The concentration of S2 in combination with reduced tyrosine
S3T	=	FI_Con	(	10	)			;	% 	The concentration of S3 in combination with reduced tyrosine
S0T	=	FI_Con	(	11	)			;	%	The concentration of S0 in combination with reduced tyrosine 
S1Tp	=	FI_Con	(	12	)		;	% 	The concentration of S1 in combination with oxidized tyrosine
S2Tp	=	FI_Con	(	13	)		;	% 	The concentration of S2 in combination with oxidized tyrosine
S3Tp	=	FI_Con	(	14	)		;	%	The concentration of S3 in combination with oxidized tyrosine
S0Tp	=	FI_Con	(	15	)		;	% 	The concentration of S0 in combination with oxidized tyrosine
QAQB	=	FI_Con	(	16	)		;	% 	The concentration of [QAQB]
QAnQB	=	FI_Con	(	17	)		;	% 	The concentration of [QA-QB];
QAQBn	=	FI_Con	(	18	)		;	%	The concentration of [QAQB-]
QAnQBn	=	FI_Con	(	19	)		;	% 	The concentration of [QA-QB-];
QAQB2n	=	FI_Con	(	20	)		;	%	The concentration of [QAQB2-]
QAnQB2n	=	FI_Con	(	21	)		;	% 	The concentration of [QA-QB2-];
PQn	=	FI_Con	(	22	)			;	%	The concentration of reduced PQ, i.e. PQH2;
										
global FI_Pool;
PQT = FI_Pool (2);  % The total concentraion of PQH2 and PQ;
global FIBF_AUX;
PQa = FIBF_AUX(2);
PQ = PQT - PQn - PQa; 

global BF_FI_com; 
if BF_FI_com ==1
    global FIBF2FI_PQ;     
    PQ = FIBF2FI_PQ;
end
 

P680PheoT	=	1	;

global ChlT2;
global ChlT;
global ChlPSI;

% ChlT	=	70	;	%	The amount of chl molecules in U in one meter square leaf area; unit: micor mole per meter square								
% ChlT2  = 290;       % The total amoutn of chlorophyll in one PSII unit, inclusing both U and A.

P680Pheo = P680PheoT - P680pPheo - P680Pheon - P680pPheon - P680ePheo;  %QF add  '- P680ePheo'

%P680ePheo	=	U /(ChlT) * P680Pheo	;	%	The amount of excitons on P680 molecules; the difference in the energy of P680 and antenna chlorophyll is not incorporated								

n	=	FI_Param(2)	;	%	n	The ratio of the number of PSI to PSII
It	=	FI_Param(1)	;	%	It	The total incident light intensity

% rate: U -> U*
Ic 	=	It * ChlT/(ChlT2 + ChlPSI)	;	%	Ic	The incident light on the core antenna; ChlT is defined in upper lines as the total amount of Chl in one U. 
% rate: A -> A*
Ia	=	It * (ChlT2 - ChlT)/(ChlT2 + ChlPSI)	;	%	Ia	The incident light on the peripheral antenna			

% Ic 	=	It * ChlT/(ChlT2 + 200*n)	;	%	Ic	The incident light on the core antenna; ChlT is defined in upper lines as the total amount of Chl in one U. 
% Ia	=	It * (ChlT2 - ChlT)/(ChlT2 + 200*n)	;	%	Ia	The incident light on the peripheral antenna			

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate the rate of different reactions %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

q = (QAQB+QAQBn+QAQB2n)/(QAQB + QAQBn+ QAQB2n + QAnQB + QAnQBn+QAnQB2n);

vA_d	=	A * kA_d	;	%	vA_d	The rate of heat dissipation from peripheral antenna
vA_f	=	A * kA_f	;	%	vA_f	The rate of fluorescence emission from peripheral antenna
vA_U	=	A *  kA_U	;	%	vA_U	The rate of exciton transfer from peripheral antenna to core antenna in open reaction center
vU_A	=	U * kU_A	    ;	%	vU_A	The rate of exciton transfer from core antenna to perpheral antenna in open center
vU_f	=	U * kU_f	;	%	vU_f	The rate of fluorescence emission from core antenna
vU_d	=	U * kU_d *(1-q)	;	%	vU_d	The rate of heat dissipation from core antenna

P = 1; 

% energy flow: 
vU_P680 = Ic + vA_U - vU_A -vU_f - vU_d; %QF add  , total energy coming to P680 and == the rate of P680 -> P680*, except f and d, energy is transport to P680 reaction center
vP680_d = P680ePheo * kU_d *(1-q);
v1	=	P680ePheo * k1 * q + P680ePheo * P * (1-q) * k1/6.2 + P680ePheo * (1-P)* (1-q) * k1;	%	v1	The rate of primary charge separation
v_r1	=	P680pPheon * k_r1*q + P680pPheon *(1-q)* k_r1*3 ;	%	v_r1	The rate of charge recombination
vP680_f = vU_P680 - (v1 - v_r1) - vP680_d;

vS1_S2	=	S1Tp * k12	;	%	vS1_S2	The rate of transition from S1 to S2 
vS2_S3	=	S2Tp * k23	;	%	vS2_S3	The rate of transition from S2 to S3
vS3_S0	=	S3Tp * k30	;	%	vS3_S0	The rate of transition from S3 to S0
vS0_S1	=	S0Tp * k01	;	%	vS0_S1	The rate of transition from S0 to S1

Coeff	=	P680pPheon/P680PheoT	;			
v1z_1	=	S1T * kz  * Coeff	;	%	v1z_1	The rate of oxidation of S1T by P680pPheon
v2z_1	=	S2T* kz* Coeff	;	%	v2z_1	The rate of oxidation of S2T  by P680pPheon
v3z_1	=	S3T* kz* Coeff	;	%	v3z_1	The rate of oxidation of S3T  by P680pPheon
v0z_1	=	S0T * kz* Coeff	;	%	v0z_1	The rate of oxidation of S0T  by P680pPheon
vz_1	=	v1z_1 + v2z_1 + v3z_1 + v0z_1	;	%	vz_1	The rate of P680pPheon reduction

Coeff	=	P680pPheo/P680PheoT	;		

v1z_2	=	S1T * kz  * Coeff	;	%	v1z_2	The rate of oxidation of S1T by P680pPheo
v2z_2	=	S2T* kz* Coeff	;	%	v2z_2	The rate of oxidation of S2T  by P680pPheo
v3z_2	=	S3T* kz* Coeff	;	%	v3z_2	The rate of oxidation of S3T  by P680pPheo
v0z_2	=	S0T * kz* Coeff	;	%	v0z_2	The rate of oxidation of S0T  by P680pPheo
vz_2	=	v1z_2 + v2z_2 + v3z_2 + v0z_2	;	%	vz_2	The rate of P680pPheo reduction 

v1z	=	v1z_1 + v1z_2	;			
v2z	=	v2z_1 + v2z_2	;			
v3z	=	v3z_1 + v3z_2	;			
v0z	=	v0z_1 + v0z_2	;			

vAB1	=	QAnQB * kAB1	;	%	vAB1	The rate of electron transfer from QA- to QB
vBA1	=	QAQBn * kBA1	;	%	vBA1	The rate of electron transfer from QB- to QA
vAB2	=	QAnQBn * kAB2	;	%	vAB2	The rate of electron transfer from QA- to QB-
vBA2	=	QAQB2n * kBA2	;	%	vBA2	The rate of electron transfer from QB2- TO QA
v3	=	QAQB2n * PQ * k3/PQT	;	%	v3	The rate of exchange of QAQBH2 with PQ
v_r3	=	QAQB * PQn * k_r3/PQT	;	%	v_r3	The rate of exchange of QAQB with PQH2

v3_n	=	QAnQB2n * PQ * k3/PQT	;	%	v3_n	The rate of exchange of QAnQBH2 with PQ
v_r3_n	=	QAnQB * PQn * k_r3/PQT	;	%	v_r3_n	The rate of exchange of QAnQB with PQH2

v_pq_ox	=	PQn * k_pq_oxy	;	%	v_pq_ox	The rate of PQH2 oxidation

v2_1	=	P680pPheon * k2 *q	;	%	v2_1	The rate of P680pPheon oxidation			
v2_2	=	P680Pheon * k2 *q	;	%	v2_1	The rate of P680pPheon oxidation			
a	=	QAQB/(QAQB+QAQBn + QAQB2n)	;	%	a				
b	=	QAQBn/(QAQB+QAQBn + QAQB2n)	;	%	b				
c	=	QAQB2n/(QAQB+QAQBn + QAQB2n)	;	%	c			

v2_00_1	=	v2_1 * a	;	%	v2_00_1	The rate of reduction of QAQB by P680pPheon			
v2_01_1	=	v2_1 * b	;	%	v2_01_1	The rate of reduction of QAQBn by P680pPheon			
v2_02_1	=	v2_1 * c	;	%	v2_02_1	The rate of reduction of QAQB2n by P680pPheon			
									
v2_00_2	=	v2_2 * a	;	%	v2_00_2	The rate of reduction of QAQB by P680Pheon			
v2_01_2	=	v2_2 * b	;	%	v2_01_2	The rate of reduction of QAQBn by P680Pheon			
v2_02_2	=	v2_2 * c	;	%	v2_02_2	The rate of reduction of QAQB2n by P680Pheon			

KE	=	1000000	;			
Coeff1	=	P680pPheo/P680PheoT	;	        %	Coeff1	
vr2_00_1	=	QAnQB * k2/KE *Coeff1	;	%	vr2_00_1	The reverse reaction of The rate of reduction of QAQB by P680pPheon
vr2_01_1	=	QAnQBn * k2/KE *Coeff1	;	%	vr2_01_1	The reverse reaction of The rate of reduction of QAQBn by P680pPheon
vr2_02_1	=	QAnQB2n * k2/KE *Coeff1	;	%	vr2_02_1	The reverse reaction of The rate of reduction of QAQB2n by P680pPheon
vr2_1	=	vr2_00_1 + vr2_01_1 + vr2_02_1	;	%	vr2_1	

Coeff2	=	P680Pheo/P680PheoT	;	        %	Coeff2	
vr2_00_2	=	QAnQB * k2/KE *Coeff2	;	%	vr2_00_2	The reverse reaction of The rate of reduction of QAQB by P680Pheon
vr2_01_2	=	QAnQBn * k2/KE *Coeff2	;	%	vr2_01_2	The reverse reaction of The rate of reduction of QAQBn by P680Pheon
vr2_02_2	=	QAnQB2n * k2/KE *Coeff2	;	%	vr2_02_2	The reverse reaction of The rate of reduction of QAQB2n by P680Pheon
vr2_2	=	vr2_00_2 + vr2_01_2 + vr2_02_2	;	%	vr2_2	

vP680qU = 10^9 * U * (P680pPheo + P680pPheon) + U * 0.15 * (kU_f + kU_d) * PQ/PQT;
vP680qA = 10^9 * A * (P680pPheo + P680pPheon) +  A * 0.15 * (kA_f + kA_d) * PQ/PQT;

%%%%%%%%%%%%%%%%%%%%
%%%  FOR TESITNG  %%
f = vA_f + vU_f;
%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Part V Output of Velocity for plot          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global FI_OLD_TIME;
global FI_TIME_N;
global FI_VEL;
global FI_CON;

if (FI_TIME_N ==0)
    FI_TIME_N = 1;
end

if (t > FI_OLD_TIME)
    FI_TIME_N = FI_TIME_N + 1;
    FI_OLD_TIME = t;
end

FI_VEL	(	FI_TIME_N	,	1	)	=	t;
FI_VEL	(	FI_TIME_N	,	2	)	=	vA_d;
FI_VEL	(	FI_TIME_N	,	3	)	=	vA_f;
FI_VEL	(	FI_TIME_N	,	4	)	=	vA_U;
FI_VEL	(	FI_TIME_N	,	5	)	=	vU_A;
FI_VEL	(	FI_TIME_N	,	6	)	=	vU_f;
FI_VEL	(	FI_TIME_N	,	7	)	=	vU_d;
FI_VEL	(	FI_TIME_N	,	8	)	=	v1;
FI_VEL	(	FI_TIME_N	,	9	)	=	v_r1;
FI_VEL	(	FI_TIME_N	,	10	)	=	vS1_S2;
FI_VEL	(	FI_TIME_N	,	11	)	=	vS2_S3;
FI_VEL	(	FI_TIME_N	,	12	)	=	vS3_S0;
FI_VEL	(	FI_TIME_N	,	13	)	=	vS0_S1;
FI_VEL	(	FI_TIME_N	,	14	)	=	vz_1;
FI_VEL	(	FI_TIME_N	,	15	)	=	v1z_1;
FI_VEL	(	FI_TIME_N	,	16	)	=	v2z_1;
FI_VEL	(	FI_TIME_N	,	17	)	=	v3z_1;
FI_VEL	(	FI_TIME_N	,	18	)	=	v0z_1;
FI_VEL	(	FI_TIME_N	,	19	)	=	vz_2;
FI_VEL	(	FI_TIME_N	,	20	)	=	v1z_2;
FI_VEL	(	FI_TIME_N	,	21	)	=	v2z_2;
FI_VEL	(	FI_TIME_N	,	22	)	=	v3z_2;
FI_VEL	(	FI_TIME_N	,	23	)	=	v0z_2;
FI_VEL	(	FI_TIME_N	,	24	)	=	v1z;
FI_VEL	(	FI_TIME_N	,	25	)	=	v2z;
FI_VEL	(	FI_TIME_N	,	26	)	=	v3z;
FI_VEL	(	FI_TIME_N	,	27	)	=	v0z;
FI_VEL	(	FI_TIME_N	,	28	)	=	vAB1;
FI_VEL	(	FI_TIME_N	,	29	)	=	vBA1;
FI_VEL	(	FI_TIME_N	,	30	)	=	vAB2;
FI_VEL	(	FI_TIME_N	,	31	)	=	vBA2;
FI_VEL	(	FI_TIME_N	,	32	)	=	v3;
FI_VEL	(	FI_TIME_N	,	33	)	=	v_r3;
FI_VEL	(	FI_TIME_N	,	34	)	=	v3_n;
FI_VEL	(	FI_TIME_N	,	35	)	=	v_r3_n;
FI_VEL	(	FI_TIME_N	,	36	)	=	v_pq_ox;
FI_VEL	(	FI_TIME_N	,	37	)	=	Ic ;
FI_VEL	(	FI_TIME_N	,	38	)	=	Ia;
FI_VEL	(	FI_TIME_N	,	39	)	=	v2_1;
FI_VEL	(	FI_TIME_N	,	40	)	=	v2_2;
FI_VEL	(	FI_TIME_N	,	41	)	=	v2_00_1;
FI_VEL	(	FI_TIME_N	,	42	)	=	v2_01_1;
FI_VEL	(	FI_TIME_N	,	43	)	=	v2_02_1;
FI_VEL	(	FI_TIME_N	,	44	)	=	v2_00_2;
FI_VEL	(	FI_TIME_N	,	45	)	=	v2_01_2;
FI_VEL	(	FI_TIME_N	,	46	)	=	v2_02_2;
FI_VEL	(	FI_TIME_N	,	47	)	=	vr2_00_1;
FI_VEL	(	FI_TIME_N	,	48	)	=	vr2_01_1;
FI_VEL	(	FI_TIME_N	,	49	)	=	vr2_02_1;
FI_VEL	(	FI_TIME_N	,	50	)	=	vr2_1;
FI_VEL	(	FI_TIME_N	,	51	)	=	vr2_00_2;
FI_VEL	(	FI_TIME_N	,	52	)	=	vr2_01_2;
FI_VEL	(	FI_TIME_N	,	53	)	=	vr2_02_2;
FI_VEL	(	FI_TIME_N	,	54	)	=	vr2_2;
FI_VEL	(	FI_TIME_N	,	55	)	=	vP680qU;
FI_VEL	(	FI_TIME_N	,	56	)	=	vP680qA;
FI_VEL	(	FI_TIME_N	,	57	)	=	vU_P680;
FI_VEL	(	FI_TIME_N	,	58	)	=	vP680_d;
FI_VEL	(	FI_TIME_N	,	59	)	=	vP680_f;
							

FI_CON(FI_TIME_N,1) = t;
FI_CON(FI_TIME_N,2) = f;
FI_CON(FI_TIME_N,3) = 1-q;
FI_CON(FI_TIME_N,4) = vA_d + vU_d;
FI_CON(FI_TIME_N,5) = vS3_S0;         
if It ==0    
    fPSII = 0; 
else
    It2 = It * 27/47;
    fPSII = (It2 - f - vA_d - vU_d)/It2;     
end
    
FI_CON(FI_TIME_N,6) = fPSII;            


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assign table
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global FI_Vel;

FI_Vel	(	1	)	=	vA_d	;	%	vA_d	The rate of heat dissipation from peripheral antenna
FI_Vel	(	2	)	=	vA_f	;	%	vA_f	The rate of fluorescence emission from peripheral antenna
FI_Vel	(	3	)	=	vA_U	;	%	vA_U	The rate of exciton transfer from peripheral antenna to core antenna in open reaction center
FI_Vel	(	4	)	=	vU_A	;	%	vU_A	The rate of exciton transfer from core antenna to perpheral antenna in open center
FI_Vel	(	5	)	=	vU_f	;	%	vU_f	The rate of fluorescence emission from core antenna
FI_Vel	(	6	)	=	vU_d	;	%	vU_d	The rate of heat dissipation from core antenna
FI_Vel	(	7	)	=	v1	;	%	v1	The rate of primary charge separation
FI_Vel	(	8	)	=	v_r1	;	%	v_r1	The rate of charge recombination
FI_Vel	(	9	)	=	vS1_S2	;	%	vS1_S2	The rate of transition from S1 to S2 
FI_Vel	(	10	)	=	vS2_S3	;	%	vS2_S3	The rate of transition from S2 to S3
FI_Vel	(	11	)	=	vS3_S0	;	%	vS3_S0	The rate of transition from S3 to S0
FI_Vel	(	12	)	=	vS0_S1	;	%	vS0_S1	The rate of transition from S0 to S1
FI_Vel	(	13	)	=	vz_1	;	%	vz_1	The rate of P680p reduction
FI_Vel	(	14	)	=	v1z_1	;	%	v1z_1	The rate of oxidation of S1T by P680pPheon
FI_Vel	(	15	)	=	v2z_1	;	%	v2z_1	The rate of oxidation of S2T  by P680pPheon
FI_Vel	(	16	)	=	v3z_1	;	%	v3z_1	The rate of oxidation of S3T  by P680pPheon
FI_Vel	(	17	)	=	v0z_1	;	%	v0z_1	The rate of oxidation of S0T  by P680pPheon
FI_Vel	(	18	)	=	vz_2	;	%	vz_2	The rate of P680pPheon reduction
FI_Vel	(	19	)	=	v1z_2	;	%	v1z_2	The rate of oxidation of S1T by P680pPheo
FI_Vel	(	20	)	=	v2z_2	;	%	v2z_2	The rate of oxidation of S2T  by P680pPheo
FI_Vel	(	21	)	=	v3z_2	;	%	v3z_2	The rate of oxidation of S3T  by P680pPheo
FI_Vel	(	22	)	=	v0z_2	;	%	v0z_2	The rate of oxidation of S0T  by P680pPheo

FI_Vel	(	23	)	=	v1z	;	%	v1z	
FI_Vel	(	24	)	=	v2z	;	%	v2z	
FI_Vel	(	25	)	=	v3z	;	%	v3z	
FI_Vel	(	26	)	=	v0z	;	%	v0z	

FI_Vel	(	27	)	=	vAB1	;	%	vAB1	The rate of electron transfer from QA- to QB
FI_Vel	(	28	)	=	vBA1	;	%	vBA1	The rate of electron transfer from QB- to QA
FI_Vel	(	29	)	=	vAB2	;	%	vAB2	The rate of electron transfer from QA- to QB-
FI_Vel	(	30	)	=	vBA2	;	%	vBA2	The rate of electron transfer from QB2- TO QA
FI_Vel	(	31	)	=	v3	;	%	v3	The rate of exchange of QAQBH2 with PQ

FI_Vel	(	32	)	=	v_r3	;	%	v_r3	The rate of exchange of QAQB with PQH2
FI_Vel	(	33	)	=	v3_n	;	%	v3_n	The rate of exchange of QAnQBH2 with PQ
FI_Vel	(	34	)	=	v_r3_n	;	%	v_r3_n	The rate of exchange of QAnQB with PQH2
FI_Vel	(	35	)	=	v_pq_ox	;	%	v_pq_ox	The rate of PQH2 oxidation
FI_Vel	(	36	)	=	Ic 	;	    %	Ic	The incident light on the core antenna

FI_Vel	(	37	)	=	Ia	;	    %	Ia	The incident light on the peripheral antenna
FI_Vel	(	38	)	=	v2_1	;	%	v2_1	The rate of P680pPheon oxidation
FI_Vel	(	39	)	=	v2_2	;	%	v2_1	The rate of P680pPheon oxidation
FI_Vel	(	40	)	=	v2_00_1	;	%	v2_00_1	The rate of reduction of QAQB by P680pPheon
FI_Vel	(	41	)	=	v2_01_1	;	%	v2_01_1	The rate of reduction of QAQBn by P680pPheon
FI_Vel	(	42	)	=	v2_02_1	;	%	v2_02_1	The rate of reduction of QAQB2n by P680pPheon

FI_Vel	(	43	)	=	v2_00_2	;	%	v2_00_2	The rate of reduction of QAQB by P680Pheon
FI_Vel	(	44	)	=	v2_01_2	;	%	v2_01_2	The rate of reduction of QAQBn by P680Pheon
FI_Vel	(	45	)	=	v2_02_2	;	%	v2_02_2	The rate of reduction of QAQB2n by P680Pheon
FI_Vel	(	46	)	=	vr2_00_1	;	%	vr2_00_1	The reverse reaction of The rate of reduction of QAQB by P680pPheon
FI_Vel	(	47	)	=	vr2_01_1	;	%	vr2_01_1	The reverse reaction of The rate of reduction of QAQBn by P680pPheon
FI_Vel	(	48	)	=	vr2_02_1	;	%	vr2_02_1	The reverse reaction of The rate of reduction of QAQB2n by P680pPheon
FI_Vel	(	49	)	=	vr2_1	;	%	vr2_1	

FI_Vel	(	50	)	=	vr2_00_2	;	%	vr2_00_2	The reverse reaction of The rate of reduction of QAQB by P680Pheon
FI_Vel	(	51	)	=	vr2_01_2	;	%	vr2_01_2	The reverse reaction of The rate of reduction of QAQBn by P680Pheon
FI_Vel	(	52	)	=	vr2_02_2	;	%	vr2_02_2	The reverse reaction of The rate of reduction of QAQB2n by P680Pheon

FI_Vel	(	53	)	=	vr2_2	;	%	vr2_2	
FI_Vel	(	54	)	=	vP680qU	;	%	vr2_2	
FI_Vel	(	55	)	=	vP680qA	;	%	vr2_2	
FI_Vel  (   56  )   =   vU_P680;
FI_Vel  (   57  )   =   vP680_d;
FI_Vel  (   58  )   =   vP680_f;

global export_mod_enabled;
if export_mod_enabled == 1
    export_mod_dict(t, "FI", FI_Con, "");
end
