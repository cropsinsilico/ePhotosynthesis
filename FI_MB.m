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


% FI_mb.m   This is the routine for calculation of the mass balance equations for the fluorescence induction model
% This routine is composed of two components;
% 1) The initialization of the rates that was transfered from the FI_Rate routine
% 2) The computation of the mass balance equations

function FI_mb = FI_Mb(t,FI_Con,FI_Param)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Calculate the rates first   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global GLight;
fini = Condition (t);
light = GLight;

FI_Param(1) = light;

FI_Vel = FI_Rate(t,FI_Con, FI_Param);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Get the rate of different reactions%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vA_d	=	FI_Vel	(	1	)	;	%	vA_d	The rate of heat dissipation from peripheral antenna
vA_f	=	FI_Vel	(	2	)	;	%	vA_f	The rate of fluorescence emission from peripheral antenna
vA_U	=	FI_Vel	(	3	)	;	%	vA_U	The rate of exciton transfer from peripheral antenna to core antenna in open reaction center
vU_A	=	FI_Vel	(	4	)	;	%	vU_A	The rate of exciton transfer from core antenna to perpheral antenna in open center
vU_f	=	FI_Vel	(	5	)	;	%	vU_f	The rate of fluorescence emission from core antenna
vU_d	=	FI_Vel	(	6	)	;	%	vU_d	The rate of heat dissipation from core antenna
v1	=	FI_Vel	(	7	)	    ;	%	v1	The rate of primary charge separation
v_r1	=	FI_Vel	(	8	)	;	%	v_r1	The rate of charge recombination
vS1_S2	=	FI_Vel	(	9	)	;	%	vS1_S2	The rate of transition from S1 to S2 
vS2_S3	=	FI_Vel	(	10	)	;	%	vS2_S3	The rate of transition from S2 to S3
vS3_S0	=	FI_Vel	(	11	)	;	%	vS3_S0	The rate of transition from S3 to S0
vS0_S1	=	FI_Vel	(	12	)	;	%	vS0_S1	The rate of transition from S0 to S1
vz_1	=	FI_Vel	(	13	)	;	%	vz_1	The rate of P680p reduction
v1z_1	=	FI_Vel	(	14	)	;	%	v1z_1	The rate of oxidation of S1T by P680pPheon
v2z_1	=	FI_Vel	(	15	)	;	%	v2z_1	The rate of oxidation of S2T  by P680pPheon
v3z_1	=	FI_Vel	(	16	)	;	%	v3z_1	The rate of oxidation of S3T  by P680pPheon
v0z_1	=	FI_Vel	(	17	)	;	%	v0z_1	The rate of oxidation of S0T  by P680pPheon
vz_2	=	FI_Vel	(	18	)	;	%	vz_2	The rate of P680pPheon reduction
v1z_2	=	FI_Vel	(	19	)	;	%	v1z_2	The rate of oxidation of S1T by P680pPheo
v2z_2	=	FI_Vel	(	20	)	;	%	v2z_2	The rate of oxidation of S2T  by P680pPheo
v3z_2	=	FI_Vel	(	21	)	;	%	v3z_2	The rate of oxidation of S3T  by P680pPheo
v0z_2	=	FI_Vel	(	22	)	;	%	v0z_2	The rate of oxidation of S0T  by P680pPheo
v1z	=	FI_Vel	(	23	)	;	%	v1z	
v2z	=	FI_Vel	(	24	)	;	%	v2z	
v3z	=	FI_Vel	(	25	)	;	%	v3z	
v0z	=	FI_Vel	(	26	)	;	%	v0z	
vAB1	=	FI_Vel	(	27	)	;	%	vAB1	The rate of electron transfer from QA- to QB
vBA1	=	FI_Vel	(	28	)	;	%	vBA1	The rate of electron transfer from QB- to QA
vAB2	=	FI_Vel	(	29	)	;	%	vAB2	The rate of electron transfer from QA- to QB-
vBA2	=	FI_Vel	(	30	)	;	%	vBA2	The rate of electron transfer from QB2- TO QA
v3	=	FI_Vel	(	31	)	;	%	v3	The rate of exchange of QAQBH2 with PQ
v_r3	=	FI_Vel	(	32	)	;	%	v_r3	The rate of exchange of QAQB with PQH2
v3_n	=	FI_Vel	(	33	)	;	%	v3_n	The rate of exchange of QAnQBH2 with PQ
v_r3_n	=	FI_Vel	(	34	)	;	%	v_r3_n	The rate of exchange of QAnQB with PQH2
v_pq_ox	=	FI_Vel	(	35	)	;	%	v_pq_ox	The rate of PQH2 oxidation
Ic 	=	FI_Vel	(	36	)	;	%	Ic	The incident light on the core antenna
Ia	=	FI_Vel	(	37	)	;	%	Ia	The incident light on the peripheral antenna
v2_1	=	FI_Vel	(	38	)	;	%	v2_1	The rate of P680pPheon oxidation
v2_2	=	FI_Vel	(	39	)	;	%	v2_1	The rate of P680pPheon oxidation
v2_00_1	=	FI_Vel	(	40	)	;	%	v2_00_1	The rate of reduction of QAQB by P680pPheon
v2_01_1	=	FI_Vel	(	41	)	;	%	v2_01_1	The rate of reduction of QAQBn by P680pPheon
v2_02_1	=	FI_Vel	(	42	)	;	%	v2_02_1	The rate of reduction of QAQB2n by P680pPheon
v2_00_2	=	FI_Vel	(	43	)	;	%	v2_00_2	The rate of reduction of QAQB by P680Pheon
v2_01_2	=	FI_Vel	(	44	)	;	%	v2_01_2	The rate of reduction of QAQBn by P680Pheon
v2_02_2	=	FI_Vel	(	45	)	;	%	v2_02_2	The rate of reduction of QAQB2n by P680Pheon
vr2_00_1	=	FI_Vel	(	46	)	;	%	vr2_00_1	The reverse reaction of The rate of reduction of QAQB by P680pPheon
vr2_01_1	=	FI_Vel	(	47	)	;	%	vr2_01_1	The reverse reaction of The rate of reduction of QAQBn by P680pPheon
vr2_02_1	=	FI_Vel	(	48	)	;	%	vr2_02_1	The reverse reaction of The rate of reduction of QAQB2n by P680pPheon
vr2_1	=	FI_Vel	(	49	)	;	%	vr2_1	
vr2_00_2	=	FI_Vel	(	50	)	;	%	vr2_00_2	The reverse reaction of The rate of reduction of QAQB by P680Pheon
vr2_01_2	=	FI_Vel	(	51	)	;	%	vr2_01_2	The reverse reaction of The rate of reduction of QAQBn by P680Pheon
vr2_02_2	=	FI_Vel	(	52	)	;	%	vr2_02_2	The reverse reaction of The rate of reduction of QAQB2n by P680Pheon
vr2_2	=	FI_Vel	(	53	)	;	%	vr2_2	
vP680qU = FI_Vel	(	54	)	;	%	vP680qU	
vP680qA = FI_Vel	(	55	)	;	%	vP680qA	
vU_P680 = FI_Vel (56);
vP680_d = FI_Vel (57);
vP680_f = FI_Vel (58);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get the mass balance equation %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This page defines the mass balance equation for the system under study											
% One problem need to be taken care of is the variables needed to transfer from FI_CalV to FI_mb											
% The Major Variables											
FI_mb				=	zeros(5,1)	;					
%FI_mb	(	1	)	=	Ia - vA_f - vA_d - vA_U + vU_A - vP680qA	;	%	A	
FI_mb	(	1	)	=	Ia - vA_f - vA_d - vA_U + vU_A;	%	A	6
FI_mb	(	2	)	=	Ic + vA_U - vU_A -vU_f - vU_d - v1 + v_r1 - vP680qU	;	%	U			
FI_mb	(	3	)	=	vU_P680 + v_r1 - v1 - vP680_d;	%	P680ePheo	QF add		
FI_mb	(	4	)	=	v1 - v_r1 - vz_1 - v2_1 + vr2_1	;	%	P680pPheon			
FI_mb	(	5	)	=	v2_1 - vr2_1 - vz_2	;	%	P680pPheo			
FI_mb	(	6	)	=	vz_1- v2_2 + vr2_2	;	%	P680Pheon			
FI_mb	(	7	)	=	vS1_S2 + vS2_S3 + vS3_S0 + vS0_S1 - vz_1 - vz_2	;	%	Yz			
FI_mb	(	8	)	=	vS0_S1 - v1z	;	%	S1T			
FI_mb	(	9	)	=	vS1_S2 - v2z	;	%	S2T			
FI_mb	(	10	)	=	vS2_S3 - v3z	;	%	S3T			
FI_mb	(	11	)	=	vS3_S0 - v0z	;	%	S0T			
FI_mb	(	12	)	=	v1z - vS1_S2	;	%	S1Tp			
FI_mb	(	13	)	=	v2z - vS2_S3	;	%	S2Tp		
FI_mb	(	14	)	=	v3z - vS3_S0	;	%	S3Tp		
FI_mb	(	15	)	=	v0z - vS0_S1	;	%	S0Tp		
FI_mb	(	16	)	=	v3 - v_r3 -v2_00_1 - v2_00_2 + vr2_00_1 + vr2_00_2	;	%	QAQB		
FI_mb	(	17	)	=	v2_00_1 + v2_00_2 - vr2_00_1 - vr2_00_2- vAB1 + vBA1 + v3_n - v_r3_n	;	%	QAnQB		
FI_mb	(	18	)	=	vAB1 -vBA1 - v2_01_1-v2_01_2 +vr2_01_1+vr2_01_2	;	%	QAQBn		
FI_mb	(	19	)	=	vBA2 - vAB2 + v2_01_1+v2_01_2 -vr2_01_1-vr2_01_2	;	%	QAnQBn		
FI_mb	(	20	)	=	vAB2 - vBA2  - v3  + v_r3  - v2_02_1-v2_02_2 +vr2_02_1+vr2_02_2	;	%	QAQB2n		
FI_mb	(	21	)	=	0-v3_n + v_r3_n + v2_02_1+v2_02_2 -vr2_02_1-vr2_02_2	;	%	QAnQB2n		
FI_mb	(	22	)	=	 v3 + v3_n - v_r3  - v_r3_n - v_pq_ox	;	%	PQn		
