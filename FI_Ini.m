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


% FI_Init.m     This is the routine that initialize the parameters, initial conditions for simulation of fluorescence induction curve.
% The following information is initialized sequentially 1) Rate constants; 2) Initial concentration ( or conditions); 3) THe maximum 
% concentration of components of photosystems.

function FI_Con = FI_Ini(begin)
global cpsii;
%cpsii=1.0169;% 1.0129;%WY201803
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initilization of the rate constant %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The rate constant used in the model									
% The rate constant used in the model									
% The rate constant used in the model									
						% Reference			
global FIVmax;                        
% kA_d	=	2*10^8	;	%	The rate constant of heat dissipation from peripheral antenna	Lazar (1999), 0.25~1 *10^(9)			
% kA_f	=	6.3 *10^6 *0.2	;	%	The rate constant of fluorescence emission from peripheral antenna; based on the fact that the natural lifetime of chlorophyll a is 15.9 ns; This value is correct and should not be changed by any.
% kA_U	=	10^10	;	%	The rate constant of exciton transfer from periphral antenna to core antenna	Reference needed, a guess			
% kU_A	=	10^10	;	%	The rate constant of exciton transfer from core antenna to peripheral antenna	Reference needed, a guess			
% kU_d	=	2*10^8	;	%	The rate constant of the heat dissipation from core antenna; Laverage and Trissl, 1995				
% kU_f	=	6.3 *10^6 *0.2		;	% 	The rate constant of the fluorescence emission from the core antenna; Laverage and Trissl 1995				
% 
% k1	=	2.5 * 10^11*cpsii;	%WY201803	The rate constant of primary charge separation for open reaction center				
% k_r1	=	3*10^8	;	%	The rate constant of charge recombination for open reactoin center; Lazar, 1999;				
% kz	=	5*10^6	;	%	The rate constant of the Tyrosine oxidation	Lazar (1999); 3.8~50 * 10^6 			
% k12	=	30000	;	%	The rate constant of the S1 to S2 transition	Lazar (1999); 0.667~33.3 * 10^3			
% k23	=	10000	;	%	The rate constant of the S2 to S3 transition	Lazar (1999); 0.667~33.3 * 10^3			
% k30	=	3000	;	%	The rate constant of the S3 to S0 transition	Lazar (1999); 0.667~33.3 * 10^3			
% k01	=	500	;	%	The rate constant of the S0 to S1 transition	Lazar (1999); 0.667~33.3 * 10^3			
% 
% k2	    =	2*10^9	;	%	The rate constant of the QA reduction by Pheo-	Lazar (1999); 2~2.3 * 10^9			
% kAB1	=	2500	;	%	The rate constant of QAQB-->QAQB-	Lazar (1999); 2.5~5 * 10^3			
% kBA1	=	200	;	%	The rate constant of the QAQB- -->QAQB	Lazar (1999)			
% kAB2	=	3300	;	%	The rate constant of the QAQB- --> QAQB2-	Lazar (1999); 1.25~3.33 * 10^3			
% kBA2	=	250	;	%	The rate constant of the QAQB2- --> QAQB- 	Lazar (1999), or same as kAB2 depend on the equilibium constant			
% k3	    =	800	;	%	The rate constant of the exchange of PQ and QBH2	Lazar (1999),0.12~1 for the fast PQ pool,  or 3~8 for the slow recycling PQ pool			
% k_r3	=	80	;	%	The rate constant of the exchange of QB and PQH2	Lazar (1999), since the equilibrium constant is 1 (205 in Lazar, 1999) 			
% k_pq_oxy=	500	;	%	The rate constant of the PQH2 oxidation	Lazar (1999),50~500			
kA_d=FIVmax(1);
kA_f=FIVmax(2);
kA_U=FIVmax(3);
kU_A=FIVmax(4);
kU_d=FIVmax(5);
kU_f=FIVmax(6);
k1=FIVmax(7)*cpsii;
k_r1=FIVmax(8);
kz=FIVmax(9);
k12=FIVmax(10);
k23=FIVmax(11);
k30=FIVmax(12);
k01=FIVmax(13);
k2=FIVmax(14);
kAB1=FIVmax(15);
kBA1=FIVmax(16);
kAB2=FIVmax(17);
kBA2=FIVmax(18);
k3=FIVmax(19);
k_r3=FIVmax(20);
k_pq_oxy=FIVmax(21);

% The rate constant used in the model											
global FI_RC;
FI_RC= zeros(5,1);
			% The rate constant used in the model									
			% The rate constant used in the model												
		FI_RC	(	1	)	=	kA_d		;	%	The rate constant of heat dissipation from peripheral antenna	Lazar (1999), 0.25~1 *10^(9)			
		FI_RC	(	2	)	=	kA_f		;	%	The rate constant of fluorescence emission from peripheral antenna	Lazar 1999, with a lifetime of 5 ns at closed reaction center			
		FI_RC	(	3	)	=	kA_U		;	%	The rate constant of exciton transfer from periphral antenna to core antenna	Reference needed, a guess			
		FI_RC	(	4	)	=	kU_A		;	%	The rate constant of exciton transfer from core antenna to peripheral antenna	Reference needed, a guess			
		FI_RC	(	5	)	=	kU_d		;	%	The rate constant of  heat emission from core antenna				
		FI_RC	(	6	)	=	kU_f		;	%	The rate constant of fluorescence emission from core antenna				
		FI_RC	(	7	)	=	k1		;	%	The rate constant of primary charge separation for open reaction center				
		FI_RC	(	8	)	=	k_r1		;	%	The rate constant of charge recombination for open reactoin center				
		FI_RC	(	9	)	=	kz		;	%	The rate constant of the Tyrosine oxidation	Lazar (1999); 3.8~50 * 10^6 			
		FI_RC	(	10	)	=	k12		;	%	The rate constant of the S1 to S2 transition	Lazar (1999); 0.667~33.3 * 10^3			
		FI_RC	(	11	)	=	k23		;	%	The rate constant of the S2 to S3 transition	Lazar (1999); 0.667~33.3 * 10^3			
		FI_RC	(	12	)	=	k30		;	%	The rate constant of the S3 to S0 transition	Lazar (1999); 0.667~33.3 * 10^3			
		FI_RC	(	13	)	=	k01		;	%	The rate constant of the S0 to S1 transition	Lazar (1999); 0.667~33.3 * 10^3			
		FI_RC	(	14	)	=	k2		;	%	The rate constant of the QA reduction by Pheo-	Lazar (1999); 2~2.3 * 10^9			
		FI_RC	(	15	)	=	kAB1		;	%	The rate constant of QAQB-->QAQB-	Lazar (1999); 2.5~5 * 10^3			
		FI_RC	(	16	)	=	kBA1		;	%	The rate constant of the QAQB- -->QAQB	Lazar (1999)			
		FI_RC	(	17	)	=	kAB2		;	%	The rate constant of the QAQB- --> QAQB2-	Lazar (1999); 1.25~3.33 * 10^3			
		FI_RC	(	18	)	=	kBA2		;	%	The rate constant of the QAQB2- --> QAQB- 	Lazar (1999), or same as kAB2 depend on the equilibium constant			
		FI_RC	(	19	)	=	k3		;	%	The rate constant of the exchange of PQ and QBH2	Lazar (1999),0.12~1 for the fast PQ pool,  or 3~8 for the slow recycling PQ pool			
		FI_RC	(	20	)	=	k_r3		;	%	The rate constant of the exchange of QB and PQH2	Lazar (1999), since the equilibrium constant is 1 (205 in Lazar, 1999) 			
		FI_RC	(	21	)	=	k_pq_oxy		;	%	The rate constant of the PQH2 oxidation	Lazar (1999),50~500			
															
							
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialization of the initial concentration of the different component  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize the leaves for a dark adapted leaves;					
% Unit		micro mol per m2			
% Initialize the leaves for a dark adapted leaves;					
%		mircomol per m2			

% Initialize the leaves for a dark adapted leaves;					
%	 Micro mol m2			

A	=	0	;	% 	The concentration of excitons in the peripheral antenna
U	=	0	;	%	The concentration fo excitons in the core antenna
P680Pheo	=	1	;	% 	The concentration of the P680Pheo
P680pPheon	=	0	;	%	The concentration for the P680+ Pheo-
P680pPheo	=	0	;	% 	The concentration of P680+ Pheo
P680Pheon	=	0	;	%	The concentration of P680Pheo-
Yz	=	1;	    %	The concentration of reduced tyrosine
S1T	=	0.8	;	% 	The concentration of S1 in combination with reduced tyrosine 
S2T	=	0	;	%	The concentration of S2 in combination with reduced tyrosine
S3T	=	0	;	% 	The concentration of S3 in combination with reduced tyrosine
S0T	=	0.2	;	%	The concentration of S0 in combination with reduced tyrosine 
S1Tp	=	0	;	% 	The concentration of S1 in combination with oxidized tyrosine
S2Tp	=	0	;	% 	The concentration of S2 in combination with oxidized tyrosine
S3Tp	=	0	;	%	The concentration of S3 in combination with oxidized tyrosine
S0Tp	=	0	;	% 	The concentration of S0 in combination with oxidized tyrosine
QAQB	=	1	;	% 	The concentration of [QAQB]
QAnQB	=	0	;	% 	The concentration of [QA-QB];
QAQBn	=	0	;	%	The concentration of [QAQB-]
QAnQBn	=	0	;	% 	The concentration of [QA-QB-];
QAQB2n	=	0	;	%	The concentration of [QAQB2-]
QAnQB2n	=	0	;	% 	The concentration of [QA-QB2-];
PQn	=	5       ;	%	The concentration of reduced PQ, i.e. PQH2;

% Assign the value to a array
% FI_ini.m						
% This is the program that initialize the major variables used in the fluorescence induction system.In this file, the n represent negative charges, _red represent that the components are associated with the closed reaction center; while _ox represent a system with open reaction center. 						
global FI_Con;
FI_Con	(	1	)	=	A	;	% 	The concentration of excitons in the peripheral antenna
FI_Con	(	2	)	=	U	;	%	The concentration fo excitons in the core antenna
FI_Con	(	3	)	=	P680Pheo	;	% 	The concentration of the P680Pheo
FI_Con	(	4	)	=	P680pPheon	;	%	The concentration for the P680+ Pheo-
FI_Con	(	5	)	=	P680pPheo	;	% 	The concentration of P680+ Pheo
FI_Con	(	6	)	=	P680Pheon	;	%	The concentration of P680Pheo-
FI_Con	(	7	)	=	Yz	;           %	The concentration of reduced tyrosine
FI_Con	(	8	)	=	S1T	;           % 	The concentration of S1 in combination with reduced tyrosine 
FI_Con	(	9	)	=	S2T	;           %	The concentration of S2 in combination with reduced tyrosine
FI_Con	(	10	)	=	S3T	;           % 	The concentration of S3 in combination with reduced tyrosine
FI_Con	(	11	)	=	S0T	;           %	The concentration of S0 in combination with reduced tyrosine 
FI_Con	(	12	)	=	S1Tp	;	% 	The concentration of S1 in combination with oxidized tyrosine
FI_Con	(	13	)	=	S2Tp	;	% 	The concentration of S2 in combination with oxidized tyrosine
FI_Con	(	14	)	=	S3Tp	;	%	The concentration of S3 in combination with oxidized tyrosine
FI_Con	(	15	)	=	S0Tp	;	% 	The concentration of S0 in combination with oxidized tyrosine
FI_Con	(	16	)	=	QAQB	;	% 	The concentration of [QAQB]
FI_Con	(	17	)	=	QAnQB	;	% 	The concentration of [QA-QB];
FI_Con	(	18	)	=	QAQBn	;	%	The concentration of [QAQB-]
FI_Con	(	19	)	=	QAnQBn	;	% 	The concentration of [QA-QB-];
FI_Con	(	20	)	=	QAQB2n	;	%	The concentration of [QAQB2-]
FI_Con	(	21	)	=	QAnQB2n	;	% 	The concentration of [QA-QB2-];
FI_Con	(	22	)	=	PQn	;	%	The concentration of reduced PQ, i.e. PQH2;

global FI_Pool;
QBt =   1;          % The total concentration of Qb site;
PQT =   8;          % The total concentration of PQ;

FI_Pool(1) = QBt;
FI_Pool(2)  =   PQT;


global FIBF_AUX;
FIBF_AUX = zeros(5,1);

global FI_RC_Reg_o; 
FI_RC_Reg_o(1) = FI_RC(11);
FI_RC_Reg_o(2) = FI_RC(12);
FI_RC_Reg_o(3) = FI_RC(13);
FI_RC_Reg_o(4) = FI_RC(7);

