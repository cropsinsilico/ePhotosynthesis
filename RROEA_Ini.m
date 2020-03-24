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




function RROEA_Con = RROEA_Ini(begin)


global RROEA_OLD_TIME;
global RROEA_TIME_N;
global RROEA_VEL;
global RROEA_CON;

RROEA_OLD_TIME = 0;
RROEA_TIME_N = 1;

RROEA_VEL = zeros(1,3);    % Clean memory
RROEA_CON = zeros(3,1);    % Clean memory

Coeff = 1; 

ke2GAPDH	=	22/60 *Coeff	;	%	The rate constant of electron transfer to GAPDH. From literature. 
ke2MDH	=	20/60	*Coeff;	        %	The rate constant of electront transfer to MDH, this rate is totally ASSUMED. 
ke2FBPase   = 1.38/60*Coeff  ;      % The rate constant of electron transfer from thioredoxin to FBPase; 1.38 default. 
ke2SBPase	=	1.65/60	*Coeff;	    %	The rate constant of electron tranfer from thioredoxin to SBPase; 1.65 default. 
ke2PRK = 59.8/60*Coeff;     % The rate constant of electron transfer from thioredoxin to PRK, Phosphoribulase kinase
ke2RubACT = 6.35/60*Coeff;  % The rate constant of electron transfer from thioredoxin to Rubisco activase
ke2Fd =   18.5*Coeff;       % The rate constant for electron transfer to ferrodoxin. This value is estimated based on the 
                            % Jmax of 180 micro mole per meter square per second. 
keFd2Thio =   10   *Coeff;      % The rate constant for electron transfer from fd to thio
keFd2Calvin =  7 * Coeff;       % The rate constant for electron transfer from fd to Calvin cycle. 
                                % Of course, this is a big assumption where the electron is transfered to NADPH 
                                % then to the Calvin cycle. This rate is much higher than the electron transfer to differnet 
                                % enzymes. 
ke2ATPGPP = 6.3/60 *Coeff;      % The transfer of electron from thioredoxin to ATPGPP

global RROEA_RC;
RROEA_RC = zeros(5,1);

% The rate constant used in the model											
RROEA_RC	(	1	)	=	ke2GAPDH	;	%	The rate constant of electron transfer to GAPDH. From literature. 
RROEA_RC	(	2	)	=	ke2MDH	;		%	The rate constant of electront transfer to MDH, this rate is totally ASSUMED. 
RROEA_RC	(	3	)	=	ke2FBPase	;	%	The rate constant of electron transfer from thioredoxin to FBPase.	
RROEA_RC	(	4	)	=	ke2SBPase	;	%	The rate constant of electron tranfer from thioredoxin to SBPase
RROEA_RC	(	5	)	=	ke2PRK	;	    %	The rate constant of electron transfer from thioredoxin to PRK, Phosphoribulase kinase
RROEA_RC	(	6	)	=	ke2RubACT	;	%	The rate constant of electron transfer from thioredoxin to Rubisco activase
RROEA_RC	(	7	)	=	ke2Fd	;	    %	The rate constant of electron transfer to fe
RROEA_RC	(	8	)	=	keFd2Thio	;	%	The rate constant of electron transfer from fd to thio
RROEA_RC	(	9	)	=	keFd2Calvin	;	    %	The rate constant of electron transfer from fd to Calvin cycle
RROEA_RC	(	10	)	=	ke2ATPGPP	;	    %	The rate constant of electron transfer to ATPGPP

% RROEA_RC = RROEA_RC  * 10;

% Here is all the equilibriun constants for the different reactions in photosystem
KEe2FBPase		=	0.311167869	;
KEe2SBPase		=	0.459194309	;
KEe2PRK		    =	0.677638775	;
KEe2ATPase		=	2.177727336	;
KEe2RuACT		=	0.677638775	;
KEe2GAPDH		=	0.044461692	;
KEe2MDH		    =	0.044461692	;
KEe2ATPGPP		=	1;
KEeFd2Thio		=	24776; 

global RROEA_KE;

RROEA_KE	(	1	)	=	KEe2FBPase	;
RROEA_KE	(	2	)	=	KEe2SBPase	;
RROEA_KE	(	3	)	=	KEe2PRK	;
RROEA_KE	(	4	)	=	KEe2ATPase	;
RROEA_KE	(	5	)	=	KEe2RuACT	;
RROEA_KE	(	6	)	=	KEe2GAPDH	;
RROEA_KE	(	7	)	=	KEe2MDH	    ;
RROEA_KE	(	8	)	=	KEe2ATPGPP	;
RROEA_KE	(	9	)	=	KEeFd2Thio	;

RROEA_KE	(	1:8 	)	=   RROEA_KE	(	1:8 	)	* 1;

% The following calculate the total concentration of different enzymes. 

global	V3	;       % Vmax of GAPDH
global	V6	;       % Vmax of FBPase
global	V9	;       % Vmax of SBPase
global	V13	;       % Vmax of PRK
global	V16	;       % Vmax of ATP synthase
global	V23	;       % Vmax of ATP glucose pyrophosphorylase
MDH_Vmax = 2;       % This value is assumed and there is no literature about it. Need to be
                    % fixed.       
global RROEA_EPS_com;

if RROEA_EPS_com == 0
    
    FC = 1;       
    fc16 = 1;      
    SC = 1;         

    V3		=	5.04 * SC	    ;	%	(Harris & Koniger, 1997)	3	GAP dehydragenase	DPGA+NADPH <->GAP + OP+NADP 
    V6		=	1.155*SC	;	%	(Harris & Koniger, 1997)	6	FBPase	FBP<->F6P+OP    1.155
    V9		=	0.168*SC * FC	;	%	(Harris & Koniger, 1997)	9	SBPase	SBP<->S7P+OP    0.168 as original value; 0.4168 was its value.
    V13		=	8.0094*SC	;	%	(Harris & Koniger, 1997)	13	Ribulosebiphosphate kinase	Ru5P+ATP<->RuBP+ADP
    V16		=	3 * SC	* fc16;	%	(Aflalo & Shavit, 1983, Davenport & McLeod, 1986)	16	ATP synthase	ADP+Pi<->ATP    1.47
    V23		=	1.68 * SC	* FC ;	%	(Latzko, Steup & Schachtele, 1981)	23	ADP-glucose pyrophosphorylase and	ADPG+Gn<->G(n+1)+ADP 0.18
end
 

NA = 100;

SA_GAPDH	=	620	;
SA_MDH	=	    184	;
SA_PRK	=	    410	;
SA_FBPase	=	119	;      
SA_SBPase	=	70	;
SA_ATPGPP	=	10	;

SA_ATPase   =   NA;
SA_RuACT	=	NA	;   

mw_GAPDH	=	147000;
mw_MDH	=	    87000;
mw_PRK	=	    40000;
mw_RuAT	=	    532000;
mw_FBPase	=	195000;
mw_SBPase	=	66000;
mw_ATPGPP	=	210000;
mw_ATPase   =   500000;

ThioT = 0.081;     
FdT =   0.081;        
RuACTT = 0.0056;  
 
global BF2RROEA_FdT;

if RROEA_EPS_com == 1;
    FdT = BF2RROEA_FdT;
end

global RROEA_Pool;		
RROEA_Pool	(	1	)	=	V3	* 1000*60/	SA_GAPDH/mw_GAPDH	;
RROEA_Pool	(	2	)	=	V6	*1000*60/	SA_FBPase/mw_FBPase	; 
RROEA_Pool	(	3	)	=	V9	*1000*60/	SA_SBPase/mw_SBPase	;
RROEA_Pool	(	4	)	=	V13	*1000*60/	SA_PRK/mw_PRK	;
RROEA_Pool	(	5	)	=	V16	*1000*60/	SA_ATPase/mw_ATPase	;
RROEA_Pool	(	6	)	=	V23	*1000*60/	SA_ATPGPP/mw_ATPGPP	;
RROEA_Pool	(	7	)	=	MDH_Vmax	*1000*60/SA_MDH/mw_MDH	;
RROEA_Pool	(	8	)	=   ThioT;
RROEA_Pool	(	9	)	=   FdT;
RROEA_Pool	(	10	)	=   RuACTT;

    
Coeff = 0.3; 

GAPDH	= RROEA_Pool	(	1	) * Coeff;	% 	The concentration of active GAPDH			
FBPase	= RROEA_Pool	(	2	)* Coeff;	%	The concentration of active FBPase
SBPase	= RROEA_Pool	(	3	)* Coeff;	% 	The concentration of active SBPase			
PRK= RROEA_Pool	(	4	)* Coeff;        %   The concentration of active PRK
ATPase = RROEA_Pool	(	5	) * Coeff;   %   The concentratino of active ATP synthase
ATPGPP = RROEA_Pool	(	6	)* Coeff;    %   The concnetratin of active ATP Glucose pyrophosphorylas 
MDH  = 0 ;     %   The concentration of active MDH
Thio = 0.081* Coeff;       % The initial concentration of reduced thioredoxin
Fd = 0.081* Coeff;          % The initial concentraiton of reduced fd
RuACT = 0.0056* Coeff;      % The concentration of active Rubisco activase


global RROEA_EPS_com;
global BF2RROEA_Fdn;

if RROEA_EPS_com == 1;
    Fd = BF2RROEA_Fdn;
end
		

RROEA_Con = zeros(3,1);

RROEA_Con	(	1	)	=	GAPDH	;	%	The initial concentration of active GAPDH
RROEA_Con	(	2	)	=	FBPase	;	%	The initial concentration of active FBPase
RROEA_Con	(	3	)	=	SBPase	;	%	The initial concentration of active SBPase
RROEA_Con	(	4	)	=	PRK	;	%	The initial concentration of actove PRK
RROEA_Con	(	5	)	=	ATPase	;	%	The initial concentration of actove ATPase
RROEA_Con	(	6	)	=	ATPGPP	;	%	The initial concentration of actove ATPGPP
RROEA_Con	(	7	)	=	MDH	;	%	The initial concentration of actove MDH
RROEA_Con	(	8	)	=	Thio	;	%	The initial concentration of reduced thioredoxin
RROEA_Con	(	9	)	=	Fd	;	%	The initial concentration of reduced ferrodoxin
RROEA_Con	(	10	)	=	RuACT	;	%	The initial concentration of active Rubisco activase


% Here defines the information for transfer between models 

global RROEA2PS_GAPDH;
global RROEA2PS_FBPase;
global RROEA2PS_SBPase;
global RROEA2PS_PRK;
global RROEA2PS_ATPase;
global RROEA2PS_ATPGPP;

 RROEA2PS_GAPDH   =   GAPDH;
 RROEA2PS_FBPase  =   FBPase;
 RROEA2PS_SBPase  =   SBPase;
 RROEA2PS_PRK     =   PRK;
 RROEA2PS_ATPase  =   ATPase;
 RROEA2PS_ATPGPP  =   ATPGPP;
