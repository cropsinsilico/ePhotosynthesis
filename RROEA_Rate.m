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



function RROEA_Vel = RROEA_Rate(t,RROEA_Con, RROEA_Param)


global RROEA_RC;

light = RROEA_Param(1);

ke2GAPDH = RROEA_RC	(	1	)			;	%	The rate constant of electron transfer to GAPDH. From literature. 
ke2MDH = RROEA_RC	(	2	)			;		%	The rate constant of electront transfer to MDH, this rate is totally ASSUMED. 
ke2FBPase = RROEA_RC	(	3	)			;	%	The rate constant of electron transfer from thioredoxin to FBPase.	
ke2SBPase = RROEA_RC	(	4	)			;	%	The rate constant of electron tranfer from thioredoxin to SBPase
ke2PRK = RROEA_RC	(	5	)			;	    %	The rate constant of electron transfer from thioredoxin to PRK, Phosphoribulase kinase
ke2RuACT = RROEA_RC	(	6	)			;	%	The rate constant of electron transfer from thioredoxin to Rubisco activase
ke2Fd = RROEA_RC	(	7	)			;	    %	The rate constant of electron transfer to fe
keFd2Thio = RROEA_RC	(	8	)			;	%	The rate constant of electron transfer from fd to thio
keFd2Calvin = RROEA_RC	(	9	)			;	    %	The rate constant of electron transfer from fd to Calvin cycle
ke2ATPGPP = RROEA_RC	(	10	)			;	    %	The rate constant of electron transfer from fd to ATPGPP


global RROEA_KE;

KEe2FBPase	=	RROEA_KE	(	1	)			;
KEe2SBPase	=	RROEA_KE	(	2	)			;
KEe2PRK	=	RROEA_KE	(	3	)			;
KEe2ATPase	=	RROEA_KE	(	4	)			;
KEe2RuACT	=	RROEA_KE	(	5	)			;
KEe2GAPDH	=	RROEA_KE	(	6	)			;
KEe2MDH	=	RROEA_KE	(	7	)			;
KEe2ATPGPP	=	RROEA_KE	(	8	)			;
KEeFd2Thio	=	RROEA_KE	(	9	)			;

GAPDH =RROEA_Con	(	1	)			;	%	The  concentration of active GAPDH
FBPase = RROEA_Con	(	2	)			;	%	The  concentration of active FBPase
SBPase = RROEA_Con	(	3	)			;	%	The  concentration of active SBPase
PRK = RROEA_Con	(	4	)			;	%	The  concentration of actove PRK
ATPase = RROEA_Con	(	5	)			;	%	The  concentration of actove ATPase
ATPGPP = RROEA_Con	(	6	)			;	%	The  concentration of actove ATPGPP
MDH = RROEA_Con	(	7	)			;	%	The  concentration of actove MDH
Thio = RROEA_Con (8    );                  %   The  concentration of 
Fd = RROEA_Con	(	9	)			;	%	The  concentration of reduced ferrodoxin
RuACT = RROEA_Con (10)  ;               % The concentration of Rubisco activase				


global RROEA_Pool;		

GAPDHT = RROEA_Pool	(	1	)	;
FBPaseT = RROEA_Pool	(	2	)		;
SBPaseT = RROEA_Pool	(	3	)		;
PRKT = RROEA_Pool	(	4		);
ATPaseT = RROEA_Pool	(	5	)	;
ATPGPPT = RROEA_Pool	(	6	)	;
MDHT = RROEA_Pool	(	7	)		;
ThioT = RROEA_Pool  (8); 
FdT = RROEA_Pool  (9); 
RuACTT = RROEA_Pool (10);


GAPDHo = GAPDHT - GAPDH;
FBPaseo = FBPaseT - FBPase;
SBPaseo = SBPaseT - SBPase;
PRKo = PRKT - PRK;
ATPaseo = ATPaseT - ATPase;
ATPGPPo = ATPGPPT - ATPGPP;
MDHo = MDHT - MDH;
Thioo = ThioT - Thio;
Fdo = FdT - Fd;
RuACTo = RuACTT - RuACT;

ve2GAPDH = ke2GAPDH * ( Thio * GAPDHo - Thioo * GAPDH/KEe2GAPDH)  ;
ve2FBPase = ke2FBPase * ( Thio * FBPaseo - Thioo * FBPase/KEe2FBPase) ;
ve2SBPase = ke2SBPase * ( Thio * SBPaseo - Thioo * SBPase/KEe2SBPase) ;
ve2PRK = ke2PRK * ( Thio * PRKo - Thioo * PRK/KEe2PRK) ;

KEe2ATPase = 1;    
ke2ATPase = 1;

ve2ATPase = ke2ATPase * ( Thio * ATPaseo - Thioo * ATPase/KEe2ATPase);
ve2ATPGPP = ke2ATPGPP * ( Thio * ATPGPPo - Thioo * ATPGPP/KEe2ATPGPP) ;
ve2MDH = ke2MDH * ( Thio * MDHo - Thioo * MDH/KEe2MDH) -MDH;

if light > 500
    ve2Fd = ke2Fd * Fdo;
else
    ve2Fd = ke2Fd * light/500 * Fdo;
end  


veFd2Thio = keFd2Thio * (Fd * Thioo - Thio * Fdo/KEeFd2Thio);

veFd2Calvin = Fd * keFd2Calvin * (FBPase/FBPaseT);      

ve2RuACT = ke2RuACT * ( Thio * RuACTo - Thioo * RuACT/KEe2RuACT);

global RROEA_OLD_TIME;
global RROEA_TIME_N;
global RROEA_VEL;
global RROEA_CON;

if (RROEA_TIME_N ==0)
    RROEA_TIME_N = 1;
end

if (t > RROEA_OLD_TIME)
    RROEA_TIME_N = RROEA_TIME_N + 1;
    RROEA_OLD_TIME = t;
end

RROEA_VEL	(	RROEA_TIME_N	,	1	)	=	t;
RROEA_VEL	(	RROEA_TIME_N	,   2	)	=	ve2GAPDH	;	
RROEA_VEL	(	RROEA_TIME_N	,   3	)	=	ve2FBPase	;
RROEA_VEL	(	RROEA_TIME_N	,   4	)	=	ve2SBPase	;
RROEA_VEL	(	RROEA_TIME_N	,   5	)	=	ve2PRK	;	
RROEA_VEL	(	RROEA_TIME_N	,   6	)	=	ve2ATPase	;	
RROEA_VEL	(	RROEA_TIME_N	,   7	)	=	ve2ATPGPP	;	
RROEA_VEL	(	RROEA_TIME_N	,   8	)	=	ve2MDH	;
RROEA_VEL	(	RROEA_TIME_N	,   9	)	=	ve2Fd	;
RROEA_VEL	(	RROEA_TIME_N	,   10	)	=	veFd2Thio	;
RROEA_VEL	(	RROEA_TIME_N	,   11	)	=	veFd2Calvin	;
RROEA_VEL	(	RROEA_TIME_N	,   12	)	=	ve2RuACT	;

RROEA_CON(RROEA_TIME_N,1) = t;
RROEA_CON(RROEA_TIME_N,2) = Thioo;

 

RROEA_Vel	(	1	)	=	ve2GAPDH	;	
RROEA_Vel	(	2	)	=	ve2FBPase	;
RROEA_Vel	(	3	)	=	ve2SBPase	;	
RROEA_Vel	(	4	)	=	ve2PRK	;
RROEA_Vel	(	5	)	=	ve2ATPase	;
RROEA_Vel	(	6	)	=	ve2ATPGPP	;
RROEA_Vel	(	7	)	=	ve2MDH	;
RROEA_Vel	(	8	)	=	ve2Fd	;
RROEA_Vel	(	9	)	=	veFd2Thio	;
RROEA_Vel	(	10	)	=	veFd2Calvin	;
RROEA_Vel	(	11	)	=	ve2RuACT	;

global RROEA2PS_GAPDH;
global RROEA2PS_FBPase;
global RROEA2PS_SBPase;
global RROEA2PS_PRK;
global RROEA2PS_ATPase;
global RROEA2PS_ATPGPP;

RROEA2PS_GAPDH = GAPDH;
RROEA2PS_FBPase =FBPase;
RROEA2PS_SBPase = SBPase;
RROEA2PS_PRK = PRK;
RROEA2PS_ATPase = ATPase;
RROEA2PS_ATPGPP = ATPGPP;
 
global RROEA2RuACT_RuAC;
RROEA2RuACT_RuAC = RuACT;

global RROEA2trDynaPS_ve2Fd;
RROEA2trDynaPS_ve2Fd = ve2Fd;

global RROEA2trDynaPS_veFd2Calvin;
RROEA2trDynaPS_veFd2Calvin = veFd2Calvin;