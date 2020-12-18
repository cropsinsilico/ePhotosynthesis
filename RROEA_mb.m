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



function RROEA_mb = RROEA_Mb(t,RROEA_Con,RROEA_Param)

global GLight;
fini = Condition (t);
light = GLight;

RROEA_Param(1) = light; 

RROEA_Vel = RROEA_Rate(t,RROEA_Con, RROEA_Param);

ve2GAPDH = RROEA_Vel	(	1	)			;	
ve2FBPase = RROEA_Vel	(	2	)			;
ve2SBPase = RROEA_Vel	(	3	)			;	
ve2PRK = RROEA_Vel	(	4	)			;
ve2ATPase = RROEA_Vel	(	5	)			;
ve2ATPGPP = RROEA_Vel	(	6	)			;
ve2MDH = RROEA_Vel	(	7	)		;
ve2Fd = RROEA_Vel	(	8	)		;
veFd2Thio = RROEA_Vel	(	9	)		;
veFd2Calvin = RROEA_Vel	(	10	)		;
ve2RuACT = RROEA_Vel	(	11	)		;

RROEA_mb				=	zeros(4,1)	;				
RROEA_mb	(	1	)	=	ve2GAPDH	;	%	GAPDH		
RROEA_mb	(	2	)	=	ve2FBPase;	%	FBPase		
RROEA_mb	(	3	)	=	ve2SBPase;	%	SBPase	
RROEA_mb	(	4	)	=	ve2PRK;	    %	PRK
RROEA_mb	(	5	)	=	ve2ATPase;	%	ATPase
RROEA_mb	(	6	)	=	ve2ATPGPP;	%	ATPGPP
RROEA_mb	(	7	)	=	ve2MDH;	    %	MDH
RROEA_mb	(	8	)	=	veFd2Thio - ve2GAPDH - ve2FBPase - ve2SBPase - ve2PRK - ve2ATPGPP -ve2RuACT;	    %	Thio
RROEA_mb	(	9	)	=	ve2Fd - veFd2Thio - veFd2Calvin;	    %	Fd
RROEA_mb    (   10  )   =   ve2RuACT;   % RuACT;