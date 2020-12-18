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
 

function SUCS_mb = SUCS_Mb(t,SUCS_Con,SUCS_Param)

global GLight;

fini = Condition (t);
light = GLight;

SUCS_Param(1) = light; 

 

SUCS_Vel = SUCS_Rate(t,SUCS_Con, SUCS_Param);
 

v51	=	SUCS_Vel	(	1	);%	DHAP+GAP --FBP
v52	=	SUCS_Vel	(	2	);%	FBP --F6P + Pi
v55	=	SUCS_Vel	(	3	);%	G1P+UTP --OPOP+UDPG 
v56	=	SUCS_Vel	(	4	);%	UDPG+F6P--SUCP + UDP
v57	=	SUCS_Vel	(	5	);%	SUCP--Pi + SUC
v58	=	SUCS_Vel	(	6	);%	F26BP--F6P + Pi
v59	=	SUCS_Vel	(	7	);%	F6P + ATP --ADP + F26BP
v60	=	SUCS_Vel	(	8	);%	ATP+UDP --UTP + ADP
% v61	=	SUCS_Vel	(	9	);                  %	POPO --2PO
v62	=	SUCS_Vel	(	10	);                      %	SUC SINK 
vdhap_in	=	SUCS_Vel	(	11	)			;   %	DHAP IN
vgap_in	=	SUCS_Vel	(	12	);                  %	GAP Export from chloroplast
vpga_in	=	SUCS_Vel	(	13	)			;       %	PGA export from chloroplast
vpga_use	=	SUCS_Vel	(	14	)			;   %	PGA utilisation in chloroplast
vatpf   =   SUCS_Vel	(	15	)			;       %	ATP synthesis rate

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get the mass balance equation%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The Major Variables										

SUCS_mb	    =	    zeros(4,1)	;				
SUCS_mb(1)	=	vdhap_in + vgap_in -  2* v51		;   %	T3Pc
SUCS_mb(2)	=	v51-v52		;                           %	FBPc
SUCS_mb(3)	=	v52 -v55 -v59 + v58	-v56	;           %	HexPc
SUCS_mb(4)	=	v59 - v58		;                       %	F26BPc
SUCS_mb(5)	=   0;                                      %   vatpf - v59 - v60;   %	ATPc
SUCS_mb(6)	=	0		        ;                       %	ADPc
SUCS_mb(8)	=	v55 - v56		;   %	UDPGc
SUCS_mb(9)	=	0;%  v60 - v55		;   %	UTPc
SUCS_mb(10)	=	v56 - v57		;   %	SUCP
SUCS_mb(11)	=	v57 - v62		;   %	SUC
SUCS_mb(12)	=	vpga_in - vpga_use	;   %	pgaC