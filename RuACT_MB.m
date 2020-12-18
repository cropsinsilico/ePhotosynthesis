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



function RuACT_mb = RuACT_Mb(t,RuACT_Con,RuACT_Param)

global GLight;
fini = Condition (t);
light = GLight;

RuACT_Param(1) = light; 

RuACT_Vel = RuACT_Rate(t,RuACT_Con, RuACT_Param);

v1	=	RuACT_Vel	(	1	)	;	%	v1	The rate of ER activation due to Rubisco activase
vn1	=	RuACT_Vel	(	2	)	;	%	vn1	The rate of E inactiavtion due to binding of RuBP
v7	=	RuACT_Vel	(	3	)	;	%	v7	The rate of formation of ECMR from ECM by binding of RuBP
vn7	=	RuACT_Vel	(	4	)	;	%	vn7	The rate of actiavtion of ECMR by Rubisco activase
v6_1	=	RuACT_Vel	(	5	)	;	%	v6_1	The rate of RuBP carboxylation
v6_2	=	RuACT_Vel	(	6	)	;	%	v6_2	The rate of RuBP oxygenation


RuACT_mb				=	zeros(4,1)	;				
RuACT_mb	(	1	)	=	vn1 - v1	;	%	ER		
RuACT_mb	(	2	)	=	v1 - v7 + vn7 + v6_1 + v6_2	- vn1;	%	EAF		
RuACT_mb	(	3	)	=	v7 - vn7 - v6_1 - v6_2	;	%	ECMR	
RuACT_mb	(	4	)	=	v6_1 + v6_2 + v1 - vn1 + vn7 - v7;	%	RuBP	
