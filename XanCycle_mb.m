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


function XanCycle_mb = XanCycle_Mb(t,XanCycle_Con,XanCycle_Param)


fini = Condition (t);

XanCycle_Vel = XanCycle_Rate(t,XanCycle_Con, XanCycle_Param);

Vva =  XanCycle_Vel	(	1	)		;	%	The velocity of v to a conversion
Vaz =  XanCycle_Vel	(	2	)		;	%	The rate of A to z
Vza =  XanCycle_Vel	(	3	)		;	%	THe rate of z to a
Vav = XanCycle_Vel	(	4	)		;	%	The rate of A to V
Vvf = XanCycle_Vel	(	5	)		;	%	The rate of V formation
Vv2ABA =  XanCycle_Vel	(	6	)	;	%	The rate of conversion from v to ABA.
VABAdg = XanCycle_Vel	(	7	)	;	%	The rate of ABA degradation										

XanCycle_mb				=	zeros(4,1)	;				

XanCycle_mb	(	1	)	=	Vvf + Vav - Vva - Vv2ABA	;	 	
XanCycle_mb	(	2	)	=	Vva - Vav + Vza - Vaz;	 	
XanCycle_mb	(	3	)	=	Vaz - Vza	;	 
XanCycle_mb	(	4	)	=	Vv2ABA - VABAdg;	 

