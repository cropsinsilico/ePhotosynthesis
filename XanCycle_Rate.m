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


function XanCycle_Vel = XanCycle_Rate(t,XanCycle_Con, XanCycle_Param)

global XanCycle_kva;
global XanCycle_kaz;
global XanCycle_kza;
global XanCycle_kav;
global XanCycle_kvf;
global XanCycle_kv2ABA;
global XanCycle_kABAdg;

Vx	=	XanCycle_Con	(	1	)			;	% The concentration of Violozanthin	
Ax	=	XanCycle_Con	(	2	)			;	%	The concentration of Anthrozanthin
Zx	=	XanCycle_Con	(	3	)		;	%	The concentration of Zeaznthin	
ABA	=	XanCycle_Con	(	4	)		;	%	The concentration of ABA
						

fini = Condition (t);
global Glight; 

light = Glight; 

global XanCycle_BF_com;
if XanCycle_BF_com ==1
    global BF2XanCycle_pHl;
    pH = BF2XanCycle_pHl;
else
    pH = 6; 
end

if pH <=5.8
    RegCof = 1;
elseif pH > 5.8 & pH < 6.5
    RegCof = (6.5-pH)/0.7;
else
    RegCof = 0;
end


Vva = Vx * XanCycle_kva * RegCof;
Vaz = Ax * XanCycle_kaz * RegCof;
Vza = Zx * XanCycle_kza;
Vav = Ax * XanCycle_kav;

Vvf     = 0;                      
Vv2ABA  = 0;         
VABAdg  = 0;        

global XanCycle_OLD_TIME;
global XanCycle_TIME_N;
global XanCycle_VEL;
global XanCycle_CON;

if (XanCycle_TIME_N ==0)
    XanCycle_TIME_N = 1;
end

if (t > XanCycle_OLD_TIME)
    XanCycle_TIME_N = XanCycle_TIME_N + 1;
    XanCycle_OLD_TIME = t;
end

XanCycle_VEL	(	XanCycle_TIME_N	,	1	)	=	t;
XanCycle_VEL	(	XanCycle_TIME_N	,   2	)	=	Vva	;	
XanCycle_VEL	(	XanCycle_TIME_N	,   3	)	=	Vaz	;	
XanCycle_VEL	(	XanCycle_TIME_N	,   4	)	=	Vza	;	
XanCycle_VEL	(	XanCycle_TIME_N	,   5	)	=	Vav	;	
XanCycle_VEL	(	XanCycle_TIME_N	,   6	)	=	Vvf	;	
XanCycle_VEL	(	XanCycle_TIME_N	,   7	)	=	Vv2ABA	;	
XanCycle_VEL	(	XanCycle_TIME_N	,   8	)	=   VABAdg;

XanCycle_CON(XanCycle_TIME_N,1) = t;
XanCycle_CON(XanCycle_TIME_N,2) = Vx;

 
XanCycle_Vel	(	1	)	=	Vva	;	%	The velocity of v to a conversion
XanCycle_Vel	(	2	)	=	Vaz	;	%	The rate of A to z
XanCycle_Vel	(	3	)	=	Vza	;	%	THe rate of z to a
XanCycle_Vel	(	4	)	=	Vav	;	%	The rate of A to V
XanCycle_Vel	(	5	)	=	Vvf	;	%	The rate of V formation
XanCycle_Vel	(	6	)	=	Vv2ABA	;	%	The rate of conversion from v to ABA.
XanCycle_Vel	(	7	)	=	VABAdg	;	%	The rate of ABA degradation
 

Xstate = Zx/(Vx + Ax+ Zx);
ABA = Vx /(Vx + Ax + Zx);

global XanCycle2FIBF_Xstate;
XanCycle2FI_Xstate = Xstate;

global Xan2Stom_ABA; 

Xan2Stom_ABA = ABA; 

global XanCycle2OUT;
XanCycle2OUT = zeros(5,1);
XanCycle2OUT(1) = Vx;
XanCycle2OUT(2) = Ax;
XanCycle2OUT(3) = Zx;
XanCycle2OUT(4) = ABA;
