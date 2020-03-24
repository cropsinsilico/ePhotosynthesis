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

% This is a function to generate output from the program. 

function OutSuc = GenOut (t)

global PS2OUT;
global PR2OUT;
global BF2OUT;
global SUCS2OUT;

RuBP = PS2OUT(1)   ;
PGA = PS2OUT(2)	   ;
DPGA = PS2OUT(3)	   ;
T3P =   PS2OUT(4)	   ;   
FBP =   PS2OUT(6)	   ;
E4P =   PS2OUT(7)	   ;
S7P =   PS2OUT(8)	   ;
SBP =   PS2OUT(9)	   ;
ATP =   PS2OUT(10)	   ;
NADPH   =   PS2OUT(11)	   ;
CO2     =   PS2OUT(12)	;    
O2      =   PS2OUT(13)	;   
HexP    =   PS2OUT(14)     ;   
PenP    =   PS2OUT(15)     ;        
Pi  =   PS2OUT(16)     ;    
ADP =   PS2OUT(17)     ;    
v1  =   PS2OUT(18)     ;

Gcea    =   PR2OUT(1)       ;
Gca     =   PR2OUT(2)       ;
Pga     =   PR2OUT(3)       ;
Pgca     =   PR2OUT(4)      ;
Gcac    =   PR2OUT(5)       ;
Goac    =   PR2OUT(6)       ;
Serc    =   PR2OUT(7)       ;
Glyc    =   PR2OUT(8)       ;
Hprc    =   PR2OUT(9)       ;
Gceac   =   PR2OUT(10)      ;
Rubp    =   PR2OUT(11)      ;
v131    =   PR2OUT(12)      ;

Fdn     =   0;
PHs     =   0;
PHl     =   0;

global O2_cond; 
O2 = O2_cond;

global FIBF_PSPR_com
if FIBF_PSPR_com ==1
    Fdn     =   BF2OUT(1);
    PHs     =   BF2OUT(2) ;
    PHl     =   BF2OUT(3) ;
    NADPH   =   BF2OUT(4);
    ATP     =   BF2OUT(5);
end


global EPS_SUCS_com;

if EPS_SUCS_com == 0
    SUCS2OUT=zeros(12,1);
end

T3Pc    =   SUCS2OUT	(	1		)	   ;	
FBPc    =   SUCS2OUT	(	2		)	   ;		
HexPc   =   SUCS2OUT	(	3		)	   ;	
F26BPc  =   SUCS2OUT	(	4		)	   ;
ATPc    =   SUCS2OUT	(	5		)	   ;
ADPc    =   SUCS2OUT	(	6		)	   ;		
OPOPc   =   SUCS2OUT	(	7		)	   ;		
UDPGc   =   SUCS2OUT	(	8		)	   ;		
UTPc   =   SUCS2OUT	    (	9		)	   ;		
SUCP    =   SUCS2OUT	(	10		)	   ;
SUC     =   SUCS2OUT	(	11		)	   ;		
PGAc    =   SUCS2OUT	(	12		)	   ;		

global XanCycle_BF_com;
global XanCycle2OUT; 

if XanCycle_BF_com ==0
    XanCycle2OUT =zeros(5,1); 
end

V = XanCycle2OUT(1); 
A = XanCycle2OUT(2); 
Z = XanCycle2OUT(3); 
ABA = XanCycle2OUT(4); 

global StomCond_TrDynaPS_com;
global StomCon2OUT;

if StomCond_TrDynaPS_com ==0
    StomCon2OUT = zeros(3,1); 
end
TurgorPressure  =   StomCon2OUT(1) ; 
Gs  =   StomCon2OUT(2) ; 
Posm    =   StomCon2OUT(3)  ; 

global CO2A;
global GLight;
global FIBF_PSPR_com;  

global PS_TIME_N

CO2A(PS_TIME_N,1) = CO2 * 3 * 10^4;   
CO2A(PS_TIME_N,2) = O2 /1.26 ;          

global V123;    
global V16;

if FIBF_PSPR_com == 1;    
     CO2A(PS_TIME_N,3) = GLight;                 
else
     CO2A(PS_TIME_N,3) = V16;                
end

global CytoPi;

CO2A(PS_TIME_N,4) = ABA; 
CO2A(PS_TIME_N,5) = ATP;               
CO2A(PS_TIME_N,6) = Gs;               


global AVR; 

CO2A	(	PS_TIME_N,7	)	=	(v1-v131) * AVR	;%   
global RuACT2RA_v61;
global RuACT_EPS_com; 

if RuACT_EPS_com ==1
    CO2A	(	PS_TIME_N,7	)	=	(RuACT2RA_v61 - v131) * AVR	;%   
end


CO2A	(	PS_TIME_N,8	)	=	PGA	;%
CO2A	(	PS_TIME_N,9	)	=	DPGA	;%
CO2A	(	PS_TIME_N,10	)	=	T3P	;%
CO2A	(	PS_TIME_N,11	)	=	FBP	;%
CO2A	(	PS_TIME_N,12	)	=	E4P	;%
CO2A	(	PS_TIME_N,13	)	=	S7P	;%
CO2A	(	PS_TIME_N,14	)	=	SBP	;%
CO2A	(	PS_TIME_N,15	)	=	ATP	;%
CO2A	(	PS_TIME_N,16	)	=	NADPH	;%
CO2A	(	PS_TIME_N,17	)	=	HexP	;%
CO2A	(	PS_TIME_N,18	)	=	PenP	;%
CO2A	(	PS_TIME_N,19	)	=	Pi	;%
CO2A	(	PS_TIME_N,20	)	=	ADP	;%
CO2A	(	PS_TIME_N,21	)	=	RuBP;%
CO2A	(	PS_TIME_N,22	)	=	Gcea	;%
CO2A	(	PS_TIME_N,23	)	=	Gca	;%
CO2A	(	PS_TIME_N,24	)	=	Pga	;%
CO2A	(	PS_TIME_N,25	)	=	Pgca	;%
CO2A	(	PS_TIME_N,26	)	=	Gcac	;%
CO2A	(	PS_TIME_N,27	)	=	Goac	;%
CO2A	(	PS_TIME_N,28	)	=	Serc	;%
CO2A	(	PS_TIME_N,29	)	=	Glyc	;%
CO2A	(	PS_TIME_N,30	)	=	Hprc	;%
CO2A	(	PS_TIME_N,31	)	=	Gceac	;%
CO2A	(	PS_TIME_N,32	)	=	Rubp	;%
CO2A	(	PS_TIME_N,33	)	=	Fdn	;%      Fdn
CO2A	(	PS_TIME_N,34	)	=	PHs	;%      PhS
CO2A	(	PS_TIME_N,35	)	=	PHl	;%      PHl
CO2A	(	PS_TIME_N,36	)	=	T3Pc	;%
CO2A	(	PS_TIME_N,37	)	=	FBPc	;%
CO2A	(	PS_TIME_N,38	)	=	HexPc	;%
CO2A	(	PS_TIME_N,39	)	=	F26BPc	;%
CO2A	(	PS_TIME_N,40	)	=	ATPc	;%
CO2A	(	PS_TIME_N,41	)	=	ADPc	;%
CO2A	(	PS_TIME_N,42	)	=	OPOPc	;%
CO2A	(	PS_TIME_N,43	)	=	UDPGc	;%
CO2A	(	PS_TIME_N,44	)	=	UTPc	;%
CO2A	(	PS_TIME_N,45	)	=	SUCP	;%
CO2A	(	PS_TIME_N,46	)	=	SUC	;   %
CO2A	(	PS_TIME_N,47	)	=	PGAc	;%
CO2A	(	PS_TIME_N,48	)	=	V	;%
CO2A	(	PS_TIME_N,49	)	=	A	;%
CO2A	(	PS_TIME_N,50	)	=	Z	;   %
CO2A	(	PS_TIME_N,51	)	=	ABA	;%
CO2A	(	PS_TIME_N,52	)	=	TurgorPressure;
CO2A	(	PS_TIME_N,53	)	=	Gs	; 
CO2A	(	PS_TIME_N,54	)	=	Posm;


CO2A    (   PS_TIME_N,100	)	=	t	;
OutSuc = 1;