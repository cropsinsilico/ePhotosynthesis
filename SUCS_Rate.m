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



function SUCS_Vel = SUCS_Rate(t,SUCS_Con, SUCS_Param)

global SUCS_RC;

light = SUCS_Param(1);

global 	KE501	;
global 	Km511	;
global 	Km512	;
global 	Km513	;
global 	Km514	;
global  KE51;
global 	Km521	;
global 	KI521	;
global 	KI522	;
global 	KI523	;
global 	KE52	;
global 	KE531	;
global 	KE541	;
global 	Km551	;
global 	Km552	;
global 	Km553	;
global 	Km554	;
global 	KE55	;
global 	Km561	;
global 	Km562	;
global 	KI561	;
global 	KI562	;
global 	KI563	;
global 	KI564	;
global 	KI565	;
global 	KE56	;
global 	Km571	;
global 	Ki572	;
global 	KE57	;
global 	Km581	;
global 	KI581	;
global 	KI582	;
global 	Km591	;
global 	Km592	;
global 	Km593	;
global 	KI591	;
global 	KI592	;
global 	KE59	;
global 	Km601	;
global 	Km602	;
global 	Km603	;
global 	Km604	;
global 	KE60	;
global 	KE61	;
global 	Km621	;


global 	V51	;
global 	V52	;
global 	V55	;
global 	V56	;
global 	V57	;
global 	V58	;
global 	V59	;
global  V60;
global  V61;
global	V62;	
global Vdhap_in  ;    %   DHAP export from chloroplast
global Vgap_in   ;    %   GAP export from chloroplast
global Vpga_in   ;    %   PGA export from chloropalst

SUCSV51	=	V51	;	%	;		DHAP+GAP --FBP
SUCSV52	=	V52	;	%	;		FBP --F6P + Pi
SUCSV55	=	V55	;	%	;		G1P+UTP --OPOP+UDPG 
SUCSV56	=	V56	;	%	;		UDPG+F6P--SUCP + UDP
SUCSV57	=	V57	;	%	;		SUCP--Pi + SUC
SUCSV58	=	V58	;	%	;		F26BP--F6P + Pi
SUCSV59	=	V59	;	%	;		F6P + ATP --ADP + F26BP
SUCSV60	=	V60	;	%	;		ATP+UDP --UTP + ADP
SUCSV61	=	V61	;	%	;		POPO --2PO
SUCSV62	=	V62	;	%	;		SUC Sink
SUCSVdhap_in  = Vdhap_in;    %   DHAP export from chloroplast
SUCSVgap_in   = Vgap_in;    %   GAP export from chloroplast
SUCSVpga_in   = Vpga_in;    %   PGA export from chloropalst


% The rate constant used in the model	


T3Pc	=	SUCS_Con	(	1		)			;	
FBPc	=	SUCS_Con	(	2		)			;	
HexPc	=	SUCS_Con	(	3		)			;	
F26BPc	=	SUCS_Con	(	4		)			;	
ATPc	=	SUCS_Con	(	5		)			;	
ADPc	=	SUCS_Con	(	6		)			;	
OPOPc	=	SUCS_Con	(	7		)			;	
UDPGc	=	SUCS_Con	(	8		)			;	
UTPc	=	SUCS_Con	(	9		)			;	
SUCP	=	SUCS_Con	(	10		)			;	
SUC	=	SUCS_Con	(	11		)			;	
PGAc	=	SUCS_Con	(	12		)			;	

global SUCS_Pool;		

ATc = SUCS_Pool	(	1	)	    ;
UTc = SUCS_Pool	(	2	)		;
PTc = SUCS_Pool	(	3	)		;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get the auxiliary variables %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% HexP
TEMP = 1 + KE541 + 1/KE531;

G6Pc = HexPc /TEMP;
F6Pc = G6Pc/KE531;
G1Pc = HexPc * KE541 /TEMP;

% T3P
GAPc = T3Pc /(1 + KE501);
DHAPc = T3Pc * KE501/(1+KE501);

% UDP
UDPc = UTc - UTPc - UDPGc;
ADPc = ATc - ATPc; 

% OP
PiTc = PTc - 2 * ( FBPc + F26BPc) - (PGAc + T3Pc + HexPc + SUCP + UTPc + ATPc);        
Pic   = ((KE61^2 + 4 * KE61 * PiTc)^0.5 - KE61)/2;
OPOPc = PiTc - Pic;

%%% Calculate the rate equations
temp51 =Km512 * Km513 * ( 1 + GAPc/Km512 + DHAPc/Km513 + FBPc/Km511 + GAPc * DHAPc/(Km512 * Km513));
v51 = SUCSV51 * (GAPc * DHAPc - FBPc/KE51)/temp51;

% Here the regulation of FBPase activity via the F26BP need to be implemented. 
Km521AP = Km521 * (1 + F26BPc/KI523);
temp52 = Km521AP * (1 + FBPc/Km521AP + Pic/KI522 + F6Pc/KI521 + Pic * F6Pc /(KI521 * KI522));
v52 = SUCSV52 * (FBPc - F6Pc  * Pic/KE52)/temp52;

temp55 = Km551 * Km552 * ( 1 + UTPc/Km552 + G1Pc/Km551 + UDPGc/Km554 + OPOPc/Km553 + UTPc * G1Pc/(Km551 * Km552) + UDPGc * OPOPc/(Km553 * Km554));
v55 = SUCSV55 * (UTPc * G1Pc - UDPGc * OPOPc/KE55)/temp55;

temp56 = (F6Pc + Km561 * (1 + FBPc/KI562))*(UDPGc + Km562 * ( 1 + UDPc/KI561)*(1+SUCP/KI563)*(1+Pic/KI564)*(1+SUC/KI565)); 
v56 = SUCSV56 * (F6Pc * UDPGc - SUCP * UDPc/KE56)/temp56 * 2 * (HexPc/(HexPc + 2));

temp57 = SUCP + Km571 * (1 + SUC/Ki572);
v57 = SUCSV57 * ( SUCP - SUC * Pic/KE57)/temp57;


temp58 = F26BPc + Km581 * (1+F6Pc/KI581) * (1+ FBPc/0.08);    
KI583 = 1.55;      
v58 = SUCSV58 * F26BPc/(temp58*(1+Pic/KI582) * (1+F6Pc/KI583)); 
         

Km593n = Km593 * (1 + PGAc/0.182);                                
Km591n = Km591 * (1 + PGAc/0.28);

SUCSV59n = SUCSV59;

Km591 = 5; 

KmF6P = 0.55;
Km593 = KmF6P ;       

temp59 = (F6Pc + Km593 ) * (ATPc + Km591 * (1 + ADPc/KI591));        % This is the orginal equation
 
v59 = SUCSV59n * (ATPc * F6Pc - ADPc * F26BPc/KE59)/temp59 ;
 
v60 = 0 ;       

v62 = SUCSV62 * SUC/(SUC + Km621);

Vmatpf  = 0;   
vatpf = 0;      

vdhap_in = SUCSVdhap_in * Pic/(Pic + 0.6) ;
vgap_in  = SUCSVgap_in * Pic/(Pic + 0.6);

global PSPR_SUCS_com;

if PSPR_SUCS_com ==0
	vpga_in = 0;
	vpga_use = 0;      
else
    global PS2SUCSV32;
    if PS2SUCSV32 ==0
       	vpga_in = 0;
        vpga_use =0;
    else
        vpga_use = PGAc * 1.05/(PGAc + 0.6);    % WY201803   
       	vpga_in = SUCSVpga_in * Pic/(Pic + 0.6);   % WY201803   
    end
end



global SUCS_OLD_TIME;
global SUCS_TIME_N;
global SUCS_VEL;
global SUCS_CON;

if (SUCS_TIME_N ==0)
    SUCS_TIME_N = 1;
end

if (t > SUCS_OLD_TIME)
    SUCS_TIME_N = SUCS_TIME_N + 1;
    SUCS_OLD_TIME = t;
end

SUCS_VEL	(	SUCS_TIME_N	,	1	)	=	t	;%	
SUCS_VEL	(	SUCS_TIME_N	,	2	)	=	v51	;%	DHAP+GAP --FBP
SUCS_VEL	(	SUCS_TIME_N	,	3	)	=	v52	;%	FBP --F6P + Pi
SUCS_VEL	(	SUCS_TIME_N	,	4	)	=	v55	;%	G1P+UTP --OPOP+UDPG 
SUCS_VEL	(	SUCS_TIME_N	,	5	)	=	v56	;%	UDPG+F6P--SUCP + UDP
SUCS_VEL	(	SUCS_TIME_N	,	6	)	=	v57	;%	SUCP--Pi + SUC
SUCS_VEL	(	SUCS_TIME_N	,	7	)	=	v58	;%	F26BP--F6P + Pi
SUCS_VEL	(	SUCS_TIME_N	,	8	)	=	v59	;%	F6P + ATP --ADP + F26BP
SUCS_VEL	(	SUCS_TIME_N	,	9	)	=	v60	;%	ATP+UDP --UTP + ADP
SUCS_VEL	(	SUCS_TIME_N	,	10	)	=	0	;%	POPO --2PO
SUCS_VEL	(	SUCS_TIME_N	,	11	)	=	v62	;%	SUC SINK 
SUCS_VEL	(	SUCS_TIME_N	,	12	)	=	vdhap_in	;%	DHAP IN
SUCS_VEL	(	SUCS_TIME_N	,	13	)	=	vgap_in	;%	GAP Export from chloroplast
SUCS_VEL	(	SUCS_TIME_N	,	14	)	=	vpga_in	;%	PGA export from chloroplast
SUCS_VEL	(	SUCS_TIME_N	,	15	)	=	vpga_use	;%	PGA utilisation in cytosol
SUCS_VEL	(	SUCS_TIME_N	,	16	)	=	vatpf	;%	ATP synthesis rate


SUCS_CON(SUCS_TIME_N,1) = t;
SUCS_CON(SUCS_TIME_N,2) = Pic;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assign table
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUCS_Vel	(	1	)	=	v51	;%	DHAP+GAP --FBP
SUCS_Vel	(	2	)	=	v52	;%	FBP --F6P + Pi
SUCS_Vel	(	3	)	=	v55	;%	G1P+UTP --OPOP+UDPG 
SUCS_Vel	(	4	)	=	v56	;%	UDPG+F6P--SUCP + UDP
SUCS_Vel	(	5	)	=	v57	;%	SUCP--Pi + SUC
SUCS_Vel	(	6	)	=	v58	;%	F26BP--F6P + Pi
SUCS_Vel	(	7	)	=	v59	;%	F6P + ATP --ADP + F26BP
SUCS_Vel	(	8	)	=	v60	;%	ATP+UDP --UTP + ADP
SUCS_Vel	(	9	)	=	0	;%	POPO --2PO
SUCS_Vel	(	10	)	=	v62	;%	SUC SINK 
SUCS_Vel	(	11	)	=	vdhap_in	;%	DHAP IN
SUCS_Vel	(	12	)	=	vgap_in	;%	GAP Export from chloroplast
SUCS_Vel	(	13	)	=	vpga_in	;%	PGA export from chloroplast
SUCS_Vel	(	14	)	=	vpga_use	;%	PGA utilisation in cytosol
SUCS_Vel	(	15	)	=	vatpf	;%	ATP synthesis rate
 

 global SUCS2PS_Pic;
 SUCS2PS_Pic = Pic;                % This is the original code. 
 
global SUCS2PS_GAPc; 
global SUCS2PS_PGAc; 
global SUCS2PS_DHAPc; 

SUCS2PS_GAPc = GAPc; 
SUCS2PS_PGAc = PGAc; 
SUCS2PS_DHAPc = DHAPc; 

global SUCS2CM_vdhap;
global SUCS2CM_vgap;
global SUCS2CM_vpga;
 
SUCS2CM_vdhap	=	vdhap_in;   %	DHAP IN
SUCS2CM_vgap	=	vgap_in	;   %	GAP Export from chloroplast
SUCS2CM_vpga	=	vpga_in	;   %	PGA export from chloroplast
 

global SUCS2OUT;
SUCS2OUT = zeros(5,1);

SUCS2OUT	(	1		)		=   T3Pc;	
SUCS2OUT	(	2		)		=   FBPc;		
SUCS2OUT	(	3		)		=   HexPc;	
SUCS2OUT	(	4		)		=   F26BPc;
SUCS2OUT	(	5		)		=   ATPc;
SUCS2OUT	(	6		)		=   ADPc;		
SUCS2OUT	(	7		)		=   OPOPc;		
SUCS2OUT	(	8		)		=   UDPGc;		
SUCS2OUT	(	9		)		=   UTPc;		
SUCS2OUT	(	10		)		=   SUCP;
SUCS2OUT	(	11		)		=   SUC;		
SUCS2OUT	(	12		)		=   PGAc;		