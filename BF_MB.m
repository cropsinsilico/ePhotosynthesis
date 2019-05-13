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



function BF_mb = BF_Mb(t,BF_Con,BF_Param)

%global GLight;
global LI;
%fini = Condition (t);
light = LI*30;

BF_Param(1) = light;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Calculate the rates BFrst   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

BF_Vel = BF_Rate(t,BF_Con, BF_Param);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Get the rate of different reactions%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%	Assign velocities of different reactions in the model					
Vbf1	=	BF_Vel	(	1	)			;
Vbf2	=	BF_Vel	(	2	)			;
Vbf3	=	BF_Vel	(	3	)			;
Vbf4	=	BF_Vel	(	4	)			;
Vbf5	=	BF_Vel	(	5	)			;

Vbf6	=	BF_Vel	(	6	)			;
Vbf7	=	BF_Vel	(	7	)			;
Vbf8	=	BF_Vel	(	8	)			;
Vbf9	=	BF_Vel	(	9	)			;
Vbf10	=	BF_Vel	(	10	)			;

Vbf11	=	BF_Vel	(	11	)			;
Vqi	    =	BF_Vel	(	12	)			;
Vipc	=	BF_Vel	(	13	)			;
Vicp	=	BF_Vel	(	14	)			;
Vinc	=	BF_Vel	(	15	)			;

Vinp	=	BF_Vel	(	16	)			;
Vdp	    =	BF_Vel	(	17	)			;
Vdc	    =	BF_Vel	(	18	)			;
Vfp	    =	BF_Vel	(	19	)			;
Vfc	    =	BF_Vel	(	20	)			;

Vsfd	=	BF_Vel	(	21	)			;
VsATP	=	BF_Vel	(	22	)			;
VgPQH2	=	BF_Vel	(	23	)			;
Vbf12	=	BF_Vel	(	24	)			;
Vbf13	=	BF_Vel	(	25	)			;

Vbf14	=	BF_Vel	(	26	)			;
Vbf15	=	BF_Vel	(	27	)			;
Vbf16	=	BF_Vel	(	28	)			;
vbfn2	=	BF_Vel	(	29	)			;
VsNADPH =   BF_Vel  (   30  )           ;
vcet    =   BF_Vel	(	31	)	        ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get the mass balance equation %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
BF_mb = zeros(5,1);

global AVR; 
CoeffVol = AVR; % This is the conversion factor between volume and area. 


BF_mb	(	1	)	=	Vbf2 - Vbf8	;	%	ISPHr	The reduced ion sulfer protein (ISPH)
BF_mb	(	2	)	=	Vbf9 - Vbf8	;	%	cytc1	The oxidized state of cytc1
BF_mb	(	3	)	=	Vbf8 - Vbf1	;	%	ISPo	The oxidized ion sulfer protein (ISP)
BF_mb	(	4	)	=	Vbf1 - Vbf2	;	%	ISPoQH2	The complex of oxidized ion sulfer protein and reduced quinone
BF_mb	(	5	)	=	Vbf2 - Vbf3	;	%	QHsemi	Semiquinone

BF_mb	(	6	)	=	Vbf4 - Vbf3	;	   %	cytbL	The oxidized cytbL
BF_mb	(	7	)	=	Vbf7-Vbf5 - vcet;  %	Qi	The quinone bound in the Qi site of cytbf complex  ????

% DEBUG
BF_mb	(	8	)	=	Vbf3 - Vbf7 - VgPQH2	;	%	Q	Quinone in thylakoid membrane in free form
BF_mb	(	9	)	=	Vbf6 - Vbf4 +Vbf5	;	%	cytbH	The oxidized form of cytbH
BF_mb	(	10	)	=	Vbf5 - Vbf6	+ vcet;	        %	Qn	Q-
BF_mb	(	11	)	=	Vbf6 - Vqi/2	;	%	Qr	Q2-
BF_mb	(	12	)	=	Vqi/2 - Vbf1 + VgPQH2;	%	QH2	The PQH2 concentration; the coefficient 2 represent the fact that 2 protons were taken up by one Q2-.

BF_mb	(	13	)	=	Vbf10 - Vbf9;	%	cytc2	oxidized cytc2
BF_mb	(	14	)	=	Vbf10 - Vbf15;	%	P700	The reduced state of P700, including both P700 and excited P700

BF_mb	(	16	)	=	0	;	%	Pi	Phosphate in stroma

BF_mb	(	15	)	=	VsATP - Vbf11	;	%	ADP	ADP in stroma
BF_mb	(	17	)	=	Vbf11 - VsATP	;	%	ATP	ATP in stroma


BF_mb	(	18	)	=	Vbf12	;	%	Ks	K ions in stroma
BF_mb	(	19	)	=	Vbf13	;	%	Mgs	Mg ions in stroma
BF_mb	(	20	)	=	Vbf14	;	%	Cls	Cl ions in stroma

BF_mb	(	21	)	=	Vicp + Vinp - Vipc - Vdp - Vfp	;	%	Aip	The number of photons in peripheral antenna
BF_mb	(	22	)	=	Vipc + Vinc - Vicp - Vdc - Vfc	;	%	Ui	The number of photons in core antenna
BF_mb	(	23	)	=	Vbf15 - Vbf16	;	%	An:	the reduced electron acceptor in PSI
%BF_mb	(	24	)	=	Vbf16 - vbfn2 * CoeffVol	- vcet;	%	Fdn	The reduced ferrodoxin; unit: mircomol m-2; Therefore,the the rate of NADPH formation need to micromole per meter square.
BF_mb	(	24	)	=	Vbf16/2 - vbfn2/2 * CoeffVol	- vcet/2;  

			
vqb	=	VgPQH2 	* 2;	%	The rate of quinone protonation			
roe	=	VgPQH2 * 2 	;	%	The rate of proton generation from oxygen evolution complex				

Hroe	=	roe/CoeffVol	;	%	Convert the unit of rate of oxygen evolution (roe) from micormole per meter square per second to mM s-1				
Hrqb	=	vqb/CoeffVol	;	%	Convert the unit of vqb from micormole per meter square per second to mM s-1; vqb is the rate of QB2- reduction in thylakoid membrane. 				
Hvqo1	=	Vbf8/CoeffVol	;	%	The rate of release of protons into lumen through Qo site				
Hvqo2	=	Vbf3/CoeffVol	;	%	The rate of proton release into lumen through Qo site				
sink	=	0	    ;	        %	The sinks of electron in stroma 
Hvqi	=	Vqi/CoeffVol  ;	    %	The rate of proton uptake from stroma at Qi site of cytbc1 complex				

global HPR; 
BF_mb	(	25	)	=	(HPR * Vbf11 -  Hrqb - Hvqi -  vbfn2);	                    %	BFHs	The proton and protonated buffer species in stroma. The proton concentration is not used in the MB procedure. The reason is that the proton concentration is buffered and therefore did not changed linerly with the generation of the protons.
BF_mb	(	26	)	=	(Hvqo1 + Hvqo2 + Hroe - HPR * Vbf11);	                    %	BFHl	The proton and protonated buffer species in lumen, similarly, we can only use the buff concentration, but, the proton concentration can not be used here. 
BF_mb	(	27	)	=	-(HPR * Vbf11 -  Hrqb - Hvqi -  vbfn2)/1000/0.015;         	%	PHs, The changes of PH in stoma, 0.03 mol /PH from Laisk et al.
BF_mb	(	28	)	=	-(Hvqo1 + Hvqo2 + Hroe - HPR * Vbf11)/1000/0.015; 	        %   PHl  The changes in PH of lumen, 0.03 is from Curz et al., 2001, Biochemistry.
BF_mb	(	29	)	=	vbfn2 - VsNADPH  ;                         