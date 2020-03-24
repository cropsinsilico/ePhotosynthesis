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



% FIBF_MB.m
% This function calculate the mass balance equation for the complete model of the light reactions.

function FIBF_mb = FIBF_MB(t, FIBF_Con, BF_Param, FI_Param)

% First Get the variables needed for the calcualtion step

    for m = 1:29
         BF_con(m)= FIBF_Con(m);
     end
     
    for m = 1:22
         FI_Con(m)=FIBF_Con(m+29);
    end
   
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Calculate auxilary variable, PQ             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global BF_Pool;
global FI_Pool;

PQn	    =   FI_Con	(	22	);	%	The concentration of reduced PQ, i.e. PQH2;
Qi	    =   BF_con	(	7	);	%	The binding Quinone
Qn	    =   BF_con	(	10	);	%	Q-
Qr	    =   BF_con	(	11	);	%	Q2-
ISPoQH2	=   BF_con	(	4	);	%	The complex of oxidized ion sulfer protein and reduced quinone
QHsemi	=   BF_con	(	5	);	%	Semiquinone

TQ	=   BF_Pool	(	8	);	%	The total concentration of plastoquinone in thylakoid membrane. ; Unit: micromole m-2 leaf area				

global FIBF_Pool;       % Since in this file, the FIBF is combined, therefore, we can use the common variable directly. 
TQ = FIBF_Pool(1); 

QBt =   FI_Pool (   1   );      %   The total concentration of QB site.
        
global FIBF_AUX;
PQ = TQ - QBt - PQn - Qi - Qn - Qr - ISPoQH2 - QHsemi; 
PQa = QBt + Qi + Qn + Qr + ISPoQH2 + QHsemi;      

FIBF_AUX(1) = PQ;
FIBF_AUX(2) = PQa;
BF_con(8) = PQ; 

global FIBF2FI_PQ; 
FIBF2FI_PQ = PQ; 

 global FI_RC;
 FI_RC(1) = FIBF_Con(52); 
 FI_RC(5) = FIBF_Con(52); 
 
 global BF_RC;
 BF_RC(19) = FIBF_Con(52);
 
 BF_mb = BF_MB(t,BF_con,BF_Param);
 FI_mb = FI_MB(t,FI_Con,FI_Param);
 
 % Assign the value of the calcualted BF_mb and FI_mb to FIBF_MB variable
 FIBF_mb = zeros(52,1);
 
    for m = 1:29
        FIBF_mb(m) = BF_mb(m);
    end
    
    for m = 1:22
        FIBF_mb(m+29) = FI_mb(m);
    end
    
% Now specially calcualte the mass balance equation for the rate constant of the heat dissipation

kd = FIBF_Con(52);          % The initialization of the initial rate constant for heat dissipation
PHl = BF_con(28);           % Get the PH value of the lumen
Hl  = 10^(PHl);            % Get the concentration of proton in lumen
QH = 10^(5.5)/(Hl + 10^(5.5));

RC = 0.1;                   % RC is the relaxation constant, which is one term borrowed from Laisk et al., 1997;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Here is the section implementing the nonphotochemical quenching. 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global XanCycle_BF_com;
global BF2FIBFMB_PHl; 
PHl  = BF2FIBFMB_PHl; 

dmax = 5 *10^8 * QH;

if XanCycle_BF_com ==1
    
    global XanCycle2FIBF_Xstate;
    if XanCycle2FIBF_Xstate <=0.3
        dmax = dmax ;
    else
        dmax = dmax * XanCycle2FIBF_Xstate/0.3;    
    end
end    
FIBF_mb(52) = RC * (dmax - kd); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate the proton generation rate from the model of FI and use that to calculate the lumen PH %            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global AVR; 
CoeffVol = AVR;      % The conversion factor between micromole per meter square per second and milimole per liter per second


global FI_Vel;

Vs3s0 = FI_Vel(11);     % This is the rate of state transition from S3 state to So state. This conversion is accompanied by splitting
                        % one molecular water molecules to release four protons. 
Hroe = 4 * Vs3s0 / CoeffVol;  % 27 is the conversion of unit from micromole per meter squre leaf area to mmol per liter.

global BF_Vel;
Vbf3	=	BF_Vel	(	3	)           ;
Vbf8	=	BF_Vel	(	8	)			;
Vbf11	=	BF_Vel	(	11	)			;

Hvqo1	=	Vbf8/CoeffVol	;	%	The rate of release of protons into lumen through Qo site				
Hvqo2	=	Vbf3/CoeffVol	;	%	The rate of proton release into lumen through Qo site		

Hroe;
global HPR; 

BF_mb	(	26	)	=	(Hvqo1 + Hvqo2 + Hroe - HPR * Vbf11);	%	BFHl	The proton and protonated buffer species in lumen, similarly, we can only use the buff concentration, but, the proton concentration can not be used here. 
BF_mb	(	28	)	=	-(Hvqo1 + Hvqo2 + Hroe - HPR * Vbf11)/1000/0.015; 	%   PHl  The changes in PH of lumen, 0.03 is from Curz et al., 2001, Biochemistry.

FIBF_mb (26) = BF_mb(26);
FIBF_mb (28) = BF_mb(28);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %          Calculate the PH of stroma        %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

v3      =   FI_Vel	(	31	)			;	%	v3	The rate of exchange of QAQBH2 with PQ; There is two proton uptake from stroma
v_r3    =	FI_Vel	(	32	)           ;	%	v_r3	The rate of exchange of QAQB with PQH2; there is two proton release into stroma
v3_n    =	FI_Vel	(	33	)           ;	%	v3_n	The rate of exchange of QAnQBH2 with PQ; there is two protons uptake from stroma
v_r3_n  =   FI_Vel	(	34	)	        ;	%	v_r3_n	The rate of exchange of QAnQB with PQH2; there is two protons release into stroma

GPQH2_qb = v3 - v_r3 + v3_n - v_r3_n;
vqb     =  2 * GPQH2_qb;  %   The total rate of proton uptake at the QB site of PSII.

Vqi     =   BF_Vel	(	12	)			;   %   The rate of proton uptake from the stroma side
Hvqi	=	Vqi/CoeffVol                ;   %	The rate of proton uptake from stroma at Qi site of cytbc1 complex	
vbfn2   =   BF_Vel	(	29	)			;   %   The rate of proton consumption by formation of NADPH 
Hrqb	=	vqb/CoeffVol	;	            %	Convert the unit of vqb from micormole per meter square per second to mM s-1; vqb is the rate of QB2- reduction in thylakoid membrane. 				

BF_mb	(	25	)	=	(HPR * Vbf11 -  Hrqb - Hvqi - vbfn2);	            %	BFHs	The proton and protonated buffer species in stroma. The proton concentration is not used in the MB procedure. The reason is that the proton concentration is buffered and therefore did not changed linerly with the generation of the protons.
BF_mb	(	27	)	=	-(HPR * Vbf11 -  Hrqb - Hvqi - vbfn2)/1000/0.015;   %	PHs, The changes of PH in stoma, 0.03 mol /PH from Laisk et al.

FIBF_mb (25) = BF_mb(25);
FIBF_mb (27) = BF_mb(27);
Vbf1    =       BF_Vel	(	1	);          % The rate of PQH2 utilization when forming the PQH2.ISP complex
GPQH2_t = GPQH2_qb - Vbf1 + Vqi/2;          % This is the total rate of PQH2 generation
Vbf7 = BF_Vel(7);               % The rate of consumption of PQ at the Qi site. 
CPQ     =  -GPQH2_qb + Vbf3 - Vbf7;
FIBF_mb (8)     =   0	;	        %	Q	Quinone in thylakoid membrane in free form
FIBF_mb	(	12	)	=	GPQH2_t;	%	QH2	The PQH2 concentration; the coefficient 2 represent the fact that 2 protons were taken up by one Q2-.
FIBF_mb	(	51	)	=	GPQH2_t;	%	QH2	The PQH2 concentration; the coefficient 2 represent the fact that 2 protons were taken up by one Q2-.