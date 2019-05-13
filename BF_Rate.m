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

function BF_Vel = BF_Rate(t,BF_con, BF_Param)

global BF_RC;

K1	=	BF_RC	(	1	)			;	%	The rate constant for formation of ISP.QH2 complex; unit:  per second 
K2	=	BF_RC	(	2	)			;	%	The rate constant for ISP.QH2-->QH(semi) + ISPH(red) ; unit:  per second	
K3	=	BF_RC	(	3	)			;	%	The rate constant for QH. + cytbL --> Q + cytbL- + H+	Unit: s-1
K4	=	BF_RC	(	4	)			;	%	The rate constant for cytbL- + cytbH --> cytbL + cytbH-	Unit: s-1
K5	=	BF_RC	(	5	)			;	%	The rate constant for CytbH- + Q --> cytbH + Q-	Unit: s-1

K6	=	BF_RC	(	6	)			;	%	The rate constant  for CytbH- + Q- --> cytbH + Q2-	Unit: s-1
K7	=	BF_RC	(	7	)			;	%	The rate constant for Q binding to Qi site; which assumed half time as 200 us, following Croft's website	Unit: s-1
K8	=	BF_RC	(	8	)			;	%	The rate constant for ISPH + CytC1 --> ISPH(ox) + CytC1+	Unit: s-1
K9	=	BF_RC	(	9	)			;	%	The rate constant for the electron transport from cytc1 to cytc2	Unit: s-1
K10	=	BF_RC	(	10	)			;	%	The rate constant for the electron transport from cytc2 to P700	Unit: s-1

Vmax11 = BF_RC	(	11	)			;	%	The maximum rate of ATP synthesis	Unit: mmol l-1 s-1; The unit for the reactions occurrs in stroma is mmol l-1 s-1
Kqi	=	BF_RC	(	12	)			;	%	The rate constant for uptake of two protons from the stroma to Q2-	s-1
PK	=	BF_RC	(	13	)			;	%	The permeability constant for K	Unit: cm s-1
PMg	=	BF_RC	(	14	)			;	%	The permeability constant for Mg	Unit: cm s-1
PCl	=	BF_RC	(	15	)			;	%	The permeability constant for Cl	Unit: cm s-1

Kau	=	BF_RC	(	16	)			;	%	The rate constant for exciton transfer from perpheral antenna to core antenna, see FI	Unit: s-1
Kua	=	BF_RC	(	17	)			;	%	The rate constant for exciton transfer from core antenna to peripheral antenna, SEE FI	Unit: s-1
Kf	=	BF_RC	(	18	)			;	%	The rate constant for fluorescence emission, see the note in FI	Unit: s-1
Kd 	=	BF_RC	(	19	)			;	%	The rate constant for heat dissipation; see the note for FI	Unit: s-1
KE8	=	BF_RC	(	20	)			;	%	ISPHr + cytc1 --> ISPHox + cytc1-	Unit: s-1

KE9	=	BF_RC	(	21	)			;	%	cytc1- + cytc2 --> cytc1 + cytc2-	Unit: s-1
K15	=	BF_RC	(	22	)			;	%	The rate constant for primary charge separation in PSI	Unit: s-1
K16	=	BF_RC	(	23	)			;	%	The rate constant for electron tranfer from electron acceptor of PSI to Fd	Unit: s-1
MemCap	=	BF_RC(	24	)			;	%	The membrane capacity	
RVA	=	BF_RC	(	25	)			;	%	The ratio of lumen volume to thylakoid membrane area	

KBs	=	BF_RC	(	26	)			;	%	The buffer equilibrium constant in stroma	
KBl	=	BF_RC	(	27	)			;	%	The buffer equilibrium constant in lumen	
KmATP	=	BF_RC	(	28	)			;	%	The michaelis menton constant for ATP for ATP synthesis
KmADP	=	BF_RC	(	29	)			;	%	The michaelis menton constant for ATP for ADP synthesis
KmPi	=	BF_RC	(	30	)			;	%	The michaelis menton constant for ATP for PI synthesis

KM2NADP	=	BF_RC	(	31	)			;	%	The michaelis menten constant for NADP	Unit: mmol l-1 s-1; The unit for the reactions occurrs in stroma is mmol l-1 s-1
KM2NADPH=   BF_RC(	32	)			;	%	The michaelis menten constant for NADPH	Unit: mmol l-1 s-1; The unit for the reactions occurrs in stroma is mmol l-1 s-1
V2M	=	BF_RC	(	33	)			;	%	The maximum rate of NADPH formation	Unit: mmol l-1 s-1; The unit for the reactions occurrs in stroma is mmol l-1 s-1
KE2	=	BF_RC	(	34	)			;	%	Equilibrium constatn	

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step II get the concentration of differnet component %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 	This is the initialization step for the module of the Q cycle, and ATP synthesis steps	

ISPHr	=	BF_con	(	1	)			;	%	The reduced ion sulfer protein (ISPH)
cytc1	=	BF_con	(	2	)			;	%	The oxidized state of cytc1
ISPo	=	BF_con	(	3	)			;	%	The oxidized ion sulfer protein (ISP)
ISPoQH2	=	BF_con	(	4	)			;	%	The complex of oxidized ion sulfer protein and reduced quinone
QHsemi	=	BF_con	(	5	)			;	%	Semiquinone

cytbL	=	BF_con	(	6	)		;	%	The oxidized cytbL
Qi	=	BF_con	(	7	)			;	%	The binding Quinone
Q	=	BF_con	(	8	)			;	%	Quinone
cytbH	=	BF_con	(	9	)		;	%	The oxidized form of cytbH
Qn	=	BF_con	(	10	)			;	%	Q-

Qr	=	BF_con	(	11	)			;	%	Q2-
QH2	=	BF_con	(	12	)			;	%	QH2
cytc2	=	BF_con	(	13	)		;	%	oxidized cytc2
P700	=	BF_con	(	14	)		;	%	The reduced state of P700, including both P700 and excited P700
ADP	=	BF_con	(	15	)			;	%	ADP in stroma

Pi	=	BF_con	(	16	)			;	%	Phosphate in stroma
ATP	=	BF_con	(	17	)			;	%	ATP in stroma
Ks	=	BF_con	(	18	)			;	%	K ions in stroma
Mgs	=	BF_con	(	19	)			;	%	Mg ions in stroma
Cls	=	BF_con	(	20	)			;	%	Cl ions in stroma

Aip	=	BF_con	(	21	)			;	%	The number of photons in peripheral antenna
U	=	BF_con	(	22	)			;	%	The number of photons in core antenna
An	=	BF_con	(	23	)			;	%	The reduced electron acceptor in PSI
Fdn	=	BF_con	(	24	)			;	%	The reduced ferrodoxin
BFHs=	BF_con	(	25  )		    ;	%	The proton concentration in stroma, put in unit: mol l-1

BFHl=	BF_con	(	26	)		    ;	%	The proton concentration in lumen; put in unit: mol l-1
PHs	=	BF_con	(	27	)			;	%	The PH value of stroma
PHl	=	BF_con	(	28	)			;	%	The PH value of lumen
NADPH=	BF_con	(	29	)		    ;	%	The NADPH concentratin in stroma, Unit: mmol l-1


global BF_Pool;

Tcyt	=	BF_Pool	(	1	)			;	%	The total amount of cytbH or cytbL; Unit: micromole m-2 leaf area				
Tcytc2	=	BF_Pool	(	2	)			;	%	The total amount of cytc; Unit: micromole m-2 leaf area				
TK	=	BF_Pool	(	3	)			;	%	The total concentration of K in both stroma and lumen. Unit: mmol l-1. In this model, it was assumed that the total concentration of K, and Mg and Cl as well, is constant. 				
TMg	=	BF_Pool	(	4	)			;	%	The total concentration of Mg in both stroma and lumen. Unit: mmol l-1. In this model, it was assumed that the total concentration of Mg, and K and Cl as well, is constant. 				
TCl	=	BF_Pool	(	5	)			;	%	The total concentration of Cl in both stroma and lumen. Unit: mmol l-1. In this model, it was assumed that the total concentration of Cl in both stroma and lumen is constant. 				
TFd	=	BF_Pool	(	6	)			;	%	The total concentration of Ferrodoxin				
TA	=	BF_Pool	(	7	)			;	%	The total concentration of the primary electron acceptor of PSI; Unit: micromole m-2 leaf area				
TQ	=	BF_Pool	(	8	)			;	%	The total concentration of plastoquinone in thylakoid membrane. ; Unit: micromole m-2 leaf area				
BFTs	=	BF_Pool	(	9	)			;	%	The total concentration of buffer in stroma; unit: mmol per liter				
BFTl	=	BF_Pool	(	10	)			;	%	The total concentration of buffer in lumen; unit: mmol per liter	
P700T	=	BF_Pool	(	11	)			;	%	The total concentration of P700; unit: micromole m-2 leaf area				
NADPHT  =   BF_Pool (12)                ;   %   The total concentration of NADPH, Unit: mmol l-1.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get the auxiliary variables %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
cytc1n	=	Tcyt - cytc1	;	%	The concentration of reduced cytc1
cytbLn	=	Tcyt - cytbL	;	%	The concentration of reduced cytbL
cytbHn	=	Tcyt - cytbH	;	%	The concentration of reduced cytbH
cytc2n	=	Tcytc2 - cytc2	;	%	The concentration of reduced cytc2
Kl	=	TK - Ks	;	%	The concentration of K in lumen
Mgl	=	TMg - Mgs	;	%	The concentration of Mg in lumen 
Cll	=	TCl - Cls	;	%	The concentration of Cl in lumen 
Fd	=	TFd - Fdn	;	%	The conncentration of oxidized Fd in stroma
A	=	TA - An	;	%	The concentration of oxidized electron acceptor in PSI
QST     =   Tcyt;   % Assuming that the total number of cytochrome is equal to the total number of quinone binding site;
% QSe     =   QST - Qi ;
QSe     =   QST - Qi - Qn - Qr;
P700p   =   P700T - P700;   % The number of positive P700;  
NADP    =   NADPHT - NADPH;

global FIBF_PSPR_com;  
global PS2BF_ADP;
global PS2BF_Pi;

if FIBF_PSPR_com ==1
    ADP = PS2BF_ADP;
    Pi  = PS2BF_Pi;
end

global AVR; 
CoeffVol = AVR;     

Hfs = 10^(-PHs)*1000;               % Calcualation of the proton concentration in stroma based on the calcualted PH value; unit: mmol l-1;
Hfl = 10^(-PHl)*1000;               % Calcualation of the proton concentration in lumen based on the calcualted PH value; unit: mmol l-1;
OHs	= 10^(-14)/(Hfs/1000)*1000	;	% Calculation of the concentration of hydroxyl ions; unit: mmol l-1;

BFs =  BFHs - Hfs;                  % The concentration of the protonated buffer speices. Notice here the variable BFHs represent the total concentration of both proton and the protonated buffer species in stroma
BFns=  BFTs - BFs;                  % The total concentration of deprotonated buffer species in stroma; 

RegPHATP = 1; 
RegPHs = 1;

Temp1 = 10^(5.5-8);   
a = Temp1/(1 + Temp1);
Temp2 = 10^(5.5-PHl);   
b = Temp2/(1+Temp2);
RegPHl = 1 - (b-a);
Iin =   BF_Param(1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate the rate of different reactions %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Vmax    =   K1 * (ISPo + ISPHr);        %   The maximum rate of formation of enzyme substrate complex
Vbf1	=	Vmax * ISPo/(ISPo+ISPHr) * QH2/TQ;	;	%	Unit: micromole s-1 m-2 leaf area
Vbf2	=	K2 * ISPoQH2* RegPHl;	%	Unit: micromole s-1 m-2 leaf area   
Vbf3	=	K3 * QHsemi * cytbL /Tcyt *RegPHl;	%	Unit: micromole s-1 m-2 leaf area
Vbf4	=	K4 * cytbLn * cytbH/Tcyt;		%	Unit: micromole s-1 m-2 leaf area
% Vbf4	=	K4 * cytbLn * cytbH/TcytK3 * QHsemi * cytbL /Tc	;	%	Unit: micromole s-1 m-2 leaf area
Vbf5	=	K5 * cytbHn * Qi/Tcyt;		%	;Unit: micromole s-1 m-2 leaf area
Vbf6	=	K6 * cytbHn * Qn/Tcyt	;	%	Unit: micromole s-1 m-2 leaf area
Vbf7	=	K7 * Tcyt * QSe/QST * Q/TQ	;	%	QSe ans QST represent the empty quinone binding site and the total number of quinone binding site respectively. Unit: micromole s-1 per meter squareleaf area;
Vbf8	=	K8 * Tcyt * (ISPHr/Tcyt * cytc1/Tcyt - ISPo/Tcyt * cytc1n/Tcyt/KE8) * RegPHl	;	    %	Unit: micromole s-1 m-2 leaf area
Vbf9	=	K9 * Tcyt * (cytc1n/Tcyt * cytc2/Tcyt - cytc1/Tcyt * cytc2n/Tcyt/KE9)	;	%	Unit: micromole s-1 m-2 leaf area
KE10 = 10;
Vbf10	=	K10 * P700p * cytc2n/P700T	- K10 * cytc2 * P700/P700T/KE10;	%	Unit: micromole s-1 m-2 leaf area

if Hfs < 0;
    Ytemp = 0;
else   
    Ytemp	=	25.359- 0.0425 * PHs ^3 + 1.0986 * PHs ^ 2 - 9.1831* PHs;		
end

Vqi	=	Kqi * Qr * 2 / 0.72 * Ytemp 	;	%	The rate of proton uptake at Qi site of bc1 complex; Unit: micromole s-1 m-2 leaf area; the coefficient 2 represents that two protons were taken up at the same time.0.72 represents the calculated proton uptake rate at Qi site when PH is 9.

% 		Calculation of the different excition transfer reactions occurred in PSI

Vipc	=	Aip * Kau	;	%	The rate of exciton transfer from peripheral antenna to core antenna; unit: micromole m-2 s-1	Unit: micromole m-2 leaf area per second
Vicp	=	U * Kua	;	%	The rate of exciton transfer from core antenna to peripheral antenna	Unit: micromole m-2 leaf area per second
 
global ChlT2;   % total Chl in PSII
% global ChlT;
global ChlPSI;  % total Chl in PSI
% Iin is the total absorbed light
 
Vinc    =   Iin * ChlPSI/(ChlT2+ChlPSI)* 95/184;    %   PPFD absorbed by core antenna of PSI    Unit: micromole m-2 leaf area per second
Vinp    =   Iin * ChlPSI/(ChlT2+ChlPSI)* 105/184;   %   PPFD absorbed by peripheral antenna of PSI  Unit: micromole m-2 leaf area per second
  
% Vinc	=	Iin * (200/(200+290))*80/200	;	%	PPFD absorbed by core antenna of PSI	Unit: micromole m-2 leaf area per second
% Vinp	=	Iin * (200/(200+290))*120/200	;	%	PPFD absorbed by peripheral antenna of PSI	Unit: micromole m-2 leaf area per second

Vdp	=	Aip * Kd	;	%	The rate of heat dissipation from peripheral antenna	Unit: micromole m-2 leaf area per second
Vdc	=	U * Kd	;	%	The rate heat dissipation from core antenna	Unit: micromole m-2 leaf area per second
Vfp	=	Aip * Kf	;	%	The fluorescence emission from peripheral antenna	Unit: micromole m-2 leaf area per second
Vfc	=	U * Kf	;	%	The fluorescence emission from core antenna	Unit: micromole m-2 leaf area per second

% 		The other empirical rates calculations				
Vsfd	=	0	;	%	The sink for Fd utilization; unit: mmol l-1 s-1;	

global CO2_cond; 
CO2 = CO2_cond * 3 * 10^4; 

MaxCO2Rate = 100 * CO2/(CO2 + 460);

VsATP	=	MaxCO2Rate * 1.5 /CoeffVol * ATP/1.5; %(ADP + ATP);	%	The sink for ATP utilizaiton, 20 represent the of CO2 assimilation, since 1 meter square amount to 27 ml, therefore, the sink capacity should be 20 * 1.5 * 1.5 mmol / 27 l-1 s-1. The 1.5 represents the 1.5 ATP consumption per CO2 fixation. 	Unit: mmol l-2 s-1
VsNADPH =  MaxCO2Rate /CoeffVol * 1 * NADPH/NADPHT; % For 6 C6 = 5 C6 + 1C6 ; 
VgPQH2  =   Q * 800 * RegPHs;    % Assuming that the rate of generation of PQH2 through QB site only depend on the PQ and PQH2 exchange capacity.
FC     	=	9.6 * 10^4	;	%	The Faraday constant;
R	    =	8.314	;	    %	J K-1 mol-1
T	    =	298	;	        %	K				
NetCharge	=	Hfs + Ks + 2 * Mgs - OHs - Cls - BFns;	%	The difference between the positive and negative charge in stroma. It was assumed that the charge is in equilibrium state in the beginning of the model, therefore, the difference in the positive and negative charges reflect the charges forming electrical potential cross the membrane. The unit is mmol l-1.
NetCharge	=	NetCharge/1000	;	                %	The unit conversion. Convert from mmol l-1 to mol l-1.
NetCharge	=	NetCharge * RVA;	                %	Convert from mol l-1 to mol cm-2.
AfC	        =	6.022*10^23	;	                    %	Avagadro's constant
UnitCharge	=	1.6*10^(-19)	;	                %	The amount of charges for one electron
NetCharge	=	NetCharge * AfC * UnitCharge;	    %	Unit: Coulomb cm-2
MPotential	=	2* NetCharge/6*10^6	;	            %	The membrane potential calculation; unit for NetCharge is : Colum per centimeter square;  Unit for MemCap is microFarady per centimeter square; the final unit is: unit is: microfarady per centimeter
temp	    =	- MPotential/0.026;	                %	0.026 pools together the RT/F; With R as 8.314 J K-1 mol-1; T: 298; F: 9.648*10^4 C mol-1. 

if temp == 0
    JK	=	PK *  (Kl - Ks )/10	;	    %	The flux of K from lumen to stroma; unit: mol/dm2/s; the unit of permeability is: cm s-1; 10 represent the conversion from cm s-1 to dm s-1; Remember that 1 liter is 1 dm3. 
    JCl	=	PCl *  (Cll - Cls )/10	;	%	The flux of Cl from lumen to stroma; unit: mol/dm2/s; unit of permeability is cm s-1; 10 represent the conversion from cm s-1 to dm s-1; remember that 1 liter amount to 1 decimeter cube
    JMg	=	PMg * (Mgl - Mgs )/10;	    %	The flux of Mg from lumen to stroma; unit: mol/dm2/s; unit of permeability is cm s-1; 10 represent the conversion of cm s-1 to dm s-1. Remember that 1 liter amount to 1 decimeter cube. 
else
    JK	=	PK * temp * (Kl - Ks * exp(-temp))/(1-exp(-temp))/10	;	    %	The flux of K from lumen to stroma; unit: mol/dm2/s; the unit of permeability is: cm s-1; 10 represent the conversion from cm s-1 to dm s-1; Remember that 1 liter is 1 dm3. 
    JCl	=	-PCl * temp * (Cll - Cls * exp(temp))/(1-exp(temp))/10	;	    %	The flux of Cl from lumen to stroma; unit: mol/dm2/s; unit of permeability is cm s-1; 10 represent the conversion from cm s-1 to dm s-1; remember that 1 liter amount to 1 decimeter cube
    JMg	=	PMg * temp/2 * (Mgl - Mgs * exp(-temp/2))/(1-exp(-temp/2))/10;	%	The flux of Mg from lumen to stroma; unit: mol/dm2/s; unit of permeability is cm s-1; 10 represent the conversion of cm s-1 to dm s-1. Remember that 1 liter amount to 1 decimeter cube. 
end

%	The ion flux needed to converted to convert to the changes in the substrate concentration changes in the stroma. 				
JKc	=	JK / RVA/100;	%	100 represent the conversion to cm-2;  RVA is the ratio between the lumen volume and thylakoid membrane area.
JClc	=	JCl / RVA/100;	%	100 represent the conversion to cm-2;  RVA is the ratio between the lumen volume and thylakoid membrane area.
JMgc	=	JMg / RVA/100;	%	100 represent the conversion to cm-2;  RVA is the ratio between the lumen volume and thylakoid membrane area.

Vbf12	=	JKc	    ;	%	Assign the variable
Vbf13	=	JMgc	;	%	Assign the variable
Vbf14	=	JClc	;	%	Assign the variable

P700e	=	U * (P700T - P700p)/120	;	%	The amount of excited P700; micromole m-2 leaf area
Vbf15	=	P700e * K15 * A/TA	;	%	The rate of PSI primary charge separation; unit: micromole m-2 leaf area per second
Vbf16	=	An * K16 * Fd/TFd	;	%	The rate of electron transport from the electron acceptor of PSI to Fd; Unit: micromole m-2 leaf area s-1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global HPR; 
RT	=	298 * 8.314	;	%	Gas constnat and the temperature
DeltaGo     =   7.3 * 4184        ;   %   The free energy change for ATP synthesis from ADP and Pi
DiffPH = PHs - PHl;
DeltaG11	=	DeltaGo -  2.3 * RT * HPR * DiffPH + HPR * MPotential * 9.6 * 10^4	;   %	The free energy change of ATP synthesis
KE11	=	exp(-DeltaG11/(RT));      %	The equilibrium constant of ATP synthesis

global DPH;

DPH = DiffPH - MPotential/0.059;

Temp	=	Vmax11 * (ADP *Pi - ATP/KE11)/((KmADP * KmPi)*(1 + ADP/KmADP + Pi/KmPi +ADP*Pi/(KmADP * KmPi)))	;	%	Unit: mmol l- s-1; The stroma volume is used as a basis for the volume
Vbf11 = Temp; 

if Vbf11 < 0
    Vbf11 = 0; 
end

global EPS_ATP_Rate;   % The EPS_ATP_Rate is used in the overall model for the calculation of the mass balance equation of ATP.
EPS_ATP_Rate = Vbf11;

FD = Fd;
FD_N = Fdn;
T_FD = TFd;

vbfn2 = 2 * V2M*(FD_N*NADP/T_FD-FD*NADPH/(T_FD*KE2))/(KM2NADP*(1+NADP/KM2NADP+NADPH/KM2NADPH));  % mmol/l/s  %QF add 2* 

vcet = V2M * Qi * FD_N/T_FD * CoeffVol;  

global BF2EPS_vbfn2;
BF2EPS_vbfn2 = vbfn2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Transfer infroamtion between models                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global FIBF_RA_Mg;
global FIBF_RA_ATP;
global BF2XanCycle_pHl;

FIBF_RA_Mg = Mgs;
FIBF_RA_ATP = ATP;
BF2XanCycle_pHl = PHl;


global BF2RROEA_Vbf16;
BF2RROEA_Vbf16 = Vbf16;

global BF2trDynaPS_vbfn2;
BF2trDynaPS_vbfn2 = vbfn2;

global BF2Stom_ATP; 
BF2Stom_ATP = ATP; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Part V Output of Velocity for plot %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global BF_OLD_TIME;
global BF_TIME_N;
global BF_VEL;
global BF_CON;

if (BF_TIME_N ==0)
    BF_TIME_N = 1;
end

if (t > BF_OLD_TIME)
    BF_TIME_N = BF_TIME_N + 1;
    BF_OLD_TIME = t;
end

BF_VEL	(	BF_TIME_N	,	1	)	=	t	    ;
BF_VEL	(	BF_TIME_N	,	2	)	=	Vbf1	;
BF_VEL	(	BF_TIME_N	,	3	)	=	Vbf2	;
BF_VEL	(	BF_TIME_N	,	4	)	=	Vbf3	;
BF_VEL	(	BF_TIME_N	,	5	)	=	Vbf4	;

BF_VEL	(	BF_TIME_N	,	6	)	=	Vbf5	;
BF_VEL	(	BF_TIME_N	,	7	)	=	Vbf6	;
BF_VEL	(	BF_TIME_N	,	8	)	=	Vbf7	;
BF_VEL	(	BF_TIME_N	,	9	)	=	Vbf8	;
BF_VEL	(	BF_TIME_N	,	10	)	=	Vbf9	;

BF_VEL	(	BF_TIME_N	,	11	)	=	Vbf10	;
BF_VEL	(	BF_TIME_N	,	12	)	=	Vbf11	;
BF_VEL	(	BF_TIME_N	,	13	)	=	Vqi	    ;
BF_VEL	(	BF_TIME_N	,	14	)	=	Vipc	;
BF_VEL	(	BF_TIME_N	,	15	)	=	Vicp	;

BF_VEL	(	BF_TIME_N	,	16	)	=	Vinc	;
BF_VEL	(	BF_TIME_N	,	17	)	=	Vinp	;
BF_VEL	(	BF_TIME_N	,	18	)	=	Vdp	    ;
BF_VEL	(	BF_TIME_N	,	19	)	=	Vdc	    ;
BF_VEL	(	BF_TIME_N	,	20	)	=	Vfp	    ;

BF_VEL	(	BF_TIME_N	,	21	)	=	Vfc	    ;
BF_VEL	(	BF_TIME_N	,	22	)	=	Vsfd	;
BF_VEL	(	BF_TIME_N	,	23	)	=	VsATP	;
BF_VEL	(	BF_TIME_N	,	24	)	=	VgPQH2	;
BF_VEL	(	BF_TIME_N	,	25	)	=	Vbf12	;

BF_VEL	(	BF_TIME_N	,	26	)	=	Vbf13	;
BF_VEL	(	BF_TIME_N	,	27	)	=	Vbf14	;
BF_VEL	(	BF_TIME_N	,	28	)	=	Vbf15	;
BF_VEL	(	BF_TIME_N	,	29	)	=	Vbf16	;
BF_VEL	(	BF_TIME_N	,	30	)	=	vbfn2	;
BF_VEL	(	BF_TIME_N	,	31	)	=	VsNADPH	;

BF_CON(BF_TIME_N,1) = t;
BF_CON(BF_TIME_N,2) = MPotential;

global BF_Vel;
BF_Vel	(	1	)	=	Vbf1	;
BF_Vel	(	2	)	=	Vbf2	;
BF_Vel	(	3	)	=	Vbf3	;
BF_Vel	(	4	)	=	Vbf4	;
BF_Vel	(	5	)	=	Vbf5	;

BF_Vel	(	6	)	=	Vbf6	;
BF_Vel	(	7	)	=	Vbf7	;
BF_Vel	(	8	)	=	Vbf8	;
BF_Vel	(	9	)	=	Vbf9	;
BF_Vel	(	10	)	=	Vbf10	;

BF_Vel	(	11	)	=	Vbf11	;
BF_Vel	(	12	)	=	Vqi	;
BF_Vel	(	13	)	=	Vipc	;
BF_Vel	(	14	)	=	Vicp	;
BF_Vel	(	15	)	=	Vinc	;

BF_Vel	(	16	)	=	Vinp;
BF_Vel	(	17	)	=	Vdp	;
BF_Vel	(	18	)	=	Vdc	;
BF_Vel	(	19	)	=	Vfp	;
BF_Vel	(	20	)	=	Vfc	;

BF_Vel	(	21	)	=	Vsfd	;
BF_Vel	(	22	)	=	VsATP	;
BF_Vel	(	23	)	=	VgPQH2	;
BF_Vel	(	24	)	=	Vbf12	;
BF_Vel	(	25	)	=	Vbf13	;

BF_Vel	(	26	)	=	Vbf14	;
BF_Vel	(	27	)	=	Vbf15	;
BF_Vel	(	28	)	=	Vbf16	;
BF_Vel	(	29	)	=	vbfn2	;
BF_Vel	(	30	)	=	VsNADPH	;
BF_Vel	(	31	)	=	vcet ;

global BF2OUT;
BF2OUT = zeros(5,1);
BF2OUT(1) = Fdn;
BF2OUT(2) = PHs;
BF2OUT(3) = PHl;
BF2OUT(4) = NADPH;
BF2OUT(5) = ATP;

global BF2FIBFMB_PHl; 

BF2FIBFMB_PHl = PHl; 

global BF2TrDynaPSMB_vcet; 
BF2TrDynaPSMB_vcet = vcet; 
