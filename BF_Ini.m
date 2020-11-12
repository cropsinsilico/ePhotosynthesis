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

function BF_con = BF_Ini(begin)
global BFVmax;
global cATPsyn;
global CPSi;
global cNADPHsyn;
% cATPsyn=1.0447;%1.01866 WY201803
% CPSi=1.0131;% 1.0237 WY201803
% cNADPHsyn=1.094468408;%1.0388 WY201803
K1	=	10^6	    ;	%	The rate constant for formation of ISP.QH2 complex; Vmax is used here. 			
K2	=	500	;	        %	The rate constant for ISP.QH2-->QH(semi) + ISPH(red)	
K3	=	5 * 10^7	;	%	The rate constant for QH. + cytbL --> Q + cytbL- + H+;
K4	=	5 * 10^7	;	%	The rate constant for cytbL- + cytbH --> cytbL + cytbH-			
K5	=	5 * 10^7	;	%	The rate constant for CytbH- + Q --> cytbH + Q-			

K6	=	5 * 10^7	;	%	The rate constant  for CytbH- + Q- --> cytbH + Q2-			
K7	=	10^4	;	%	The rate constant for Q binding to Qi site			
K8	=	1000	;	%	The rate constant for ISPH + CytC1 --> ISPH(ox) + CytC1+			
K9	=	8.3 * 10^6	;	%	The rate constant for the electron transport from cytc1 to cytc2			
K10	=	8 * 10 ^ 8	;	%	The rate constant for the electron transport from cytc2 to P700			

Vmax11	=	6 *cATPsyn;	%	WY201803 The maximum rate of ATP synthesis; unit: mmol l-1 s-1; default 1.47; 			
Kqi	=	10^ 3	;	    %	The rate constant for the uptake of 2 protons from stroma; Unit: s-1; the calculation is done based on a simplified assumption of first order kinetics for these reaction. The half time for the uptake of two protons is assumed, following Crofts "http://pop.life.uiuc.edu/~a-crofts/ahab/qcycle1.html"--The modified Q cycle, as <150 microsecond. 150 microsecond is used as an estimate of the halftime for this electron transfer. 			

PMODTEM = 1; 
PK	=	3.6 * 10^(-8)	* PMODTEM;	%	The permeability constant for K, cm s-1			   
PMg	=	3.6 * 10^(-8)	* PMODTEM;	%	The permeability constant for Mg, cm s-1			
PCl	=	1.8 * 10^(-8)	* PMODTEM;	%	The permeability constant for Cl, cm s-1			

Kau	=	10^10	;	%	The rate constant for exciton transfer from perpheral antenna to core antenna, see FI, unit: s-1			
Kua	=	10^10	;	%	The rate constant for exciton transfer from core antenna to peripheral antenna, SEE FI, unit: s-1			
Kf	=	6.3 * 10 ^6	;	%	The rate constant for fluorescence emission, see the note in FI, unit: s-1			
Kd 	=	2 * 10^8	;	%	The rate constant for heat dissipation; see the note for FI, unit: s-1			

K15	=	10^10*CPSi	;	%	WY 201803 The rate constant for the primary charge sepration in PSI, assuming the half time 30ps, unit: s-1			
K16	=	10^5	;	%	The rate constant for the electron transfer from electron acceptor A- to Fd, unit: s-1; following reference of Setif PQ and Bottin H, 1995, Biochemistry.			

Em_ISP	    =	0.31	;	%	The midpoint potential fo ISP; unit: V			
Em_CytC1	=	0.27	;	%	The midpoint potential for cytc1; unit: V			
Em_CytbL	=	-0.09	;	%	The midpoint potential for cytbL; unit: V			
Em_CytbH	=	0.05	;	%	The midpoint potential for cytbH; UNIT: V			
Em_CytC2	=	0.35	;	%	The midpoint potential for cytc2, unit: V			

% ISPHr + cytc1 --> ISPHox + cytc1-					
DeltaEm	 = 	Em_CytC1 - Em_ISP	;	
DeltaG	=	DeltaEm *(-9.649) * 10^4	;	
RT = 8.314 * 298;       			
KE8	=	exp(-DeltaG/RT)	;		

% cytc1- + cytc2 --> cytc1 + cytc2-					
DeltaEm	 = 	Em_CytC2 - Em_CytC1	;
DeltaG	=	DeltaEm *(-9.649) * 10^4 	;		
RT = 8.314 * 298;      				
KE9	=	exp(-DeltaG/RT)	;		

MemCap	=	0.6 * 10^(-6)	;	%	The membrane capacity, microFarady/cm2.
RVA	    =	8*10^(-10)	    ;	%	THe ratio of the lumen volume to thylakoid membrane; unit: L/cm2. The data from Cruz 2001.

KBs	=	1.1*10^(-8)	;	    %	The buffer equilibrium constant
KBl	=	5.1 * 10^(-6)	;	%	The buffer equilibrium constant

KM1ATP  = 0.12;
KM1ADP  = 0.014;            % Originally 
KM1PI   = 0.3;              % Originally 


KM2NADP = 0.05;     % From Fridlyand and Scheibe 1999
KM2NADPH= 0.035;    % From Fridlyand and Scheibe 1999
V2M = 27.8*cNADPHsyn;%WY 201803     % Calcualted based on 6.4 mmol (mg chl)-1h-1; Unit: mmol/l/s;
KE2 = 495;          % From Fridlyand paper , 1999, BBA, 1413, 1, 31-42


K1=BFVmax(1);
K2=BFVmax(2);
K3=BFVmax(3);
K4=BFVmax(4);
K5=BFVmax(5);
K6=BFVmax(6);
K7=BFVmax(7);
K8=BFVmax(8);
K9=BFVmax(9);
K10=BFVmax(10);
Vmax11=BFVmax(11)*cATPsyn;;
Kau=BFVmax(12);
Kua=BFVmax(13);
Kf=BFVmax(14);
Kd=BFVmax(15);
K15=BFVmax(16)*CPSi;
K16=BFVmax(17);
V2M=BFVmax(18)*cNADPHsyn;


% The rate constant used in the model											
global BF_RC;
BF_RC= zeros(5,1);
%	Assign values to the array for rate constant								

BF_RC	(	1	)	=	K1	;	%	The rate constant for formation of ISP.QH2 complex; unit:  per second	
BF_RC	(	2	)	=	K2	;	%	The rate constant for ISP.QH2-->QH(semi) + ISPH(red) ; unit:  per second	
BF_RC	(	3	)	=	K3	;	%	The rate constant for QH. + cytbL --> Q + cytbL- + H+	Unit: s-1
BF_RC	(	4	)	=	K4	;	%	The rate constant for cytbL- + cytbH --> cytbL + cytbH-	Unit: s-1
BF_RC	(	5	)	=	K5	;	%	The rate constant for CytbH- + Q --> cytbH + Q-	Unit: s-1

BF_RC	(	6	)	=	K6	;	%	The rate constant  for CytbH- + Q- --> cytbH + Q2-	Unit: s-1
BF_RC	(	7	)	=	K7	;	%	The rate constant for Q binding to Qi site; which assumed half time as 200 us, following Croft's website	Unit: s-1
BF_RC	(	8	)	=	K8	;	%	The rate constant for ISPH + CytC1 --> ISPH(ox) + CytC1+	Unit: s-1
BF_RC	(	9	)	=	K9	;	%	The rate constant for the electron transport from cytc1 to cytc2	Unit: s-1
BF_RC	(	10	)	=	K10	;	%	The rate constant for the electron transport from cytc2 to P700	Unit: s-1

BF_RC	(	11	)	=	Vmax11	;	%	The maximum rate of ATP synthesis	Unit: mmol l-1 s-1; The unit for the reactions occurrs in stroma is mmol l-1 s-1
BF_RC	(	12	)	=	Kqi	;	%	The rate constant for uptake of two protons from the stroma to Q2-	s-1
BF_RC	(	13	)	=	PK	;	%	The permeability constant for K	Unit: cm s-1
BF_RC	(	14	)	=	PMg	;	%	The permeability constant for Mg	Unit: cm s-1
BF_RC	(	15	)	=	PCl	;	%	The permeability constant for Cl	Unit: cm s-1

BF_RC	(	16	)	=	Kau	;	%	The rate constant for exciton transfer from perpheral antenna to core antenna, see FI	Unit: s-1
BF_RC	(	17	)	=	Kua	;	%	The rate constant for exciton transfer from core antenna to peripheral antenna, SEE FI	Unit: s-1
BF_RC	(	18	)	=	Kf	;	%	The rate constant for fluorescence emission, see the note in FI	Unit: s-1
BF_RC	(	19	)	=	Kd 	;	%	The rate constant for heat dissipation; see the note for FI	Unit: s-1
BF_RC	(	20	)	=	KE8	;	%	ISPHr + cytc1 --> ISPHox + cytc1-	Unit: s-1

BF_RC	(	21	)	=	KE9	;	%	cytc1- + cytc2 --> cytc1 + cytc2-	Unit: s-1
BF_RC	(	22	)	=	K15	;	%	The rate constant for primary charge separation in PSI	Unit: s-1
BF_RC	(	23	)	=	K16	;	%	The rate constant for electron tranfer from electron acceptor of PSI to Fd	Unit: s-1
BF_RC	(	24	)	=	MemCap	;	%	The membrane capacity	
BF_RC	(	25	)	=	RVA	;	%	The ratio of lumen volume to thylakoid membrane area	

BF_RC	(	26	)	=	KBs	;	%	The buffer equilibrium constant in stroma	
BF_RC	(	27	)	=	KBl	;	%	The buffer equilibrium constant in lumen	
BF_RC	(	28	)	=	KM1ATP	;	%	The michaelis menton constant for ATP for ATP synthesis
BF_RC	(	29	)	=	KM1ADP	;	%	The michaelis menton constant for ATP for ADP synthesis
BF_RC	(	30	)	=	KM1PI	;	%	The michaelis menton constant for ATP for PI synthesis
															
BF_RC	(	31	)	=	KM2NADP	;	%	The michaelis menten constant for NADP	Unit: mmol l-1 s-1; The unit for the reactions occurrs in stroma is mmol l-1 s-1
BF_RC	(	32	)	=	KM2NADPH	;	%	The michaelis menten constant for NADPH	Unit: mmol l-1 s-1; The unit for the reactions occurrs in stroma is mmol l-1 s-1
BF_RC	(	33	)	=	V2M	;	%	The maximum rate of NADPH formation	Unit: mmol l-1 s-1; The unit for the reactions occurrs in stroma is mmol l-1 s-1
BF_RC	(	34	)	=	KE2	;	    %	Equilibrium constatn	

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialization of the initial concentration of the different component  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize the leaves for a dark adapted leaves;					
% Unit		micro mol per m2	or mmol l-2 stroma volume		
    
%	The initialization of the model with concentration of each substrate in dark-adapted leaves				
%	The total concentration of PSII is assumed to be 1 micromole per meter square		

ISPHr	=	0	;	%	The reduced ion sulfer protein (ISPH); unit: micromole per m2
cytc1	=	1	;	%	The oxidized state of cytc1; unit: micromole per meter square
ISPo	=	1	;	%	The oxidized ion sulfer protein (ISP); unit: micromole per meter square
ISPoQH2	=	0	;	%	The complex of oxidized ion sulfer protein and reduced quinone; unit: micromole per meter square
QHsemi	=	0	;	%	Semiquinone; micromole per meter square
cytbL	=	1	;	%	The oxidized cytbL; micromole per meter square
Qi	=	0	;	    %	The binding quinone on the quinone site; micromole per meter square
Q	=	1	;	    %	Quinone; micromole per meter square
cytbH	=1	;	    %	The oxidized form of cytbH; micromole per meter square
Qn	=	0	;	    %	Q- ; unit: micromole per meter square 
Qr	=	0	;	    %	The reduced quinone Q2- ; micromole per meter square
QH2	=	5	;	    %	The reduced quinone PQH2; micromole per meter square
cytc2	=	1	;	%	oxidized cytc2; micromole per meter square
P700	=	0.5	;	%	The reduced state of P700, including both P700 and excited P700; micromole per meter square
ADP	=	0.82	;	%	ADP in stroma, from the earlier photorespiration model; mmol l-1
Pi	=	0.9	;	    %	Phosphate in stroma, from the photorespiration model; mmol l-1
ATP	=	0.68	;	%	ATP in stroma, from the photorespiration model; mmol l-1

Ks	=	10	;	    %	K ions in stroma, mM, from the literature; mmol l-1; 90 might be an default;
Mgs	=	5;%2.5	;	    %	Mg ions in stroma, mM, from the literature of the ion estimate
Cls	=	1	;	    %	Cl ions in stroma, mM, from the literature of the ion estimate


Aip	=	0	;	    %	The number of photons in peripheral antenna; micromole per meter square
U	=	0	;	    %	The number of photons in core antenna; micromole per meter square
An	=	0	;	    %	The reduced electron acceptor in PSI; micromole per meter square
Fdn 	=	0.3	;	%	The reduced ferrodoxin; micromole per meter square leaf area

BFHs	=	19.0001; %50.0001	;	%	The protonated buffer species  and free proton together in stroma; mmol l-1; The value follows Laisk and Walker, 1989. But they did not give reference about the source of this number.; default 25
BFHl	=	19.0001;	%	The protonated buffer species and free proton together in lumen; mmol l-1; The value follows Laisk and Walker, 1989. But they did not give reference about the source of this number. ; default 5

PHs     =   7   ;   %   PH of stroma
PHl     =   7   ;   %   PH of lumen
NADPH   =   0.21;   %   The NADPH concentration in stroma at dark

% Assign the value to a array
% BF_ini.m						
% 	This is the initialization step for the module of the Q cycle, and ATP synthesis steps				
global BF_con;
BF_con	(	1	)	=	ISPHr	;	%	The reduced ion sulfer protein (ISPH)
BF_con	(	2	)	=	cytc1	;	%	The oxidized state of cytc1
BF_con	(	3	)	=	ISPo	;	%	The oxidized ion sulfer protein (ISP)
BF_con	(	4	)	=	ISPoQH2	;	%	The complex of oxidized ion sulfer protein and reduced quinone
BF_con	(	5	)	=	QHsemi	;	%	Semiquinone

BF_con	(	6	)	=	cytbL	;	%	The oxidized cytbL
BF_con	(	7	)	=	Qi	;	%	The binding Quinone
BF_con	(	8	)	=	Q	;	%	Quinone
BF_con	(	9	)	=	cytbH	;	%	The oxidized form of cytbH
BF_con	(	10	)	=	Qn	;	%	Q-

BF_con	(	11	)	=	Qr	;	%	Q2-
BF_con	(	12	)	=	QH2	;	%	QH2
BF_con	(	13	)	=	cytc2	;	%	oxidized cytc2
BF_con	(	14	)	=	P700	;	%	The reduced state of P700, including both P700 and excited P700
BF_con	(	15	)	=	ADP	;	%	ADP in stroma

BF_con	(	16	)	=	Pi	;	%	Phosphate in stroma
BF_con	(	17	)	=	ATP	;	%	ATP in stroma
BF_con	(	18	)	=	Ks	;	%	K ions in stroma
BF_con	(	19	)	=	Mgs	;	%	Mg ions in stroma
BF_con	(	20	)	=	Cls	;	%	Cl ions in stroma

BF_con	(	21	)	=	Aip	;	%	The number of photons in peripheral antenna
BF_con	(	22	)	=	U	;	%	The number of photons in core antenna
BF_con	(	23	)	=	An	;	%	The reduced electron acceptor in PSI
BF_con	(	24	)	=	Fdn	;	%	The reduced ferrodoxin
BF_con	(	25	)	=	BFHs	;	%	The total concentration of proton and protonated buffer species in stroma, put in unit: mmol l-1
BF_con	(	26	)	=	BFHl	;	%	The total concentration of proton and protonated buffer species in lumen, unit: mmol l-1
										
BF_con	(	27	)	=	PHs	;	%	The PH value of the stroma
BF_con	(	28	)	=	PHl	;	%	The PH value of the lumen
BF_con	(	29	)	=	NADPH	;	%	The NADPH concentration in stroma, Unit: mmol l-1;

% Assigning the pool variables
%	The sizes of different pools in the model						
Tcyt	=	1		;	%	Unit: micromole m-2	The total concentration of cytochrome. It is assumed that the concentration of cytbL, cytbH, and cytc1 is equal as Tcyt.
Tcytc2	=	1		;	%	Unit: micromole m-2	The total concentration of cytc2, as Tcyt , the unit is micromole per meter square leaf area
TK	=	20		;	%	Unit: mmol l-1	The total concentration of potassium, 180 mM, in the bigining assuming the concentration of ions in stroma and lumen are same
TMg	=	10		;	%	Unit: mmol l-1	The total concentraiton of Mg2+, 18 mM, in the beginning assuming the concentration of the ions in the stroma and lumen are same
TCl	=	2;%1		;	%	Unit: mmol l-1	The total concentraiton of Cl-, 1 mM, assuming equal concentrations of ion concentrations in stroma and lumen in the beginning. 
TFd	=	1		;	%	Unit: micromole m-2	The total concentration of Ferrodoxin, this assumed that only 1 Fd is associated with one PSI unit.
TA	=	1		;	%	Unit: micromole m-2	The total concentration of electron acceptor of PSI
TQ	=	8		;	%	Unit: micromole m-2	The total concentration of quinone at all different states
BFTs	=	38		;	%	Unit: mmol l-1	The total concentration of buffer in stroma; 
BFTl	=	38		;	%	Unit: mmol l-1	The total concentration of buffer in lumen


P700T   =   1       ;   %   %   The total number of P700 unit: micromole m-2 leaf area				
NADPHT  =   1;      %   The total concentration of NADPH in stroma; 1 is an guessed value;

%	Assign the pools to the global pool variables		

global BF_Pool;
BF_Pool = zeros(5,1);
BF_Pool	(	1	)	=	Tcyt;	%	The total amount of cytbH or cytbL; Unit: micromole m-2 leaf area				
BF_Pool	(	2	)	=	Tcytc2;	%	The total amount of cytc; Unit: micromole m-2 leaf area				
BF_Pool	(	3	)	=	TK	;	%	The total concentration of K in both stroma and lumen. Unit: mmol l-1. In this model, it was assumed that the total concentration of K, and Mg and Cl as well, is constant. 				
BF_Pool	(	4	)	=	TMg	;	%	The total concentration of Mg in both stroma and lumen. Unit: mmol l-1. In this model, it was assumed that the total concentration of Mg, and K and Cl as well, is constant. 				
BF_Pool	(	5	)	=	TCl	;	%	The total concentration of Cl in both stroma and lumen. Unit: mmol l-1. In this model, it was assumed that the total concentration of Cl in both stroma and lumen is constant. 				
BF_Pool	(	6	)	=	TFd	;	%	The total concentration of Ferrodoxin				
BF_Pool	(	7	)	=	TA	;	%	The total concentration of the primary electron acceptor of PSI; Unit: micromole m-2 leaf area				
BF_Pool	(	8	)	=	TQ	;	%	The total concentration of plastoquinone in thylakoid membrane. ; Unit: micromole m-2 leaf area				
BF_Pool	(	9	)	=	BFTs	;	%	The total concentration of buffer in stroma; unit: mmol per liter				
BF_Pool	(	10	)	=	BFTl	;	%	The total concentration of buffer in lumen; unit: mmol per liter				
BF_Pool	(	11	)	=	P700T	;	%	The total number of P700; unit: micromole m-2 leaf area				
BF_Pool	(	12	)	=	NADPHT	;	%   The total concentration of NADPH in stroma; 1 is an guessed value;

global HPR ; 
HPR = 4.66; 

global BF2RROEA_Fdn; 
global BF2RROEA_FdT; 

BF2RROEA_Fdn = Fdn;
BF2RROEA_FdT = TFd;
