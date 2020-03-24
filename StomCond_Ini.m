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

% StomCond_Init.m     This is the routine that initialize the parameters, initial conditions for simulation of the stomatal conductance model

function StomCond_Con = StomCond_Ini(begin)


global StomCond_OLD_TIME;
global StomCond_TIME_N;
global StomCond_VEL;
global StomCond_CON;

StomCond_OLD_TIME = 0;
StomCond_TIME_N = 1;

StomCond_VEL = zeros(1,3);    % Clean memory
StomCond_CON = zeros(3,1);    % Clean memory								

krs	=	0.02	;      
k_turg2gs = 0.105 ;              % The constant to convert the turgor pressure to stomatal conductance. From Buckley et al 2003 pce. mol m-2 s-1 Mpa-1
Liquid2GasCO2 = 30000;          % The conversion of liquid to gas cocnentration. For CO2, 0.012 mmol corresponds to 360 ppm, i.e. multiply a factor of 360/0.012 = 360/(12 * 10-3) = 10^4; UNIT: ppm
Liquid2GasO2 = 0.81 * 1000;     % The conversion of liquid to gas cocnentration. For O2, 0.23 mM corresponds to 20%, i.e. multiply a factor of 0.2/0.23 = 0.81; unit: mmol mol-1
GLConv  = 1/30000;     % The conversion from gas to liquid when expressing the flux of CO2 in. The basic idea is to convert the micro mole per meter square to mmol per liter. Based on the previous calcualtions, 
                       % one micromole per meter square leaf area has is converted to a unit of a volume of stroma as 1/27 mmol l-1

Vcmax = 70;         % Vcmax;
Jmax = 150;         % Jmax

global StomCond_RC;
StomCond_RC = zeros(1,1);

% The rate constant used in the model											
StomCond_RC	(	1	)	=	krs	;	% The rate constant of osmotic potential change following the chnages in the stromal concentration of ATP
StomCond_RC	(	2	)	=	k_turg2gs	 ;		% % The constant to convert the turgor pressure to stomatal conductance
StomCond_RC	(	3	)	=	Liquid2GasCO2;		%% The conversion of liquid to gas cocnentration
StomCond_RC	(	4	)	=	GLConv	;	 	% % The conversion from gas to liquid when expressing the flux of CO2 in. 
StomCond_RC	(	5	)	=	Vcmax	;	 	% Vcmax
StomCond_RC	(	6	)	=	Jmax	;	 	% Jmax
							

Posm	=	0.75	;	 
CO2 = 0.008 ;            


global StomCond_TrDynaPS_com;
global TrDynaPs_StomCond_ATP;

global CO2_cond;
global CO2_PS2StomCond;

if StomCond_TrDynaPS_com ==1
    CO2A = CO2_cond;
    CO2 = CO2_PS2StomCond;             
end

 		

StomCond_Con = zeros(2,1);
StomCond_Con	(	1	)	=	Posm	;	 
StomCond_Con	(	2	)	=	CO2	    ;	 

 

Posm_MAX	=	3	    ;	   
WaterPg =   0.9       ;        
MechAdv =   1.98       ;        
C02A = 370 ;                    
ATPm  = 1.5;                   
Pe = 0.3;                     

global StomCond_Pool;		
StomCond_Pool = zeros(2,1);		
StomCond_Pool (1) = Posm_MAX;	
StomCond_Pool (2) = WaterPg;	
StomCond_Pool (3) = MechAdv;	
StomCond_Pool (4) = C02A;	
StomCond_Pool (5) = ATPm ;	
StomCond_Pool (6) = Pe ;	

global DropHum;
DropHum = 0;