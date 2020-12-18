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


function StomCond_Vel = StomCond_Rate(t,StomCond_Con, StomCond_Param)

Light = StomCond_Param(1);

global StomCond_RC;

% The rate constant used in the model		

krs	=	StomCond_RC	(	1	);	        %  The rate constant of osmotic potential change following the chnages in the stromal concentration of ATP
k_turg2gs	=	StomCond_RC	(	2	);	%  The rate constant of osmotic potential change following the chnages in the stromal concentration of ATP
Liquid2Gas = StomCond_RC	(	3	);	%  The conversion of liquid to gas cocnentration
GLConv  =   StomCond_RC	(	4	);		%  The conversion from gas to liquid when expressing the flux of CO2 in. 
Vcmax = StomCond_RC	(	5	)	;       %  The maximum rate of Rubisco-limited photosynthesis
Jmax = StomCond_RC	(	6	)	;       %  The maximum rate of RuBP - regeneration limited CO2 fixation rate


Posm	=	StomCond_Con	(	1	);			
CO2	=	StomCond_Con	(	2	);	
							
global StomCond_Pool  ;		

Posm_MAX    = StomCond_Pool(1) ;	
WaterPg     = StomCond_Pool(2);	
MechAdv     = StomCond_Pool (3);	
CO2A        = StomCond_Pool (4);	
ATPm        = StomCond_Pool (5);	
Pe          = StomCond_Pool (6);	

global DropHum;

if DropHum ==0
    if t>450000         % Here is a place to add in certin regulations. 
        
        StomCond_Pool (6)  =       2;
        StomCond_Pool(2)   =       0.6;
        
        Pe          = StomCond_Pool (6);	
        WaterPg     = StomCond_Pool(2);	    
    end
    DropHum = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               Regulations                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global StomCond_TrDynaPS_com;

% global TrDynaPs_StomCond_ATP;
% global TrDynaPS_StomCond_CO2_consum;

global Temp_cond;
Temp = Temp_cond;
global CO2_cond; 
CO2A = CO2_cond * 3 * 10000; 

CO2I = CO2 * Liquid2Gas;
global AVR;         % The conversion factor between micro mole per meter square per second and milimole per liter per second. 

if StomCond_TrDynaPS_com ==1
    global BF2Stom_ATP; 
    global PS2Stom_CO2_consum;
    
    ATP = BF2Stom_ATP;
    CO2_consum = PS2Stom_CO2_consum;
    CO2A = CO2_cond * 3 * 10^4;
else
    [ATP, CO2_consum] = ssPSFun(Vcmax, Jmax, Temp, CO2I, Light); 
    CO2_consum = CO2_consum / AVR;    % Convert the unit from micromol m-2 s-1 to mmol l-1
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get the auxiliary variables %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The turgor pressure
ATP = ATP + 0; 

TurgorPressure = Posm + WaterPg - Pe * MechAdv;
Gs = TurgorPressure * k_turg2gs;        % Unit: mol m-2 s-1.    
Influx_CO2 = Gs * (CO2A - CO2I);        % Unit: mol m-2 s-1 * micromole mol-1  
Influx_CO2 = Influx_CO2 /AVR;            % Unit conversion

% Calcualte the rate equation 

% ATP  = 0.8; 

Poten_Osm = ATP/ATPm * Posm_MAX; 

% Notice here that the ABA concentration will modify this potential osmotic 
% potential that can be obtained. The higher the ABA concentration, the lower the
% maximum potential osmotic potential it can support. Therefore, the rate equation is modified as 

if ATP < 1.5
   Poten_Osm = ATP/ATPm * Posm_MAX; 
end

if ATP > 1.5
   Poten_Osm = Posm_MAX; 
end

global Xan2Stom_ABA; 

ABA = Xan2Stom_ABA;

%%%%%%%%%%%%%%%%%%
%%% Hypothesis testing  %%%????????????????????????/
%%%%%%%%%%%%%%%%%%

ABA = 0; 

%%%%%%%%%%%%%%%%%%
%%%     Hypothesis testing %%%      ????????????????????????/
%%%%%%%%%%%%%%%%%%

if StomCond_TrDynaPS_com == 1
    
    % Poten_Osm = Poten_Osm /(ABA * 100 + 1); 
    
    PercentBinding = ABA * 0.25;         % Notice here that the ABA concentration is in a ratio format. The default is 0.25. Now, let's use 0.5. 
        
    Poten_Osm = Poten_Osm * (1-PercentBinding);
    % Notice that the concentration of ABA might be way below one, therefore, one 
    % pseudonumber 1 is given here. 
end

Osm_Change = krs * (Poten_Osm - Posm);

% Calculate the rate of CO2 fixation based on steady state model

global StomCond_OLD_TIME;
global StomCond_TIME_N;
global StomCond_VEL;
global StomCond_CON;

if (StomCond_TIME_N ==0)
    StomCond_TIME_N = 1;
end

if (t > StomCond_OLD_TIME)
    StomCond_TIME_N = StomCond_TIME_N + 1;
    StomCond_OLD_TIME = t;
end

StomCond_VEL	(	StomCond_TIME_N	,	1	)	=	t;
StomCond_VEL	(	StomCond_TIME_N	,   2	)	=	Osm_Change	;	%	Osm_Change		THe rate of osmotic change
StomCond_VEL	(	StomCond_TIME_N	,   3	)	=	Influx_CO2	;	%	Influx_CO2  The rate of CO2 influx	
StomCond_VEL	(	StomCond_TIME_N	,   4	)	=   CO2_consum  ;   %   The CO2 consumption flux
StomCond_VEL	(	StomCond_TIME_N	,   5	)	=   ATP  ;   %   The CO2 consumption flux


StomCond_CON(StomCond_TIME_N,1) = t;
StomCond_CON(StomCond_TIME_N,2) = TurgorPressure;
StomCond_CON(StomCond_TIME_N,3) = Gs;
StomCond_CON(StomCond_TIME_N,4) = ATP;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           Assign table            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

StomCond_Vel	(	1	)	=	Osm_Change	;	    %	%	Osm_Change		THe rate of osmotic change							
StomCond_Vel	(	2	)	=	Influx_CO2	;	    %	Influx_CO2  The rate of CO2 influx									
StomCond_Vel	(	3	)	=	CO2_consum	;	    %	Influx_CO2  The rate of CO2 influx									

% hold on;figure(1); 
% plot(t,ATP);
%hold on; 

global StomCon2OUT; 
StomCond2OUT = zeros(5,1); 

StomCon2OUT(1) = TurgorPressure; 
StomCon2OUT(2) = Gs; 
StomCon2OUT(3) = Posm; 

% GenOut(t); 


