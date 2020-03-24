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



clear; 

time1 = clock;
suc = PreProcess;

Begin = 1;
fin = SYSInitial(Begin);

global options1;

global tglobal;
time = tglobal;

%%%%%%%%%%%%%%%%%%%%%%%%
%   Calculation  step %
%%%%%%%%%%%%%%%%%%%%%%%%

global ATPActive;
ATPActive = 0;

global EPS_ATP_Rate;        
EPS_ATP_Rate = 0;

ModelComb = IniModelCom;        

global BF_FI_com;            
BF_FI_com = 1;

global PR_PS_com;    
PR_PS_Com = 1;

global FIBF_PSPR_com;  
FIBF_PSPR_com = 1;    

global RuACT_EPS_com;      
RuACT_EPS_com = 1;        

global RedoxReg_RA_com;    
RedoxReg_RA_com = 0;        

global XanCycle_BF_com;
XanCycle_BF_com = 1;

global RROEA_EPS_com;
RROEA_EPS_com = 1;       

global EPS_SUCS_com;
EPS_SUCS_com = 1;

global PSPR_SUCS_com;    % This is a variable indicating whether the PSPR model is actually need to be combined with SUCS or not. If 1 then means combined; 0 means not. 
PSPR_SUCS_com = 1;

global StomCond_TrDynaPS_com;
StomCond_TrDynaPS_com = 1;

global CO2A;
CO2A = zeros(5,1);

% Next is to initialize the vector. 

StomPS_Con = StomPS_Ini;

suc = ParamSet;
 
global RuACT2PS_Percent; 
RuACT2PS_Percent = 0.001; 

[Tt,d] = ode15s(@StomPS_mb,[0,time],StomPS_Con); 
 
ATPActive = 0;
   
global BF_VEL;
global FI_VEL;
global PS_VEL;
global PR_VEL; 
global StomCond_VEL;
global XanCycle_VEL;
global RROEA_VEL;
global RuACT_VEL;


time2 = clock;
TotalRunTime = time2-time1
    