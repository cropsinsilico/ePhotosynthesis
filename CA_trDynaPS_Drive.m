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




global CA_FLAG;
CA_FLAG = 1;

Begin = 1;
fin = SSI(Begin);
global options1;
global tglobal;
time = tglobal;

global ATPActive;
ATPActive = 0;

global EPS_ATP_Rate;        
EPS_ATP_Rate = 0;

ModelComb = IniModelCom;       

global BF_FI_com;            
BF_FI_com = 1;

global PR_PS_com;   
PR_PS_com = 1;

global FIBF_PSPR_com;  
FIBF_PSPR_com = 1;    

global RuACT_EPS_com;     
RuACT_EPS_com = 1;         

global RedoxReg_RA_com;     
RedoxReg_RA_com = 1;         

global XanCycle_BF_com;
XanCycle_BF_com = 1;

global RROEA_EPS_com;
RROEA_EPS_com = 1;

global CO2A;
CO2A = zeros(5,1);
 

trDynaPS_Con = trDynaPS_Ini;
 

va1 = 0;
global PS12ratio;  
BF_Param = zeros(5,1);
BF_Param(1) = va1;
BF_Param(2) = PS12ratio;

FI_Param = zeros(5,1);
FI_Param(1) = va1;
FI_Param(2) = PS12ratio;

PS_PR_Param = 0;
EPS_Param = 0;

RuACT_Param = zeros(2,1);
RuACT_Param(1) = va1;
RuACT_Param(2) = PS12ratio;

XanCycle_Param = zeros(2,1);
XanCycle_Param(1) = va1;
XanCycle_Param(2) = PS12ratio;


RROEA_Param = zeros(2,1);
RROEA_Param(1) = va1;
RROEA_Param(2) = PS12ratio;

RedoxReg_Param = 0;  

[Tt,d] = ode15s(@trDynaPS_mb,[0,time],trDynaPS_Con,options1,BF_Param, FI_Param, PS_PR_Param, RuACT_Param, RedoxReg_Param, XanCycle_Param, RROEA_Param);

    M = size(d);
    a = M(1);
    
    FIBF = zeros(a,52);
    PSPR = zeros(a,23);
    
    for m = 1:52
         FIBF(:,m)= d(:,m);
    end
    
    for m = 1:23
         PSPR(:,m)=d(:,m+52);
    end
    
    global FI;
    global BF;
    
    global PS;
    global PR;
    
    global RuACT_VEL;
    global PR_VEL;
    global CarbonRate2;
    CarbonRate2 = RuACT_VEL(:,6) * 37.5;
    


