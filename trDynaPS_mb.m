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


% trDynaPS_mb.m  This model includes the mass balance equations for the full model of photosynthesis.

function trDynaPS_DYDT = trDynaPS_mb(t, trDynaPS_Con, BF_Param, FI_Param, PS_PR_Param, SUCS_Param, RuACT_Param, RedoxReg_Param, XanCycle_Param, RROEA_Param)

global trDynaPS2RedReg_cal
trDynaPS2RedReg_cal = 0;

for m = 1:96
    DynaPS_Con(m) = trDynaPS_Con (m);
end 


for m = 1:10
    RROEA_Con(m) = trDynaPS_Con(m + 110); 
end

light = Condition (t);

FI_Param(1) = light;
BF_Param(1) = light;

RROEA_Param= 1; 
RROEA_DYDT = RROEA_mb(t, RROEA_Con, RROEA_Param);

DynaPS_DYDT = DynaPS_mb(t, DynaPS_Con, BF_Param, FI_Param, PS_PR_Param, RuACT_Param, RedoxReg_Param, XanCycle_Param, SUCS_Param);

trDynaPS_DYDT = zeros(10,1);

for index = 1:96    
    trDynaPS_DYDT(index) = DynaPS_DYDT(index);
end

for index = 1:10
    trDynaPS_DYDT(index+110) = RROEA_DYDT(index); 
end

global RROEA2trDynaPS_ve2Fd;
global BF2RROEA_Vbf16;
global RROEA2trDynaPS_veFd2Calvin;
global BF2trDynaPS_vbfn2;

global BF2TrDynaPSMB_vcet; 
global AVR;        
%global PRGlu;
%%%WY201804
Temp = RROEA_DYDT(9) - RROEA2trDynaPS_ve2Fd + BF2RROEA_Vbf16/AVR + RROEA2trDynaPS_veFd2Calvin - BF2trDynaPS_vbfn2 - BF2TrDynaPSMB_vcet/AVR;
%Temp = RROEA_DYDT(9) - RROEA2trDynaPS_ve2Fd + BF2RROEA_Vbf16/AVR + RROEA2trDynaPS_veFd2Calvin - BF2trDynaPS_vbfn2 - BF2TrDynaPSMB_vcet/AVR-PRGlu/AVR;

trDynaPS_DYDT(119) = Temp * AVR;
trDynaPS_DYDT(24) = Temp * AVR; 
GenOut (t); 
