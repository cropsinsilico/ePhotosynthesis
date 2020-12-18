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

 
function StomPS_DYDT = StomPS_mb(t, StomPS_Con); 

[FI_Param, BF_Param, PS_PR_Param, SUCS_Param, EPS_Param, RuACT_Param, XanCycle_Param, RROEA_Param, RedoxReg_Param, StomCond_Param] = ParamSet;

global trDynaPS2RedReg_cal
trDynaPS2RedReg_cal = 0;

 
for m = 1:120           % TEST 
    trDynaPS_Con(m) = StomPS_Con (m);
end

for m = 1:2
    StomCond_Con(m) = StomPS_Con(m + 130);
end

 
light = Condition (t);

FI_Param(1) = light;
BF_Param(1) = light;

trDynaPS_DYDT = trDynaPS_mb(t, trDynaPS_Con, BF_Param, FI_Param, PS_PR_Param, SUCS_Param, RuACT_Param, RedoxReg_Param, XanCycle_Param, RROEA_Param);
StomCond_DYDT =  StomCond_Mb(t,StomCond_Con,StomCond_Param);

StomPS_DYDT = zeros(10,1);

for index = 1:120        
    StomPS_DYDT(index) = trDynaPS_DYDT(index);
end

for index = 1:2
    StomPS_DYDT(index+130) = StomCond_DYDT(index);
end

 StomPS_DYDT(63) = StomCond_DYDT(2); 

GenOut(t); 

DONE = 1; 