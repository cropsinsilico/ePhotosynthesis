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


% DynaPS_mb.m  This model includes the mass balance equations for the full model of photosynthesis.

function DynaPS_DYDT = DynaPS_mb(t, DynaPS_Con, BF_Param, FI_Param, PS_PR_Param, RuACT_Param, RedoxReg_Param, XanCycle_Param, SUCS_Param)

% Try out one new way of calculating the mass balance equation.
% In this new way, all the previous calcuations of mass balance equation is preserved and only the necessary changes are made.

%% Step One: Get the initialization of the concentrations for the RedoxReg model which will be used in the calculation of mb of RedoxReg.

for m = 1:92
    RA_Con(m) = DynaPS_Con (m);
end


for m = 1:4
    XanCycle_Con(m) = DynaPS_Con(m + 92);
end

% This is a sensitivity test to show that the model is stable udner fluctuating light

light = Condition (t);

FI_Param(1) = light;
BF_Param(1) = light;

RA_DYDT = RA_mb(t, RA_Con, BF_Param, FI_Param, PS_PR_Param, RuACT_Param, SUCS_Param);
XanCycle_DYDT = XanCycle_mb(t, XanCycle_Con, XanCycle_Param);

% Here get the rate of Thioredoxin reduction and oxidation and use it to construct the differential equation for both thio and fd. 

DynaPS_DYDT = zeros(10,1);

for index = 1:92
    DynaPS_DYDT(index) = RA_DYDT(index);
end

for index = 1:4
    DynaPS_DYDT(index+92) = XanCycle_DYDT(index);
end