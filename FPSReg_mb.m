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



function FPSReg_DYDT = FPSReg_mb(t, FPSReg_Con, BF_Param, FI_Param, PS_PR_Param, RuACT_Param, RedoxReg_Param, XanCycle_Param, RROEA_Param, SUCS_Param)

for m = 1:120       % TEST 104
    trDynaPS_Con(m) = FPSReg_Con (m); 
end

for m = 1:12
    SUCS_Con(m) = FPSReg_Con(m + 120);
end

light = Condition (t);

FI_Param(1) = light;
BF_Param(1) = light;

SUCS_DYDT = SUCS_mb(t, SUCS_Con, SUCS_Param);

trDynaPS_DYDT = trDynaPS_mb(t, trDynaPS_Con, BF_Param, FI_Param, PS_PR_Param, RuACT_Param, RedoxReg_Param, XanCycle_Param, RROEA_Param);


FPSReg_DYDT = zeros(10,1);

for index = 1:120
    FPSReg_DYDT(index) = trDynaPS_DYDT(index);  
end

for index = 1:12
    FPSReg_DYDT(index+120) = SUCS_DYDT(index); 
end

global PS2CM_vdhap;
vdhap = PS2CM_vdhap;

global PS2CM_vgap;
vgap = PS2CM_vgap;

global SUCS2CM_vdhap;
global SUCS2CM_vgap;
vdhap_ins = SUCS2CM_vdhap		; 
vgap_ins = SUCS2CM_vgap			;   

SUCS_DYDT(1)	=	SUCS_DYDT(1) + vdhap + vgap -(vdhap_ins + vgap_ins) 		;   
FPSReg_DYDT(121) = SUCS_DYDT(1);

global PS2CM_vpga;
vpga = PS2CM_vpga;
global SUCS2CM_vpga;
vpga_ins = SUCS2CM_vpga	;    

SUCS_mb(12)	=	SUCS_DYDT(12) - vpga_ins + vpga 		;   %	pgaC
FPSReg_DYDT(132) = SUCS_DYDT(12);


