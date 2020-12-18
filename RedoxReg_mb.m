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





function RedoxReg_DYDT = RedoxReg_mb(t, RedoxReg_Con, BF_Param, FI_Param, PS_PR_Param, RuACT_Param, RedoxReg_Param, SUCS_Param)


global trDynaPS2RedReg_cal
trDynaPS2RedReg_cal = 1;

for m = 1:92
    RA_Con(m) = RedoxReg_Con(m);
end

ThioRe = RedoxReg_Con(93);

RedoxReg_Vel = RedoxReg_Rate (t, RedoxReg_Con, RedoxReg_Param);

RA_DYDT = RA_mb(t, RA_Con, BF_Param, FI_Param, PS_PR_Param, RuACT_Param, SUCS_Param);

RedoxReg_DYDT = zeros(10,1);

for index = 1:92
    RedoxReg_DYDT(index) = RA_DYDT(index);
end

vred = RedoxReg_Vel(1);
vox  = RedoxReg_Vel(2);

RedoxReg_DYDT(93) = vred - vox;
RedoxReg_DYDT(93) = 0;


Temp = RedoxReg_DYDT(24);
RedoxReg_DYDT(24) = Temp;  