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


function FPSReg_Con = FPSReg_Ini

BEGIN = 1;


global FPSReg_OLD_TIME;
global FPSReg_TIME_N;
global FPSReg_VEL;
global FPSReg_CON;

FPSReg_OLD_TIME = 0;
FPSReg_TIME_N = 1;
FPSReg_VEL = zeros(1,3);    % Clean memory
FPSReg_CON = zeros(3,1);    % Clean memory


trDynaPS_Con = trDynaPS_Ini;
FPSReg_Con = zeros(5,1);

for m = 1:120
    FPSReg_Con(m) = trDynaPS_Con (m);
end

Begin = 1;
SUCS_Con = SUCS_Ini(Begin);

for m = 1:12
    FPSReg_Con(m+120) = SUCS_Con (m);
end
