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

Begin = 1;
 
fin = SYSInitial(Begin);
global options1;
global tglobal;
time = tglobal;

ModelComb = IniModelCom;      

global PSPR_SUCS_com;
PSPR_SUCS_com = 0;  
SUCS_Con = SUCS_Ini; 
va1 = 1000; 
n = 1;  

SUCS_Param = zeros(2,1);
SUCS_Param(1) = va1;
SUCS_Param(2) = n;

global SUCS_VEL;

[Tt,d] = ode15s(@SUCS_MB,[0,time],SUCS_Con,options1,SUCS_Param);
 

done = SUCS_Graph(Tt, d);
clock;

IniModelCom;

AOUT = SUCS_VEL';