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


function [FI_Param, BF_Param, PS_PR_Param, SUCS_Param, EPS_Param, RuACT_Param, XanCycle_Param, RROEA_Param, RedoxReg_Param, StomCond_Param] = ParamSet

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

SUCS_Param = zeros(2,1);
SUCS_Param = 0;

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

StomCond_Param = zeros(2,1);
StomCond_Param = 1;
