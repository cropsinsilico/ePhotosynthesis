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



function ModelComb = IModelCom;

global RuACT_EPS_com;
RuACT_EPS_com = 0;

global BF_FI_com;            
BF_FI_com = 0;

global PR_PS_com;   
PR_PS_com = 0;

global FIBF_PSPR_com;  
FIBF_PSPR_com = 0;    

global ATPActive;
ATPActive = 0;

global RedoxReg_RA_com;
RedoxReg_RA_com = 0;

global XanCycle_BF_com;
XanCycle_BF_com = 0;

global RROEA_EPS_com;
RROEA_EPS_com = 0;

global StomCond_TrDynaPS_com;
StomCond_TrDynaPS_com = 0;

global PSPR_SUCS_com;
PSPR_SUCS_com = 0;  

global trDynaPS_SUCS_com;
trDynaPS_SUCS_com = 0;

global EPS_SUCS_com;
EPS_SUCS_com = 0;

ModelComb = 1;