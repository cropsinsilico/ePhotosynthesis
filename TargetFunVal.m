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

function FunV = TargetFunVal;

global PS_VEL;
global PR_VEL;

PSVCOEFF = 30;
PS_VEL = PS_VEL * PSVCOEFF; 

PS_VEL(1,:)=PS_VEL(1,:)/(PSVCOEFF);     % The time;
PSVEL = PS_VEL';

global VolRatioStCyto
if VolRatioStCyto ==1
    PR_VEL = PR_VEL * PSVCOEFF;
else
    PR_VEL = PR_VEL * PSVCOEFF*4/9;
end
    

n = size(PS_VEL);
a = PS_VEL(2,n(2)); 
b = PR_VEL(n(2),9); 

CO2AR = a - b;

FunV = CO2AR;