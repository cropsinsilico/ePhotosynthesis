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

function DynaPS_Con = DynaPS_Ini

BEGIN = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Clear up memory for simulation       %  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global DynaPS_OLD_TIME;
global DynaPS_TIME_N;
global DynaPS_VEL;
global DynaPS_CON;

DynaPS_OLD_TIME = 0;
DynaPS_TIME_N = 1;
DynaPS_VEL = zeros(1,3);    % Clean memory
DynaPS_CON = zeros(3,1);    % Clean memory


RA_Con = RA_Ini(BEGIN);
DynaPS_Con = zeros(5,1);

for m = 1:92
    DynaPS_Con(m) = RA_Con (m);
end


XanCycle_Con = XanCycle_Ini;

for m = 1:4
    DynaPS_Con(m+92) =  XanCycle_Con (m);
end
