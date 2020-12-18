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



function DynaPS_Con = trDynaPS_Ini

BEGIN = 1;

global trDynaPS_OLD_TIME;
global trDynaPS_TIME_N;
global trDynaPS_VEL;
global trDynaPS_CON;

trDynaPS_OLD_TIME = 0;
trDynaPS_TIME_N = 1;
trDynaPS_VEL = zeros(1,3);    % Clean memory
trDynaPS_CON = zeros(3,1);    % Clean memory

% Now get the combined total concentration of different concentration variables. 
DynaPS_Con = DynaPS_Ini;
trDynaPS_Con = zeros(5,1);

for m = 1:96
    trDynaPS_Con(m) = DynaPS_Con (m);
end

RROEA_Con = RROEA_Ini;

for m = 1:10
    DynaPS_Con(m+110) = RROEA_Con (m);
end
