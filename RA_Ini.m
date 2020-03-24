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


function RA_Con = RA_Ini(BEGIN)

BEGIN = 1;
EPS_Con = EPS_Ini(BEGIN);
RuACT_Con = RuACT_Ini(BEGIN);


global RuACT_OLD_TIME;
global RuACT_TIME_N;
global RuACT_VEL;
global RuACT_CON;

RuACT_OLD_TIME = 0;
RuACT_TIME_N = 1;

RuACT_VEL = zeros(1,3);    % Clean memory
RuACT_CON = zeros(3,1);    % Clean memory

% Now get the combined total concentration of different concentration variables. 

RA_Con = zeros(5,1);

for m = 1:88
    RA_Con (m) = EPS_Con(m);
end

for m = 1:4
    RA_Con(m+88) = RuACT_Con(m);         
end