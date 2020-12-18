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




function RA_DYDT = RA_mb(t, RA_Con, BF_Param, FI_Param, PS_PR_Param, RuACT_Param, SUCS_Param)

for m = 1:88
    EPS_Con(m) = RA_Con(m);
end

for m = 1:4
    RuACT_Con(m) = RA_Con(m+88);
end

EPS_DYDT = EPS_mb(t,EPS_Con,BF_Param, FI_Param, PS_PR_Param, SUCS_Param);
RuACT_DYDT = RuACT_MB(t, RuACT_Con, RuACT_Param);

RA_DYDT = zeros(92,1);

for m = 1:88
    RA_DYDT(m) = EPS_DYDT(m);
end

for m = 1:4
    RA_DYDT(m+88) = RuACT_DYDT(m);
end

global PSPR2RA_v1;
global PSPR2RA_v13;
global PSPR2RA_v111;


global RuACT2RA_v61;
global RuACT2RA_v62;
global RuACT2RA_v1;
global RuACT2RA_vn1;
global RuACT2RA_v7;
global RuACT2RA_vn7;

DYDT_RuBP =  RuACT2RA_v1 + PSPR2RA_v13 - RuACT2RA_vn1 + RuACT2RA_vn7 - RuACT2RA_v7;  
 RA_DYDT(53) = DYDT_RuBP;
 RA_DYDT(92) = DYDT_RuBP;

DYDT_PGA = EPS_DYDT(54) - 2 * PSPR2RA_v1 + 2 * RuACT2RA_v61 - PSPR2RA_v111 + RuACT2RA_v62;      % Originally it is pspr(2), now use EPS_DYDT(54).
 RA_DYDT(54) = DYDT_PGA;


DYDT_PGCA = EPS_DYDT(69) - PSPR2RA_v111 + RuACT2RA_v62;           
RA_DYDT(69) = DYDT_PGCA;