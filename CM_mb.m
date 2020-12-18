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


function CM_DYDT = CM_mb(t, CM_Con, PS_PR_Param, SUCS_Param)

for m = 1:23
    PSPR_Con(m) = CM_Con(m);
end

for m = 1:12
    SUCS_Con(m) = CM_Con(23 + m);
end

PSPR_Con(24) = CM_Con(36);


SUCS_DYDT = SUCS_mb(t, SUCS_Con, SUCS_Param);
PSPR_DYDT = PS_PRmb(t,PSPR_Con,PS_PR_Param);

CM_DYDT = zeros(36,1);
for m = 1:23
    CM_DYDT(m) = PSPR_DYDT(m);
end

for m = 1:12
    CM_DYDT(m+23) = SUCS_DYDT(m);
end

CM_DYDT(36) = PSPR_DYDT(24);

global PS2CM_vdhap;
vdhap = PS2CM_vdhap;        % The rate of export out of chloroplast

global PS2CM_vgap;          % The rate of export out of chloroplast
vgap = PS2CM_vgap;

global SUCS2CM_vdhap;       % The rate of import into the cytosol
global SUCS2CM_vgap;        % The rate of import into the cytosol
vdhap_ins = SUCS2CM_vdhap		    ;   %	DHAP IN
vgap_ins  = SUCS2CM_vgap			;   %	GAP IN

SUCS_DYDT(1)   =	SUCS_DYDT(1) + vdhap + vgap -(vdhap_ins + vgap_ins) 		;   %	T3Pc 
CM_DYDT(24) = SUCS_DYDT(1);   


global PS2CM_vpga;
vpga = PS2CM_vpga;

global SUCS2CM_vpga;
vpga_ins = SUCS2CM_vpga	;                                       %	PGA export from chloroplast 

SUCS_DYDT(12)	=	SUCS_DYDT(12) - vpga_ins + vpga 		;   %	pgaC
CM_DYDT(35)     =   SUCS_DYDT(12);