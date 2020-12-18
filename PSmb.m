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

function PSdydt = PSmb(t,PSs, Param)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Modifying KM, KI, KE VMAX for different reactions as the regulation%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Regulations first. 
fini = Condition (t);

% Get the rate for the reactions in the photosynthesis sytem

PSr = PSRate(t,PSs,Param);

% Get the rate

v1	=	PSr(1)	;
v2	=	PSr(2)	;
v3	=	PSr(3)	;
NONE=	PSr(4)	;
v5	=	PSr(5)	;
v6	=	PSr(6)	;
v7	=	PSr(7)	;
v8	=	PSr(8)	;
v9	=	PSr(9)	;
v10	=	PSr(10)	;
v13	=	PSr(11)	;
v16	=	PSr(12)	;
v23	=	PSr(13)	;
v31	=	PSr(14)	;
v32	=	PSr(15)	;
v33	=	PSr(16)	;
v24 =   PSr ( 17);
v25 =   PSr (18); 

% Implement the mass balance equations

tmp = zeros(15,1);


tmp(1) = v13-v1;


tmp(2) = 2*v1-v2-v32;


tmp(3) = v2-v3;


tmp(4) = v3 - 2 * v5 - v7 - v8 - v10 - v31 - v33;


tmp(5) = v23 - v24;


tmp(6) = v5-v6;


tmp(7) = v7-v8;


tmp(8) = v9-v10;


tmp(9) = v8 - v9;


tmp(10) = v16 - v2 - v23 - v13 - v25;           

global PS2EPS_NADPH;
global PS_C_CN;

tmp(11) = 0;



tmp(12) = 0;


tmp(13) = 0;


tmp(14) = v6 - v7 - v23 + v25;


tmp(15) = v7 + v10 * 2 - v13;

PSdydt = zeros(15,1);
PSdydt = tmp;