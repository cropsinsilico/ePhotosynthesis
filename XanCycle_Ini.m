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

function [XanCycle_Con] = XanCycle_Ini

global XanCycle_kva;
global XanCycle_kaz;
global XanCycle_kza;
global XanCycle_kav;
global XanCycle_kvf;
global XanCycle_kv2ABA;
global XanCycle_kABAdg;

XanCycle_kva = 0.163/60;     % Ruth Frommolt et a; 2001; Planta
XanCycle_kaz = 0.691/60;     % Ruth Frommolt et a; 2001; Planta
XanCycle_kza = 0.119/60;     % Ruth Frommolt et a; 2001; Planta
XanCycle_kav = 0.119/60;     % Ruth Frommolt et a; 2001; Planta. This is not given in the paper. Therefore, teh value is really an educated guess. 

XanCycle_kvf = 0;            % This is the rate of formation of v from its precursors, reprenting the net generation of new V.
XanCycle_kv2ABA = 0;         % This represent the rate constant of conversion from v to ABA. This is just a guess. 
XanCycle_kABAdg = 0;         % This represent the rate constant of conversion ABA degradation

Vx = 160;
Ax = 10;
Zx = 5; 
ABA = 1;   
XanCycle_Con = zeros(4,1);
Coeff = 0.37;

XanCycle_Con(1) = Vx * 0.37;
XanCycle_Con(2) = Ax  * 0.37;
XanCycle_Con(3) = Zx  * 0.37;
XanCycle_Con(4) = ABA;

global XanCycle_OLD_TIME;
global XanCycle_TIME_N;
global XanCycle_VEL;
global XanCycle_CON;

XanCycle_OLD_TIME = 0;
XanCycle_TIME_N = 1;

XanCycle_VEL = zeros(1,4);    % Clean memory
XanCycle_CON = zeros(1,4);    % Clean memory


global XanCycle2FIBF_Xstate;
XanCycle2FIBF_Xstate = Zx/(Ax + Vx + Zx);

