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

% SYSInitial.m
% This routine initialze the parameters used for all the routines.

function fin = SYSInitial(Begin)

% The total runing time is set as tglobal
global tglobal;
global Primer; 
Primer = 300; 

tglobal = 1000; 

 
global options1
options1 = odeset('RelTol',1e-4);

global PS12ratio;        
global input_PSIIcore;
global input_PSI;
global PSIIantennaSize;
global PSIantennaSize;
global input_LHCII;
global input_LHCI; 

input_PSIIcore=1;
input_PSI=1;
PS12ratio = input_PSI/input_PSIIcore;

PSIIantennaSize = 37;
PSIantennaSize = 95;
input_LHCII=13;
input_LHCI=6; 
global ChlT2;
global ChlT;
global ChlPSI;
ChlT2 = input_PSIIcore * (PSIIantennaSize + 13 * input_LHCII);  % U and A, PSII and LHCII
ChlT  = PSIIantennaSize * input_PSIIcore; % U , PSII
ChlPSI = input_PSI * (PSIantennaSize + 13 * input_LHCI);  % U and A of PSI, total Chl in PSI


global AVR;         
AVR = 30; 

global GP; 
GP = 0; 

fin = 1; 
