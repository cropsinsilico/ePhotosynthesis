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


% PSPR_Out.m    This program provide the output for the model of the PS and PR.

function CA_PSPR_Out(Tt, d)

global PS_VEL;
global PR_VEL;
global CarbonRate;
global size;
global Metabolites;

Metabolites = d;
global AVR;        
 

    PR_VEL(:,2) = PR_VEL(:,2)*AVR;
    PR_VEL(:,3) = PR_VEL(:,3)*AVR;
    PR_VEL(:,4) = PR_VEL(:,4)*AVR;
    PR_VEL(:,5) = PR_VEL(:,5)*AVR;
    PR_VEL(:,6) = PR_VEL(:,6)*AVR;
    PR_VEL(:,7) = PR_VEL(:,7)*AVR;
    PR_VEL(:,8) = PR_VEL(:,8)*AVR;
    PR_VEL(:,9) = PR_VEL(:,9)*AVR;
    PR_VEL(:,10) = PR_VEL(:,10)*AVR;
    PR_VEL(:,11) = PR_VEL(:,11)*AVR;

% Get the fluxes for reactions in photosynthesis pathway

PS_VEL = PS_VEL * AVR;
PS_VEL(1,:)=PS_VEL(1,:)/(AVR);

% The following plot the net carbon assimilation rate based on PS_VEL and PR_VEL.

CarbonRate = PS_VEL(2,:);
CO2Release = PR_VEL(:,9);
NetCarbonT = CarbonRate' - CO2Release;
LengthT  = length(NetCarbonT);
CarbonRate = NetCarbonT(LengthT);


PSVEL = PS_VEL';
PRVEL = PR_VEL';

global CO2A;

