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



% EPS_mb.m  This model includes the mass balance equations for the full model of the light reactions.

function EPS_DYDT = EPS_mb(t, EPS_Con, BF_Param, FI_Param, PS_PR_Param, Sucs_Param)

% Try out one new way of calculating the mass balance equation.
% In this new way, all the previous calcuations of mass balance equation is preserved and only the necessary changes are made.

%% Step One: Get the initialization of the concentrations for the PSPR model which will be used in the calculation of mb of CM.

for m = 1:52
    FIBF_Con(m) = EPS_Con(m);
end

for m = 1:35
    CMs(m) = EPS_Con(m+52);
end

% This is a sensitivity test to show that the model is stable udner fluctuating light
% the condition are is distributed into the separate mb file. 

% Step II, Calculate the PYPT using the existing modules.
CM_DYDT = CM_mb(t,CMs,PS_PR_Param, Sucs_Param);
FIBF_DYDT = FIBF_MB(t, FIBF_Con, BF_Param, FI_Param);

% Dump a copy of the model parameters if export_mod_enabled is set
global export_mod_enabled;
if (export_mod_enabled)
    if (t == 0)
        dump_vars("EPS_init.txt", EPS_Dict(EPS_Con));
    else
        EPS_Con_dump = zeros(87,1);
        for m = 1:52
            EPS_Con_dump(m) = FIBF_Con(m);
        end
        for m = 1:35
            EPS_Con_dump(m+52) = CMs(m);
        end
        dump_vars("EPS_final.txt", EPS_Dict(EPS_Con_dump));
    end
    export_mod_enabled = 0;
end

% Step III: Calculate the mass balanec equation for the EPS model. This basically need to make sure that the variables 
% used in the mass balance equation should be in exact sequence with the sequence used in the inialization.

EPS_DYDT = zeros(87,1);
for m = 1:52
    EPS_DYDT(m) = FIBF_DYDT(m);
end

for m = 1:35
    EPS_DYDT(m+52) = CM_DYDT(m);
end

global EPS_ATP_Rate;   % The EPS_ATP_Rate is used in the overall model for the calculation of the mass balance equation of ATP.
global PS2EPS_V16;     

EPS_DYDT(61) = CM_DYDT(9) - PS2EPS_V16 + EPS_ATP_Rate;
EPS_DYDT(17) = EPS_DYDT(61);

global PS2EPS_v3;
global BF2EPS_vbfn2;
global PS2EPS_NADPH;

global BF_RC; 
Vmax11 = BF_RC	(	11	)			;	%	The maximum rate of ATP synthesis	Unit: mmol l-1 s-1; The unit for the reactions occurrs in stroma is mmol l-1 s-1

%EPS_DYDT(62) = BF2EPS_vbfn2 - PS2EPS_v3 - 1 * PS2EPS_NADPH/(PS2EPS_NADPH + 0.5) ;
EPS_DYDT(62) = BF2EPS_vbfn2/2 - PS2EPS_v3;% - 1 * PS2EPS_NADPH/(PS2EPS_NADPH + 0.5) ;  %QF changed /2 and ;% - 1 * PS2EPS_NADPH/(PS2EPS_NADPH + 0.5)
EPS_DYDT(29) = EPS_DYDT(62);           
