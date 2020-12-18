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

% function DynaPS_Graph This function draw the graphs for DynaPS

function done = DynaPS_Graph (Tt, d)
M = size(d);
    a = M(1);
    
    FIBF = zeros(a,52);
    PSPR = zeros(a,23);
    
    for m = 1:52
         FIBF(:,m)= d(:,m);
    end
    
    for m = 1:24
         PSPR(:,m)=d(:,m+52);
    end
    
    global FI;
    global BF;
    
    global PS;
    global PR;
    
    OutPutLight = FIBF_Out(Tt, FIBF);

    OutputPSPR = PSPR_Out(Tt, PSPR);
    
    % output the rate of CO2 UTPAKE.
    global AVR;         % The conversion factor between micro mole per meter square per second and milimole per liter per second. 

    global RuACT_VEL;
    global PR_VEL;
    
    temp = RuACT_VEL(:,6) ;
    CarbonRate = temp * AVR; 
    
    CO2Release = PR_VEL(:,9) * AVR; 
    Assim = CarbonRate - CO2Release;
    
    newplot;
    p = plot((PR_VEL(:,1)), Assim,'.');ylabel('Assimilation Rate');xlabel('second');title('Net Assimilation Rate');
    pause;

    newplot;
    p = plot((PR_VEL(:,1)), CarbonRate,'.');ylabel('Assimilation Rate');xlabel('second');title('Gross Assimilation Rate');
    pause;
    
    % SUCS
    SUCS = d(:,77:88);
    done = SUCS_Graph(Tt, SUCS);    
    
    % RuACT
    
    d_RuACT= d(:,89:92);
    done = RuACT_Graph(Tt, d_RuACT);
    
    % Xanthophyll cycle
     for m = 1:4
        subplot(1,4,m);p = plot(Tt,d(:,m+92),'.');ylabel('mM');xlabel(' second');suc = XanCycle_AddTitle(m,p,1);
    end
        done = 1;