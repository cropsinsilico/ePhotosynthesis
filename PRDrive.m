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

Begin = 1;
fin = SYSInitial(Begin);
global options1;
global tglobal;
time = tglobal;



PrS = PRinitial(Begin);
backsim = 1;


global PR_OLD_TIME;
global PR_TIME_N;
global PR_VEL;
global PR_CON;

PR_OLD_TIME = 0;
PR_TIME_N = 1;

PR_VEL = zeros(1,11);  
PR_CON = zeros(1,1);    

ModelComb = IniModelCom;        

if (backsim == 1)
    PR_Param = zeros(2,1);
    PR_Param(1) = 1;                       
    [Tt,d] = ode15s(@PRmb,[0,time],PrS,options1,PR_Param);
    
    for m = 1:6
        subplot(2,3,m);p = plot(Tt,d(:,m),'.');suc = PRespTitle(m,p,1);xlabel('time (s)');ylabel('mM');
    end
    pause;

    for m = 5:10
        subplot(2,3,m-4);p = plot(Tt,d(:,m),'.');suc = PRespTitle(m,p,1);xlabel('time (s)');ylabel('mM');
    end
    pause;
    
    % Plot the information about flux:
    % Since the flux were calculated using unit: mM/s; it is converted into micro mole per meter square per second.
    global AVR; 
    
    VCOEFF = AVR;
    VCOEFF2 = AVR;
    
    PR_VEL(:,2) = PR_VEL(:,2)*VCOEFF;
    PR_VEL(:,3) = PR_VEL(:,3)*VCOEFF;
    PR_VEL(:,4) = PR_VEL(:,4)*VCOEFF;
    PR_VEL(:,5) = PR_VEL(:,5)*VCOEFF2;
    PR_VEL(:,6) = PR_VEL(:,6)*VCOEFF2;
    PR_VEL(:,7) = PR_VEL(:,7)*VCOEFF2;
    PR_VEL(:,8) = PR_VEL(:,8)*VCOEFF2;
    PR_VEL(:,9) = PR_VEL(:,9)*VCOEFF2;
    PR_VEL(:,10) = PR_VEL(:,10)*VCOEFF2;
    PR_VEL(:,11) = PR_VEL(:,11)*VCOEFF2;
    
    for m = 1:6
        subplot(2,3,m); p = plot(PR_VEL(:,1),PR_VEL(:,m+1),'.'); suc = PRespTitle(m,p,2);xlabel('time (s)');ylabel('uM m-2 s-1');
    end
    pause;
    
    for m = 5:10
        subplot(2,3,m-4); p = plot(PR_VEL(:,1),PR_VEL(:,m+1),'.'); suc = PRespTitle(m,p,2);xlabel('time (s)');ylabel('uM m-2 s-1');
    end
    pause;
    
else

    d = zeros(10,100);
    
    for t = 1:time
        DSDT = PRrate(t,PrS);
        PrS = PrS + DSDT;
        d(:,t)= PrS;
    end

    for m = 1:10
        subplot(2,5,m);plot(1:100,d(m,:),'.');
    end
end