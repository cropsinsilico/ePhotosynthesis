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

% FIBF_Out.m    This model is the part that will be responsible for output of result.

function FinishOut = FIBF_Out(Tt, d)

% Assign the concentrations from d to the concentration variable in BF and FI.

global FI_VEL;
global FI_CON;
global BF_VEL;
global BF_CON;

global FI;
global BF;


    M = size(d);
    a = M(1);
    
    BF = zeros(a,29);
    FI = zeros(a,22);
    
    for m = 1:29
         BF(:,m)= d(:,m);
    end
    
    for m = 30:51
         FI(:,m-29)=d(:,m);
    end
    
    for m = 1:10
        subplot(2,5,m);p = plot((Tt),FI(:,m),'.');ylabel('micro mol/m2');xlabel('log 10 second');suc = FI_AddTitle(m,p,1);
    end
    pause;

    for m = 11:20
        subplot(2,5,m-10);p = plot((Tt),FI(:,m),'.');ylabel('micro mol/m2');xlabel('log 10 second');suc = FI_AddTitle(m,p,1);
    end

    pause;

    for m = 1:10
        subplot(2,5,m);p = plot((FI_VEL(:,1)),FI_VEL(:,m+1),'.');ylabel('micro mol/m2/s');xlabel('log 10 second');suc = FI_AddTitle(m,p,2);
    end
    pause;

    for m = 11:20
        subplot(2,5,m-10);p = plot((FI_VEL(:,1)),FI_VEL(:,m+1),'.');ylabel('micro mol/m2/s');xlabel('log 10 second');suc = FI_AddTitle(m,p,2);
    end
    pause;

    for m = 21:30
        subplot(2,5,m-20);p = plot((FI_VEL(:,1)),FI_VEL(:,m+1),'.');ylabel('micro mol/m2/s');xlabel('log 10 second');suc = FI_AddTitle(m,p,2);
    end
    pause;

    for m = 31:40
        subplot(2,5,m-30);p = plot((FI_VEL(:,1)),FI_VEL(:,m+1),'.');ylabel('micro mol/m2/s');xlabel('log 10 second');suc = FI_AddTitle(m,p,2);
    end
    pause;
    

    for m = 41:50
        subplot(2,5,m-40);p = plot((FI_VEL(:,1)),FI_VEL(:,m+1),'.');ylabel('micro mol/m2/s');xlabel('log 10 second');suc = FI_AddTitle(m,p,2);
    end
    pause;
    
    for m = 51:53
        subplot(1,3,m-50);p = plot((FI_VEL(:,1)),FI_VEL(:,m+1),'.');ylabel('micro mol/m2/s');xlabel('log 10 second');suc = FI_AddTitle(m,p,2);
    end
    pause;
    
    for m = 1:5
        subplot(1,5,m);p = plot((FI_VEL(:,1)),FI_CON(:,m+1),'.');ylabel('micro mol/m2/s');xlabel('log 10 second');suc = FI_AddTitle(m,p,3);
    end
    pause;

    for m = 1:10
        subplot(2,5,m);p = plot((Tt),BF(:,m),'.');ylabel('micro mol/m2');xlabel('second');suc = BF_AddTitle(m,p,1);
    end
    pause;

    for m = 11:14
        subplot(1,4,m-10);p = plot((Tt),BF(:,m),'.');ylabel('micro mol/m2');xlabel(' second');suc = BF_AddTitle(m,p,1);
    end
    pause;
   
    for m = 16:20
        subplot(1,5,m-15);p = plot((Tt),BF(:,m),'.');ylabel('mmol/l');xlabel(' second');suc = BF_AddTitle(m,p,1);
    end
    pause;

    for m = 21:24
        subplot(2,2,m-20);p = plot((Tt),BF(:,m),'.');ylabel('micromol/m2');xlabel('second');suc = BF_AddTitle(m,p,1);
    end
    pause;
    
    for m = 25:26
        subplot(1,2,m-24);p = plot((Tt),BF(:,m),'.');ylabel('milimole/l');xlabel('second');suc = BF_AddTitle(m,p,1);
    end
    pause;
    
    for m = 27:28
        subplot(1,2,m-26);p = plot((Tt),BF(:,m),'.');ylabel('PH');xlabel('second');suc = BF_AddTitle(m,p,1);
    end
    pause;

    newplot;
    p = plot((Tt),d(:,29),'.');ylabel('mmol l-1');xlabel('second');suc = BF_AddTitle(29,p,1);
    pause;
    
    % Output the Velocity
    for m = 1:10
        subplot(2,5,m);p = plot((BF_VEL(:,1)),BF_VEL(:,m+1),'.');ylabel('micro mol/m2/s');xlabel('second');suc = BF_AddTitle(m,p,2);
    end
    pause;
    
    newplot;
    p = plot((BF_VEL(:,1)),BF_VEL(:,12),'.');ylabel('mmol/l/s');xlabel('second');suc = BF_AddTitle(11,p,2);
    pause;
    
    for m = 12:20
        subplot(2,5,m-11);p = plot((BF_VEL(:,1)),BF_VEL(:,m+1),'.');ylabel('micro mol/m2/s');xlabel('second');suc = BF_AddTitle(m,p,2);
    end
    pause;

    for m = 21:22
        subplot(2,4,m-20);p = plot((BF_VEL(:,1)),BF_VEL(:,m+1),'.');ylabel('mmol/l/s');xlabel('second');suc = BF_AddTitle(m,p,2);
    end
    pause;

    newplot;
    p = plot((BF_VEL(:,1)),BF_VEL(:,23),'.');ylabel('mircomol/m2/s');xlabel('second');suc = BF_AddTitle(23,p,2);
    pause;
   
    for m = 24:26
        subplot(1,3,m-23);p = plot((BF_VEL(:,1)),BF_VEL(:,m+1),'.');ylabel('mmol/l/s');xlabel('second');suc = BF_AddTitle(m,p,2);
    end
    pause;
    
    for m = 27:28
        subplot(1,2,m-26);p = plot((BF_VEL(:,1)),BF_VEL(:,m+1),'.');ylabel('micromol /m2 leaf/s');xlabel('second');suc = BF_AddTitle(m,p,2);
    end
    pause;

    newplot;
    p = plot((BF_VEL(:,1)),BF_VEL(:,30),'.');ylabel('mmol/l/s');xlabel('second');suc = BF_AddTitle(29,p,2);
    pause;
    
    newplot;
    FinishOut = 1;
    