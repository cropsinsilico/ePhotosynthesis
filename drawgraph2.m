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

function AC = drawgraph2 (d, Tt)

%subplot(2,2,1);
%plot(PS_VEL(1,:),1.5-CO2A(:,1));
%title('ATP vs Time')

%subplot(2,2,2);
%plot(PS_VEL(1,:),PR_VEL(:,4));
%title('Glycerate kinase vs Time')

%subplot(2,2,3);
%plot(PS_VEL(1,:),CO2A(:,3));
%title('CO2 vs T')

%subplot(2,2,4);
%plot(PS_VEL(1,:),CO2A(:,4));
%title('A vs T')

%pause;
%newplot;
%subplot(2,2,1);
%plot(PS_VEL(1,:),PS_VEL(3,:));
%title('PGA kinase vs Time')

%subplot(2,2,2);
%plot(PS_VEL(1,:),PS_VEL(12,:));
%title('PRK activase vs Time')

%subplot(2,2,3);
%plot(PS_VEL(1,:),PS_VEL(14,:));
%title('Starch synthesis vs Time')

%subplot(2,2,4);
%plot(PS_VEL(1,:),PS_VEL(13,:));
%title('ATP synthesis vs Time')

%pause;
%hold on;

% M = PR_VEL(:,1);

TotalTime =12000;
TimeInterval = 1200;

M = Tt;

Total = TotalTime/TimeInterval;

global ACIArray;
ACIArray = zeros(2,5);

for index = 1: (Total-1)
    A = find(M > TimeInterval * index);
    D = A(1);
    B = D -2;
    ACIArray(1,index) = d(B, 20);
    ACIArray(2,index) = d(B, 21);
    ACIArray(3,index) = 0.04 * index;
        
end


hold off;
A = size(d);
RowNumber = A(1);
ACIArray(1,Total) = d(RowNumber, 20);
ACIArray(2,Total) = d(RowNumber, 21);
ACIArray(3,Total) = 0.04*Total;

%newplot;
subplot(2,2,1);
plot(ACIArray(3,:),ACIArray(1,:),'ro');
title('Serine');
xlabel('Percentage Enzyme Activity(mM s-1)');ylabel('Serine Concentration');
hold off;

subplot(2,2,2);
plot(ACIArray(3,:),ACIArray(2,:),'ro');
title('Glycine');
xlabel('Percentage Enzyme Activity(mM s-1)');ylabel('Glycine Concentration');
hold off;

AC = ACIArray';


%subplot(2,2,2);
%plot(ACIArray(1,:),ACIArray(3,:),'ro');
%title('ATP vs Activity');
%xlabel('Enzyme Activity(mM s-1)');ylabel('[ATP] mM');
%hold off;

done = 1; %ACIArray';
%pause;    % Do the light curve. 
%M = PR_VEL(:,1);
%Total = TotalTime/TimeInterval;
%LightArray = zeros(2,5);

%for index = 1: (Total-1)
%    A = find(M > TimeInterval * index);
%    D = A(1);
%    B = D-2;
%    LightArray(1,index) = CO2A(B, 6);
%    LightArray(2,index) = CO2A(B, 4);
%    LightArray(3,index) = 1.5 - CO2A(B, 1);     % Plot ATP
%    
%end


%A = size(CO2A);
%RowNumber = A(1);
%LightArray(1,Total) = CO2A(RowNumber, 6);
%LightArray(2,Total) = CO2A(RowNumber, 4);
%LightArray(3,Total) = 1.5 - CO2A(RowNumber, 1);     % Plot ATP

%newplot;
%subplot(2,2,3);
%plot(LightArray(1,:),LightArray(2,:),'bo');
%title('AQ curve');
%xlabel('Q(microE)');ylabel('A micromole s-1 s-1');


%subplot(2,2,4);
%plot(LightArray(1,:),LightArray(3,:),'bo');
%title('ATP vs Light');
%xlabel('Q(microE)');ylabel('[ATP] mM');

%pause;

%subplot(2,2,3);
%plot(ACIArray(2,:),ACIArray(3,:),'ro');
%title('A vs [ATP]'); xlabel('A (micomole m-2 s-1)');ylabel('ATP (mM)');
%hold on;