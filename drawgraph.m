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


% function done = drawgraph (PS_VEL, PR_VEL, CO2A)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     These two variables were taken from condition and SystemInfor
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global tglobal; 
global NumInter_draw; 

TotalTime = tglobal;
IntervalNumber = NumInter_draw;
IntervalNumber = 10; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Total = IntervalNumber;
TimeInterval = TotalTime/IntervalNumber;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global PR_VEL;

% M = PS_VEL(1,:);
M = CO2A(:,100);

Total = TotalTime/TimeInterval;

global ACIArray;
ACIArray = zeros(2,2);

for index = 1: (Total-1)
    A = find(M > TimeInterval * index);
    D = A(1);
    B = D -2;
    for index2 = 1:54       % 47
        ACIArray(index2,index) = CO2A(B, index2);
    end
end

A = size(CO2A);
RowNumber = A(1); 

for index3 = 1:54
    ACIArray(index3,Total) = CO2A(RowNumber, index3);
end

%  Here is to generate the graph
% Here is the place to deterimine whether it CO2 response, O2 response, Light response, or enzyme change. 

Response = 'light';       % The other options are: O2, light, enz;
% 1 is CO2, 3 is light

for n =1: 6
    newplot;
	for m=1:6
        subplot(2,3,m);
        % plot(ACIArray(4,:),ACIArray(m+ 6 + 6 * (n-1),:),'ro');
        plot(ACIArray(3,:),ACIArray(m+ 6 + 6 * (n-1),:),'ro');
        suc=GenOut_Addtitle(Response, m + (n-1)*6);
        if m ==1
            ylabel('umol CO2 m-2 s-1')
        end
	end
    pause;
end

for m=43:54
    subplot(3,4,m-42);
    plot(ACIArray(3,:),ACIArray(m,:),'ro');
    sub=GenOut_AddTitle(Response, m-6);
end

done = ACIArray';
