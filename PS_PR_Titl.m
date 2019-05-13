%   Copyright   Xin-Guang Zhu and Stephen P. Long, University of Illinois 
%   Copyright ©  2007

%   This file is part of CarbonMetabolism.

%    CarbonMetabolism is free software; you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation; either version 3 of the License, or
%    (at your option) any later version.

%    CarbonMetabolism is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.

%    You should have received a copy of the GNU General Public License (GPL)
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.

% PS_PR_Titl.m
% This is the routine that give titles for each graph. 

function suc = PS_PR_Titl(m,p,n)

if m == 1
    title('RuBP');
elseif m==2
    title('PGA');
elseif m==3
    title('DPGA');
elseif m==4
    title('T3P');
elseif m==5
    title('FBP');
elseif m==6
    title('E4P');
elseif m==7
    title('S7P');
elseif m==8
    title('SBP');
elseif m==9
    title('ATP');
elseif m==10
    title('NADPH');
elseif m==11
    title('CO2');
elseif m==12
    title('O2');
elseif m==13
    title('HexP');
elseif m==14
    title('PenP');
elseif m==15
    title('GCEA');
elseif m==16
    title('GCA');
elseif m==17
    title('PGCA');
elseif m==18
    title('GCA_c');
elseif m==19
    title('GOA_c');
elseif m==20
    title('SER_c');
elseif m==21
    title('GLY_c');
elseif m==22
    title('HPR_c');
elseif m==23
    title('GCEA_c');
end
    xlabel('time(s)');
    ylabel('Concentration(mM)');
    suc = 1;
