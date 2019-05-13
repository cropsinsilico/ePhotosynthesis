% This function is used to store all the required environmental variables, such as light, CO2, O2, humidity as such. 
%   Copyright ? 2007

%   Copyright   Xin-Guang Zhu and Stephen P. Long, University of Illinois 

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



function fini = cdn (t)

global RUBISCOTOTAL;        % Total concentration of Rubisco.
global RUBISCOMETHOD;       % The method for calculation of Rubisco catalyzed reaction
RUBISCOTOTAL = 3;           
RUBISCOMETHOD = 2;          % 1: Use enzyme concentration for calculation
                            % 2: Use the michaelis menton and enzyme concentration together for calculation
global VolRatioStCyto;
VolRatioStCyto = 1;

global GP;
global gp2condition_RuCon;

if GP ==1
    RUBISCOTOTAL = gp2condition_RuCon;
end


% First get the generic environmental conditions
                            
%global CO2_cond;
global O2_cond;
global GLight;
% global V16;


%CO2Temp = 550*0.7;%375*0.7;%270;  % default is 270 ppm, which corresponds to an atmospheric CO2 concentration of 360 ppm. 
O2Temp = 0.21;  % default is 0.21, i.e. 21%. 

%CO2_cond = CO2Temp /(3 * 10^4);
O2_cond = O2Temp *1.26;

light = 1000; 
GLight = light;

% Here the time dependent variable is regulated. 
% Tinter = 1000000;         % The number of intervals in the simulations
% FirstMet = 0;
% global tglobal;

% TotalIndex = tglobal/Tinter;

%     for index = 1:TotalIndex 
%             b = index * Tinter;
%         if t <= b & FirstMet == 0
%             
%             % Light regulation
%             light = 1000;
%             
%             % CO2 regulation
%             temp =  70 * index;
%             temp = 270; 
%             CO2_cond = temp/(3 * 10 ^ 4);
% 
%             % O2 regualtion
%             O2Temp = 0.21;  % default is 0.21, i.e. 21%. 
%             O2_cond = O2Temp*1.26;
%             
%             FirstMet = 1;                   
%         end
%     end


% To simulate a perturbation experiment, set PLimited to be 1. 

% PLimited = 0; 
% if PLimited ==1
%   if t < 2000
%       O2 = 0.21 * 1.26; 
%   elseif t>2000 & t<4000
%       O2 = 0.021 * 1.26; 
%   else
%       O2 = 0.21 * 1.26; 
%   end
% 
%   O2_cond = O2; 
% 
% % The enzyme activities can also be modified here. For example, to simulate a lower triose phosphate export capacity, we can use the following code. 
% %  global V31; 
% %  global V32; 
% %  global V33; 
% %  V31 = 0.01; 
% %  V32 = 0.01; 
% %  V33 = 0.01; 
%  
% end     % End the P limited testing case. 

% Users can define other kinds of expeirmental scenario in this function. 

fini = 1;
