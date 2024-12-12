
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


% Function [light] = condition; This function is used to store all the required 
% environmental variables, such as light, CO2, O2, humidity as such. This function 
% contains two parts. Part a includes the generic (default) conditions and the 
% second part contains the detailed conditions for different time period.

function fini = Condition (t)
% Allow early abort if environmental inputs handled via ePhotosynthesis
% function
global dontCondition;
if (dontCondition)
    global NumInter_draw;
    NumInter_draw = 10;
    fini = 1;
    return;
end
global TestCa;
global TestLi;
global RUBISCOMETHOD;         % The method for calculation of Rubisco catalyzed reaction
RUBISCOMETHOD = 2;          % 1: Use enzyme concentration for calculation
                            % 2: Use the michaelis menton and enzyme together for calculation
global VolRatioStCyto
VolRatioStCyto =1; 

                            
% First get the generic conditions
                            
global CO2_cond;
global O2_cond;
global GLight;
global V16;
global Temp_cond;

global Cond_V16;        % This variable is transfered from PSInitial for modificatin of V16, the rate of ATP synthesis. 
CO2Temp = TestCa*0.7;%280;          % CO2 concentation  % ppm
O2Temp = 0.21;          % O2 concentration  %default is 0.21, i.e. 21%. 

CO2_cond = CO2Temp /(3 * 10^4);
O2_cond = O2Temp*1.26;
Temp_cond = 25;

light =  TestLi*0.85*0.85;  % light umol m-2 s-1      

% Here the time dependent variable is regulated. 
global tglobal;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%   Here define how many interval needed for the experiments  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NumberInterval = 10;

global NumInter_draw; 
NumInter_draw = NumberInterval;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is an experimental protocol for doing repeatative experiment

Tinter = tglobal/NumberInterval;

FirstMet = 0;

%     for index = 1:NumberInterval 
%             b = index * Tinter;
%         if t <= b & FirstMet == 0
%             
%             modifier = 1 * index; 
%             
%             global BF2XanCycle_pHl;
% 
%             global Xan2Stom_ABA;
% 
%                 
%             % Light regulation
%             light =  100 * index; 
%             light = 2000; 
%             
%             % CO2 regulation
%             temp = 280;
%             %temp = 1000-100 * (index-1); 
%             
%             % O2 regualtion
%             O2Temp = 0.21 ;   
%             
%             CO2_cond = temp /(3 * 10 ^ 4);   
%             O2_cond = O2Temp*1.26;
%             
%             % Regulation of Vmax
%             % Global RuACT_RC;
%             % RuACT_RC(9) = 25 + 25  * index/2;
%             
%             FirstMet = 1;                
%             
%         end
%     end
        
        

% % % PAM MEASUREMENT
% if t < 20
%     light = 7;
% else
%     light=500;
% end
% % % PAM MEASUREMENT
%   StepL = 20; 
%     if t<1
%         light = 8000;
%     elseif t>1 & t < 20
%         light = 1;
%     else 
%         index = floor(t/StepL);
%         if t > (StepL*index) & t < (StepL*index+1)
%             light = 8000; 
%         else 
%             light = 500; 
%         end
%     end   


% if t<200
%     light = 1000;
% elseif t>200 & t<400
%     light = 100;
% else
%     light = 1000;
% end

GLight   = light;
fini = 1;   
