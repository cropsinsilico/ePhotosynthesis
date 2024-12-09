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
function Resulta=DynaPS_Drive(ParaNum, Ratio)
% trDynaPS_Drive.m
% This part include the function to begin the simulation.

% The time information is set in a global variable called tglobal in SYSInitial. 
InitRatios(ParaNum, Ratio);

% DynaPS_Drive.m
% This part include the function to begin the simulation.

% The time information is set in a global variable called tglobal in SYSInitial. 
Begin = 1;
fin = SYSInitial(Begin);
global options1;
global tglobal;
time = tglobal;

%%%%%%%%%%%%%%%%%%%%%%%%
%   Calculation  step %
%%%%%%%%%%%%%%%%%%%%%%%%

global ATPActive;
ATPActive = 0;

global EPS_ATP_Rate;        % Indicate in the beginning there is no ATP synthesis activity.
EPS_ATP_Rate = 0;

ModelComb = IniModelCom;        % Initialize the structure of the model, i.e. Is this model separate or combined with others. 

global BF_FI_com;            % The combination of BF and FI model 
BF_FI_com = 1;

global PR_PS_com;    % This is a variable indicating whether the PR model is actually need to be combined with PS or not. If 1 then means combined; 0 means not. 
PR_PS_com = 1;

global FIBF_PSPR_com; % 1 means that the overall EPS model is used. 0 means partial model of FIBF is used. 
FIBF_PSPR_com = 1;    

global RuACT_EPS_com;     % A global variable to indicate whether the RuACT is run by itself or combined with others. 
RuACT_EPS_com = 1;        % Since this is run within this program, it is combinbed, therefore, it is assigned value 1, otherwise, assign value 0. 

global RedoxReg_RA_com;     % This is the connection between Redox and RA.
RedoxReg_RA_com = 0;        % This means that the connection is not provided there. 

global XanCycle_BF_com;
XanCycle_BF_com = 1;

global EPS_SUCS_com;
EPS_SUCS_com = 1;

global PSPR_SUCS_com;    % This is a variable indicating whether the PSPR model is actually need to be combined with SUCS or not. If 1 then means combined; 0 means not. 
PSPR_SUCS_com = 1;

SUCS_Param = 0;

global CO2A;
CO2A = zeros(5,1);


% Next is to initialize the vector. 

DynaPS_Con = DynaPS_Ini;

va1 = 0;
global PS12ratio; % The ratio of the PSI unit to the PSII unit
BF_Param = zeros(5,1);
BF_Param(1) = va1;
BF_Param(2) = PS12ratio;

FI_Param = zeros(5,1);
FI_Param(1) = va1;
FI_Param(2) = PS12ratio;

PS_PR_Param = 0;
EPS_Param = 0;

RuACT_Param = zeros(2,1);
RuACT_Param(1) = va1;
RuACT_Param(2) = PS12ratio;

XanCycle_Param = zeros(2,1);
XanCycle_Param(1) = va1;
XanCycle_Param(2) = PS12ratio;

RedoxReg_Param = 0; % This parameter is just used here as a future storage tool. NOt used now. 

[Tt,d] = ode15s(@DynaPS_mb,[0,time],DynaPS_Con,options1,BF_Param, FI_Param, PS_PR_Param, RuACT_Param, RedoxReg_Param, XanCycle_Param, SUCS_Param);


%done = DynaPS_Graph(Tt,d);
global BF_VEL;
global FI_VEL;
global BF_CON; 
global PS_VEL;
global PR_VEL; 
global FI_CON; 
global PS_CON; 
global PR_CON;
global SUCS_VEL;
global RuACT_VEL; 
global XanCycle_VEL;
global RedoxReg_VEL;
global RROEA_VEL;
global AVR; 
[row,col]=size(RuACT_VEL);
    PSIIabs=FI_VEL( :,  57 ) ;
    PSIabs=BF_VEL(:, 11	);
    %PSIabs2=BF_VEL(:, 14)+BF_VEL(:, 16);
    temp = RuACT_VEL(:,6) ;
    CarbonRate = temp * AVR; 
    VPR=RuACT_VEL(:,7)* AVR;
%     CO2Release = PR_VEL(:,9) * AVR; 
%     Assim = CarbonRate - CO2Release;
    Vpgasink=SUCS_VEL(:,15)'*AVR;
    VStarch=(PS_VEL(14,:)-PS_VEL(20,:))*AVR;
    Vsucrose=SUCS_VEL(:,11)'*AVR;
    Resulta=[0;0;0;0;0;0;0];
    Resulta(1)=PSIIabs(row);
    Resulta(2)=PSIabs(row);
    %Resulta(3)=PSIabs2(row);
    Resulta(3)=CarbonRate(row);
    Resulta(4)=VPR(row);
    Resulta(5)=Vpgasink(row);
    Resulta(6)=Vsucrose(row);
    Resulta(7)=VStarch(row);
    global FluxTR;
    FluxTR=zeros(81,1);
    FluxTR(1)=RuACT_VEL(row,6);%PS
    FluxTR(2)=RuACT_VEL(row,7);%PR
    FluxTR(3)=PS_VEL(3,row);% v2
    FluxTR(4)=PS_VEL(4,row);% v3
    FluxTR(5)=PS_VEL(6,row);% v5
    FluxTR(6)=PS_VEL(7,row);% v6
    FluxTR(7)=PS_VEL(8,row);% v7
    FluxTR(8)=PS_VEL(9,row);% v8
    FluxTR(9)=PS_VEL(10,row);% v9
    FluxTR(10)=PS_VEL(11,row);% v10
    FluxTR(11)=PS_VEL(12,row);% v13
    FluxTR(12)=PS_VEL(14,row);% v23
    FluxTR(13)=PS_VEL(19,row);% v24
    FluxTR(14)=PS_VEL(20,row);% v25
    FluxTR(15)=PR_VEL(row,3);%v112
    FluxTR(16)=PR_VEL(row,4);%v113
    FluxTR(17)=PR_VEL(row,5);%v121
    FluxTR(18)=PR_VEL(row,6);%v122
    FluxTR(19)=PR_VEL(row,7);%v123
    FluxTR(20)=PR_VEL(row,8);%v124
    FluxTR(21)=PR_VEL(row,9);%v131
    FluxTR(22)=PR_VEL(row,10);%vlin
    FluxTR(23)=PR_VEL(row,11);%v2out
    FluxTR(24)=SUCS_VEL(row,2);%v51	;%	DHAP+GAP --FBP
    FluxTR(25)=SUCS_VEL(row,3);%v52	;%	FBP --F6P + Pi
    FluxTR(26)=SUCS_VEL(row,4);%v55	;%	G1P+UTP --OPOP+UDPG 
    FluxTR(27)=SUCS_VEL(row,5);%v56	;%	UDPG+F6P--SUCP + UDP
    FluxTR(28)=SUCS_VEL(row,6);%v57	;%	SUCP--Pi + SUC
    FluxTR(29)=SUCS_VEL(row,7);%v58	;%	F26BP--F6P + Pi
    FluxTR(30)=SUCS_VEL(row,8);%v59	;%	F6P + ATP --ADP + F26BP
    FluxTR(31)=SUCS_VEL(row,9);%v60	;%	ATP+UDP --UTP + ADP
    FluxTR(32)=SUCS_VEL(row,11);%v62	;%	SUC SINK 
    FluxTR(33)=SUCS_VEL(row,12);%vdhap_in	;%	DHAP export from chloroplast
    FluxTR(34)=SUCS_VEL(row,13);%vgap_in	;%	GAP Export from chloroplast
    FluxTR(35)=SUCS_VEL(row,14);%vpga_in	;%	PGA export from chloroplast
    FluxTR(36)=SUCS_VEL(row,15);%vpga_use	;%	PGA utilisation in cytosol
    FluxTR(37:66)=BF_VEL(row,2:31);
    FluxTR(67:124)=FI_VEL(row,2:59);
    FluxTR(125:131)=XanCycle_VEL(row,2:8);
    %FluxTR(132:142)=RROEA_VEL(row,2:12);
    FluxTR(1:36)=FluxTR(1:36)*AVR;
    FluxTR(47)=FluxTR(47)*AVR;
    FluxTR(65)=FluxTR(65)*AVR/2;
% This is to set the regualtions to be as beginning. 
    ATPActive = 0;
    BF_FI_com = 0;
    PR_PS_com = 0;
    FIBF_PSPR_com = 0;    
    RuACT_EPS_com = 0;
    RedoxReg_RA_com = 0;
    XanCycle_BF_com = 0;
% global BF_VEL;
% global FI_VEL;
% global PS_VEL;
clock
    
IniModelCom;
%save FDC2
