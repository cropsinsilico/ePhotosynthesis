%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%   Copyright   Xin-Guang Zhu, Yu Wang, Donald R. ORT and Stephen P. LONG
%CAS-MPG Partner Institute for Computational Biology, Shanghai Institutes for Biological Sciences, CAS, Shanghai,200031 
%China Institute of Genomic Biology and Department of Plant Biology, Shanghai Institutes for Biological Sciences, CAS, Shanghai,200031 
%University of Illinois at Urbana Champaign
%Global Change and Photosynthesis Research Unit, USDA/ARS, 1406 Institute of Genomic Biology, Urbana, IL 61801, USA.
 
%   This file is part of e-photosynthesis
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
function Resulta=trDynaPS_Drive(Cai,Lii,ATPCost,SucPath,ParaNum, Ratio)
global TestCa;
global TestLi;
global TestSucPath;
global TestATPCost;
TestCa=Cai;
TestLi=Lii;
TestSucPath=SucPath;
TestATPCost=ATPCost;
% trDynaPS_Drive.m
% This part include the function to begin the simulation.

% The time information is set in a global variable called tglobal in SYSInitial. 
global PSRatio;
PSRatio=ones(103,1); 
if ParaNum <= 103
PSRatio(ParaNum)=Ratio;
end
global SUCRatio;
SUCRatio=ones(66,1);
if ParaNum > 103&&ParaNum<=169
SUCRatio(ParaNum-103)=Ratio;
end
global PRRatio;
PRRatio=ones(48,1);
if ParaNum > 169&&ParaNum<=217
PRRatio(ParaNum-169)=Ratio;
end

global RacRatio;
RacRatio=ones(16,1);
if ParaNum > 217&&ParaNum<=233
RacRatio(ParaNum-217)=Ratio;
end

global FIRatio;
FIRatio=ones(23,1);
if ParaNum > 233&&ParaNum<=256
FIRatio(ParaNum-233)=Ratio;
end

global BFRatio;
BFRatio=ones(49,1);
if ParaNum > 256&&ParaNum<=305
BFRatio(ParaNum-256)=Ratio;
end    


global XanRatio;
XanRatio=ones(4,1);
if ParaNum > 305&&ParaNum<=309
XanRatio(ParaNum-305)=Ratio;
end 

%clear; 
time1 = clock;
suc = PreProcess;

Begin = 1;
fin = SYSInitial(Begin);
global options1;
global tglobal;
time = tglobal;


global ATPActive;
ATPActive = 0;

global EPS_ATP_Rate;        % Indicate in the beginning there is no ATP synthesis activity.
EPS_ATP_Rate = 0;

ModelComb = IniModelCom;        % Initialize the structure of the model, i.e. Is this model separate or combined with others. 

global BF_FI_com;            % The combination of BF and FI model 
BF_FI_com = 1;

global PR_PS_com;    % This is a variable indicating whether the PR model is actually need to be combined with PS or not. If 1 then means combined; 0 means not. 
PR_PS_Com = 1;

global FIBF_PSPR_com; % 1 means that the overall EPS model is used. 0 means partial model of FIBF is used. 
FIBF_PSPR_com = 1;    

global RuACT_EPS_com;     % A global variable to indicate whether the RuACT is run by itself or combined with others. 
RuACT_EPS_com = 1;        % Since this is run within this program, it is combinbed, therefore, it is assigned value 1, otherwise, assign value 0. 

global RedoxReg_RA_com;     % This is the connection between Redox and RA.
RedoxReg_RA_com = 0;        % This means that the connection is there. 

global XanCycle_BF_com;
XanCycle_BF_com = 1;

global RROEA_EPS_com;
RROEA_EPS_com = 1;

global EPS_SUCS_com;
EPS_SUCS_com = 1;

global PSPR_SUCS_com;    % This is a variable indicating whether the PSPR model is actually need to be combined with SUCS or not. If 1 then means combined; 0 means not. 
PSPR_SUCS_com = 1;

global CO2A;
CO2A = zeros(5,1);

% Next is to initialize the vector. 
global d;
trDynaPS_Con = trDynaPS_Ini;

suc = ParamSet;

[FI_Param, BF_Param, PS_PR_Param, SUCS_Param, EPS_Param, RuACT_Param, XanCycle_Param, RROEA_Param, RedoxReg_Param] = ParamSet;

[Tt,d] = ode15s(@trDynaPS_mb,[0,time],trDynaPS_Con,options1,BF_Param, FI_Param, PS_PR_Param, SUCS_Param, RuACT_Param, RedoxReg_Param, XanCycle_Param, RROEA_Param);
ATPActive = 0;
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
    %Vsucrose=SUCS_VEL(:,11)'*AVR;
    Vt3p=(PS_VEL(15,:)+PS_VEL(17,:))*AVR;
   % Vt3p=(SUCS_VEL(:,2))'*AVR*2;%+SUCS_VEL(:,13)
    Resulta=[0;0;0;0;0;0;0];
%     Resulta(1)=PSIIabs(row);
%     Resulta(2)=PSIabs(row);
    %Resulta(3)=PSIabs2(row);

    Resulta(1)=CarbonRate(row);%Vc
    Resulta(2)=VPR(row);%Vo
    Resulta(3)=Vpgasink(row);%PGA
    Resulta(4)=Vt3p(row);%VT3P
    Resulta(5)=VStarch(row);%Vstarch
    Resulta(6)=PR_VEL(row,10)*AVR;%Vt_glycerate
    Resulta(7)=PR_VEL(row,11)*AVR;%Vt_glycolate
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
    FluxTR(132:142)=RROEA_VEL(row,2:12);
    FluxTR(1:36)=FluxTR(1:36)*AVR;
    FluxTR(47)=FluxTR(47)*AVR;
    FluxTR(65)=FluxTR(65)*AVR/2;
    
    
%FluxTR(132:133)=RedoxReg_VEL(row,2:3);    
%RedoxReg_VEL(2,RedoxReg_TIME_N) = Vred;
%RedoxReg_VEL(3,RedoxReg_TIME_N) = Vox;
% RROEA_VEL	(	RROEA_TIME_N	,   2	)	=	ve2GAPDH	;	
% RROEA_VEL	(	RROEA_TIME_N	,   3	)	=	ve2FBPase	;
% RROEA_VEL	(	RROEA_TIME_N	,   4	)	=	ve2SBPase	;
% RROEA_VEL	(	RROEA_TIME_N	,   5	)	=	ve2PRK	;	
% RROEA_VEL	(	RROEA_TIME_N	,   6	)	=	ve2ATPase	;	
% RROEA_VEL	(	RROEA_TIME_N	,   7	)	=	ve2ATPGPP	;	
% RROEA_VEL	(	RROEA_TIME_N	,   8	)	=	ve2MDH	;
% RROEA_VEL	(	RROEA_TIME_N	,   9	)	=	ve2Fd	;
% RROEA_VEL	(	RROEA_TIME_N	,   10	)	=	veFd2Thio	;
% RROEA_VEL	(	RROEA_TIME_N	,   11	)	=	veFd2Calvin	;
% RROEA_VEL	(	RROEA_TIME_N	,   12	)	=	ve2RuACT	;

% FI_VEL	(	FI_TIME_N	,	2	)	=	vA_d;
% FI_VEL	(	FI_TIME_N	,	3	)	=	vA_f;
% FI_VEL	(	FI_TIME_N	,	4	)	=	vA_U;
% FI_VEL	(	FI_TIME_N	,	5	)	=	vU_A;
% FI_VEL	(	FI_TIME_N	,	6	)	=	vU_f;
% FI_VEL	(	FI_TIME_N	,	7	)	=	vU_d;
% FI_VEL	(	FI_TIME_N	,	8	)	=	v1;
% FI_VEL	(	FI_TIME_N	,	9	)	=	v_r1;
% FI_VEL	(	FI_TIME_N	,	10	)	=	vS1_S2;
% FI_VEL	(	FI_TIME_N	,	11	)	=	vS2_S3;
% FI_VEL	(	FI_TIME_N	,	12	)	=	vS3_S0;
% FI_VEL	(	FI_TIME_N	,	13	)	=	vS0_S1;
% FI_VEL	(	FI_TIME_N	,	14	)	=	vz_1;
% FI_VEL	(	FI_TIME_N	,	15	)	=	v1z_1;
% FI_VEL	(	FI_TIME_N	,	16	)	=	v2z_1;
% FI_VEL	(	FI_TIME_N	,	17	)	=	v3z_1;
% FI_VEL	(	FI_TIME_N	,	18	)	=	v0z_1;
% FI_VEL	(	FI_TIME_N	,	19	)	=	vz_2;
% FI_VEL	(	FI_TIME_N	,	20	)	=	v1z_2;
% FI_VEL	(	FI_TIME_N	,	21	)	=	v2z_2;
% FI_VEL	(	FI_TIME_N	,	22	)	=	v3z_2;
% FI_VEL	(	FI_TIME_N	,	23	)	=	v0z_2;
% FI_VEL	(	FI_TIME_N	,	24	)	=	v1z;
% FI_VEL	(	FI_TIME_N	,	25	)	=	v2z;
% FI_VEL	(	FI_TIME_N	,	26	)	=	v3z;
% FI_VEL	(	FI_TIME_N	,	27	)	=	v0z;
% FI_VEL	(	FI_TIME_N	,	28	)	=	vAB1;
% FI_VEL	(	FI_TIME_N	,	29	)	=	vBA1;
% FI_VEL	(	FI_TIME_N	,	30	)	=	vAB2;
% FI_VEL	(	FI_TIME_N	,	31	)	=	vBA2;
% FI_VEL	(	FI_TIME_N	,	32	)	=	v3;
% FI_VEL	(	FI_TIME_N	,	33	)	=	v_r3;
% FI_VEL	(	FI_TIME_N	,	34	)	=	v3_n;
% FI_VEL	(	FI_TIME_N	,	35	)	=	v_r3_n;
% FI_VEL	(	FI_TIME_N	,	36	)	=	v_pq_ox;
% FI_VEL	(	FI_TIME_N	,	37	)	=	Ic ;
% FI_VEL	(	FI_TIME_N	,	38	)	=	Ia;
% FI_VEL	(	FI_TIME_N	,	39	)	=	v2_1;
% FI_VEL	(	FI_TIME_N	,	40	)	=	v2_2;
% FI_VEL	(	FI_TIME_N	,	41	)	=	v2_00_1;
% FI_VEL	(	FI_TIME_N	,	42	)	=	v2_01_1;
% FI_VEL	(	FI_TIME_N	,	43	)	=	v2_02_1;
% FI_VEL	(	FI_TIME_N	,	44	)	=	v2_00_2;
% FI_VEL	(	FI_TIME_N	,	45	)	=	v2_01_2;
% FI_VEL	(	FI_TIME_N	,	46	)	=	v2_02_2;
% FI_VEL	(	FI_TIME_N	,	47	)	=	vr2_00_1;
% FI_VEL	(	FI_TIME_N	,	48	)	=	vr2_01_1;
% FI_VEL	(	FI_TIME_N	,	49	)	=	vr2_02_1;
% FI_VEL	(	FI_TIME_N	,	50	)	=	vr2_1;
% FI_VEL	(	FI_TIME_N	,	51	)	=	vr2_00_2;
% FI_VEL	(	FI_TIME_N	,	52	)	=	vr2_01_2;
% FI_VEL	(	FI_TIME_N	,	53	)	=	vr2_02_2;
% FI_VEL	(	FI_TIME_N	,	54	)	=	vr2_2;
% FI_VEL	(	FI_TIME_N	,	55	)		=	vP680qU	;	%	vr2_2	
% FI_VEL	(	FI_TIME_N	,	56	)	=	vP680qA	;	%	vr2_2	
% FI_VEL	(	FI_TIME_N	,	57	)	 =   vU_P680;
% FI_VEL	(	FI_TIME_N	,	58	)	  =   vP680_d;
% FI_VEL	(	FI_TIME_N	,	59	)	   =   vP680_f;		
% BF_VEL	(	BF_TIME_N	,	2	)	=	Vbf1	;
% BF_VEL	(	BF_TIME_N	,	3	)	=	Vbf2	;
% BF_VEL	(	BF_TIME_N	,	4	)	=	Vbf3	;
% BF_VEL	(	BF_TIME_N	,	5	)	=	Vbf4	;
% BF_VEL	(	BF_TIME_N	,	6	)	=	Vbf5	;
% BF_VEL	(	BF_TIME_N	,	7	)	=	Vbf6	;
% BF_VEL	(	BF_TIME_N	,	8	)	=	Vbf7	;
% BF_VEL	(	BF_TIME_N	,	9	)	=	Vbf8	;
% BF_VEL	(	BF_TIME_N	,	10	)	=	Vbf9	;
% BF_VEL	(	BF_TIME_N	,	11	)	=	Vbf10	;
% BF_VEL	(	BF_TIME_N	,	12	)	=	Vbf11	;
% BF_VEL	(	BF_TIME_N	,	13	)	=	Vqi	    ;
% BF_VEL	(	BF_TIME_N	,	14	)	=	Vipc	;
% BF_VEL	(	BF_TIME_N	,	15	)	=	Vicp	;
% BF_VEL	(	BF_TIME_N	,	16	)	=	Vinc	;
% BF_VEL	(	BF_TIME_N	,	17	)	=	Vinp	;
% BF_VEL	(	BF_TIME_N	,	18	)	=	Vdp	    ;
% BF_VEL	(	BF_TIME_N	,	19	)	=	Vdc	    ;
% BF_VEL	(	BF_TIME_N	,	20	)	=	Vfp	    ;
% BF_VEL	(	BF_TIME_N	,	21	)	=	Vfc	    ;
% BF_VEL	(	BF_TIME_N	,	22	)	=	Vsfd	;
% BF_VEL	(	BF_TIME_N	,	23	)	=	VsATP	;
% BF_VEL	(	BF_TIME_N	,	24	)	=	VgPQH2	;
% BF_VEL	(	BF_TIME_N	,	25	)	=	Vbf12	;
% BF_VEL	(	BF_TIME_N	,	26	)	=	Vbf13	;
% BF_VEL	(	BF_TIME_N	,	27	)	=	Vbf14	;
% BF_VEL	(	BF_TIME_N	,	28	)	=	Vbf15	;
% BF_VEL	(	BF_TIME_N	,	29	)	=	Vbf16	;
% BF_VEL	(	BF_TIME_N	,	30	)	=	vbfn2	;
% BF_VEL	(	BF_TIME_N	,	31	)	=	VsNADPH	;

% PS_VEL(13,PS_TIME_N) = v16;
% PS_VEL(15,PS_TIME_N) = v31;
% PS_VEL(16,PS_TIME_N) = v32;
% PS_VEL(17,PS_TIME_N) = v33;
% PS_VEL(18,PS_TIME_N) = Pi;

%     TSave = PR_VEL(:,1);
% time2 = clock;
% TotalRunTime = time2-time1
    
IniModelCom;
%save LRC
end