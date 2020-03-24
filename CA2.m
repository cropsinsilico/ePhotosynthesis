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



% CA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program calls a PS_PRDrive routine several times.  A graph is 
% constructed of the steady-state values of "The net rate of carbon 
% assimilation rate" for different enzyme concentration values.  If the 
% program finishes without doing anything, check the VmaxVary matrix in 
% CA_User_Input.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


CA_User_Input 
global NewCalculate;


global CarbonRate;                      %   From CA_PSPR_Out
global CarbonRate2;                     %   From CA_trDynaPS_Drive
global Metabolites;                     %   From CA_PSPR_Out
global StepSize;                        %   From CA_User_Input
global VmaxMatrix;                      %   From CA_User_Input
global VmaxXValue;
global VmaxVary;                        %   From CA_User_Input
global MetConCoeff;                     %   From CA_User_Input - metabolite to be measured against
global DriveSelector;                   %   From CA_User_Input
global J;                               %   Counts to determine if you should use the 1st or
                                        %    2nd column of VmaxXValue
                                       
global PS_VEL;
global InExcel;
global UseInputVmax;
                                        
VmaxLabels = ['V1      ';
              'V2      ';
              'V3      ';
              'V4      ';  
              'V5      ';
              'V6      ';
              'V7      ';
              'V8      ';
              'V9      ';
              'V10     ';
              'V11     ';
              'V12     ';
              'V13     ';
              'V16     ';
              'V21     ';
              'V22     ';
              'V23     ';
              'V31     ';
              'V32     ';
              'V33     ';
              'V111    ';
              'V112    ';
              'V113    ';
              'V121    ';
              'V122    ';
              'V123    ';
              'V124    ';
              'V131    ';
              'V1T     ';
              'V2T     ';
              'FIBF_PQT';
              'ChlT    ';
              'Chl2T   ';
              'V2M     ';
              'Vmax11  ';
              'Tcyt    ';
              'BF_An   ';
              'BF_Cn   ';
              '        ';
              '        ';
			'	KM11   '	;
			'	KM12   '	;
			'	KM13 	 '	;
			'	KI11 	 '	;
			'	KI12 	 '	;
			'	KI13   '	;
			'	KI14   '	;
			'	KI15   '	;
			'	KM21  	'	;
			'	KM22  	'	;
			'	KM23  	'	;
			'	KM31a 	'	;
			'	KM32b 	'	;
			'	KM41   '	;
			'	KM42  	'	;
			'	KE4    '	    ;
			'	KM51   '	;
			'	KM52   '	;
			'	KM53   '	;
			'	KE5    '	;
			'	KM61   '	;
			'	KI61   '	;
			'	KI62   '	;
			'	KE6    '	;
			'	KM71   '	;
			'	KM72   '	;
			'	KM73   '	;
			'	KM74   '	;
			'	KE7    '	;
			'	KM8	   '	;
			'	KM81  	'	;
			'	KM82  	'	;
			'	KE8    '	;
			'	KM9	   '	;
			'	KI9	   '	;
			'	KE9	   '	;
			'	KM10   '	;
			'	KM101	 '	;
			'	KM102	 '	;
			'	KM103	 '	;
			'	KE10	  '	;
			'	KE11	  '	;
			'	KE12	  '	;
			'KM131  	'	;
			'	KM132 	'	;
			'	KI131	 '	;
			'	KI132 	'	;
            
            
			'	KI133 	'	;
			'	KI134 	'	;
			'	KI135 	'	;
			'	KE13  	'	;
			'	KM161 	'	;
			'	KM162 	'	;
			'	KM163 	'	;
			'	KM16  	'	;
			'	KE21  	'	;
			'	KE22  	'	;
			'	KM231 	'	;
			'	KM232 	'	;
			'	KA231 	'	;
			'	KA232 	'	;
			'	KA233 	'	;
			'	KI23  	'	;
			'	KM311 	'	;
			'	KM312 	'	;
			'	KM313 	'	;
			'	KM32  	'	;
			'	KM33  	'	;		
              '        '    ;
              '        ';
              '        ';
              '        ';
              '        ';
              '        ';
              '        ';
              '        ';
              '        ';
              '        ';
              '        ';
              '        ';
            
            
            
			'	KO	    '	;
			'	KC	    '	;
			'	KR	    '	;
			'	KM112	 '	;
			'	KI1122	'	;
			'	KI1121	'	;
			'	KM1131	'	;
			'	KM1132	'	;
			'	KI113 	'	;
			'	KE113 	'	;
			'	KM121 	'	;
			'	KM1221	'	;
			'	KM1222	'	;
			'	KI1221	'	;
			'	KE122 	'	;
			'	KM123 	'	;
			'	KI123 	'	;
			'	KE123 	'	;
			'	KM1241	'	;
			'	KM1242	'	;
			'	KI124 	'	;
			'	KE124 	'	;
			'	KM1311	'	;
			'	KI1311	'	;
			'	KI1312	'	;
			'	KM1312	'	;
			'	KM1011	'	;
			'	KI1011	'	;
			'	KM1012	'	;
			'	KI1012	'	;  
            
             '		      '	
 '		      '	
 '		      '	
 '		      '	
 '		      '	
 '		      '	
 '		      '	
 '		      '	
 '		      '	
 '		      '	
 '	KE501	 '	;
 '	Km511	 '	;
 '	Km512	 '	;
 '	Km513	 '	;
 '	KE51	  '	;
 '	Km514	 '	;
 '	Km521	 '	;
 '	KI521	 '	;
 '	KI522	 '	;
 '	KI523	 '	;
 '	KE52	  '	;
 '	KE531	 '	;
 '	KE541	 '	;
 '	Km551	 '	;
 '	Km552	 '	;
 '	Km553	 '	;
 '	Km554	 '	;
 '	KE55 	 '	;
 '	Km561	 '	;
 '	Km562	 '	;
 '	KI561	 '	;
 '	KI562	 '	;
 '	KI563	 '	;
 '	KI564	 '	;
 '	KI565	 '	;
 '	KE56 	 '	;
 '	Km571	 '	;
 '	Ki572	 '	;
 '	KE57 	 '	;
 '	Km581	 '	;
 '	KI581	 '	;
 '	KI582	 '	;
 '	Km591	 '	;
 '	Km592	 '	;
 '	Km593	 '	;
 '	KI591	 '	;
 '	KI592	 '	;
 '	KE59 	 '	;
 '	Km601	 '	;
 '	Km602	 '	;
 '	Km603	 '	;
 '	Km604	 '	;
 '	KE60 	 '	;
 '	KE61 	 '	;
 '	Km621	 '	;
 '	V51  	 '	;
 '	V52  	 '	;
 '	V55  	 '	;
 '	V56  	 '	;
 '	V57  	 '	;
 '	V58  	 '	;
 '	V59  	 '	;
 '	V60  	 '	;
 '	V61	   '	;
 '	V62	   '	;
 'Vdhap_in'	;
 '	Vgap_in'	;
 '	Vpga_in'	;

];
          
MetaboliteLabels = ['RuBP        ';'PGA         ';'DPGA        ';  %   Labels for metabolite used to calculate
                    'T3P         ';'FBP         ';'E4P         ';  %    control coefficient.
                    'S7P         ';'SBP         ';'ATP         ';
                    'NADHP       ';'CO2         ';'O2          ';
                    'HexP        ';'PenP        ';'GCEA        ';
                    'GCA         ';'PGCA        ';'GCAc        ';
                    'GOAc        ';'SERc        ';'GLYc        ';
                    'HPRc        ';'GCEAc       ';'Carbon Rate '];

MetLength  = length(MetConCoeff);       % The different metabolites
VmaxLength = length(VmaxVary);          % The different Vmax used
                    
% CA_User_Input;

OriginalValues  = VmaxMatrix(:,1);                                        
VmaxXValue      = zeros(length(VmaxVary),2);   
VmaxYValue                  = zeros(1,2);
ControlCoefficientsVector   = zeros(VmaxLength,1);
ControlCoefficientsTable    = zeros(MetLength,VmaxLength); 
MetConCoeff = CA_DriveSelection(MetConCoeff,DriveSelector);

for f=1:MetLength
    if MetConCoeff(f,1) == 1
        
        for inum=1:VmaxLength  
            VmaxXValue = zeros(length(VmaxMatrix),2);
            VmaxXValue(:,1) = VmaxMatrix;
            VmaxXValue(:,2) = VmaxMatrix;
            if  VmaxVary(inum,1) == 1
                    Vstep = VmaxMatrix(inum,1) * StepSize;
                    VmaxXValue(inum,1) = VmaxMatrix(inum,1) ;
                    VmaxXValue(inum,2) = VmaxMatrix(inum,1) - Vstep ;      
                for Jn=1:2
                    global CA_Vmax;
                    CA_Vmax = zeros(2,1);
                    CA_Vmax(1) = VmaxXValue(inum,Jn);          
                    CA_Vmax(2) = inum;                       
                    if NewCalculate ==1
                            global PSPR_SUCS_com;
                            if DriveSelector == 1                   
                                PSPR_SUCS_com = 0;                  
                                CA_PSDrive;                                             
                            elseif DriveSelector == 2
                                PSPR_SUCS_com = 0;
                                CA_PS_PRDrive;                    
                            elseif DriveSelector == 3
                                CA_trDynaPS_Drive;
                            elseif DriveSelector == 4
                                CA_CMDrive
                            end
                            
                    end
                    
                    if Jn ==2
                           global PSPR_SUCS_com;
                            if DriveSelector == 1                   
                                PSPR_SUCS_com = 0;                  
                                CA_PSDrive;                                             
                            elseif DriveSelector == 2
                                PSPR_SUCS_com = 0;
                                CA_PS_PRDrive;                    
                            elseif DriveSelector == 3
                                CA_trDynaPS_Drive;
                            elseif DriveSelector == 4
                                CA_CMDrive
                            end
                    end
                    
                    
                    
                 % Obtain the metabolite concentration or flux value
                    if f < MetLength
                        VmaxYValue(1,Jn) = Metabolites(length(Metabolites(:,f)),f);
                    else
                        if DriveSelector == 3
                            VmaxYValue(1,Jn) = CarbonRate2(length(CarbonRate2),1);
                        else
                            Jn
                            VmaxYValue(1,Jn) = CarbonRate(1,length(CarbonRate)); 
                            
                        end
                        
                        if Jn ==1
                            if NewCalculate ==1
                                VmaxYControl = VmaxYValue(1,1);
                                NewCalculate =0;
                            end
                        end
                    end
                    
                end
                
                Ychange = (VmaxYValue(1,2)-VmaxYControl);
                YPercentChange = Ychange/VmaxYControl;
                ControlCoeff = YPercentChange/StepSize;
                ControlCoefficientsTable(f,inum) = ControlCoeff;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                
            end
        end       
    end
end
CA_Table(MetLength,VmaxLength,VmaxVary,VmaxLabels,MetConCoeff,MetaboliteLabels,OriginalValues,ControlCoefficientsTable,DriveSelector)
