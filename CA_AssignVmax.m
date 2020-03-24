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


function fin = CA_AssignVmax

global CA_Vmax;

Vmax = CA_Vmax(1);
index = CA_Vmax(2);

    
	switch index
        case 1
            global V1;
            V1 = Vmax;
        case 2
            global V2;
            V2 = Vmax;
        case 3
            global V3;
            V3 = Vmax;
        case 4
            global V4;
            V4 = Vmax;
        case 5
            global V5;
            V5 = Vmax;
        case 6
            global V6;
            V6 = Vmax;
        case 7
            global V7;
            V7 = Vmax;
        case 8
            global V8;
            V8 = Vmax;
        case 9
            global V9;
            V9 = Vmax;
        case 10
            global V10;
            V10 = Vmax;
        case 11
            global V11;
            V11 = Vmax;
        case 12
            global V12;
            V12 = Vmax;
        case 13 
            global V13;
            V13 = Vmax;
        case 14
            global V16;
            V16 = Vmax;
        case 15
            global V21;
            V21 = Vmax; 
        case 16
            global V22;
            V22 = Vmax;
        case 17
            global V23;
            V23 = Vmax;
        case 18 
            global V31;
            V31 = Vmax;
        case 19
            global V32;
            V32 = Vmax;
        case 20
            global V33;
            V33 = Vmax;
        case 21
            global V111;
            V111 = Vmax;
        case 22
            global V112;
            V112 = Vmax;
        case 23
            global V113;
            V113 = Vmax;
        case 24
            global V121;
            V121 = Vmax;
        case 25
            global V122;
            V122 = Vmax;
        case 26
            global V123;
            V123 = Vmax;
        case 27
            global V124;
            V124 = Vmax;
        case 28
            global V131;
            V131 = Vmax;
        case 29 
            global V1T;
            V1T = Vmax;
        case 30
            global V2T;
            V2T = Vmax;
        case 31
            global FIBF_PQT;
            FIBF_PQT = Vmax;
        case 32
            global ChlT;
            ChlT = Vmax;
        case 33
            global ChlT2
            ChlT2 = Vmax;
        case 34
            global V2M;
            V2M = Vmax;
        case 35
            global Vmax11;
            Vmax11 = Vmax;
        case 36
            global Tcyt;
            Tcyt = Vmax;
        case 37
            global BF_An;
            BF_An = Vmax;
        case 38
            global BF_Cn;
            BF_Cn = Vmax;
    
    % First get set the variable to be global
case 41    
global	KM11	;
KM11 = Vmax;

case 42
global	KM12	;
KM12= Vmax;

case 43
global	KM13	;
KM13= Vmax;

case 44
global	KI11	;
KI11= Vmax;

case 45
global	KI12	;
KI12= Vmax;

case 46
global	KI13	;
KI13= Vmax;

case 47
global	KI14	;
KI14= Vmax;

case 48
global KI15	;
KI15= Vmax;

case 49
global	KM21	;
KM21= Vmax;

case 50
global	KM22	;
KM22= Vmax;

case 51
global	KM23	;
KM23= Vmax;

case 52
global	KM31a	;
KM31a= Vmax;

case 53
    global	KM32b	;
    KM32b= Vmax;
    
    
case 54
global	KM41	;
KM41= Vmax;

case 55
    global	KM42	;
    KM42= Vmax;
    
case 56
global	KE4	;
KE4= Vmax;

case 57

  global	KM51	;
  KM51= Vmax;
  
case 58
global	KM52	;
KM52= Vmax;

case 59
    global	KM53	;
    KM53= Vmax;
    
case 60
global	KE5	;
KE5= Vmax;

case 61
global	KM61	;
KM61= Vmax;

case 62
global	KI61	;
KI61= Vmax;

case 63
global	KI62	;
KI62= Vmax;


case 64
global	KE6	;
KE6= Vmax;


case 65
global	KM71	;
KM71= Vmax;


case 66
global	KM72	;
KM72= Vmax;


case 67
global	KM73	;
KM73= Vmax;


case 68 
global	KM74	;
KM74= Vmax;


case 69
global	KE7	;
KE7= Vmax;


case 70
global	KM8	;
KM8= Vmax;


case 71
global	KM81	;
KM81= Vmax;


case 72
global	KM82	;
KM82= Vmax;


case 73
global	KE8	;
KE8= Vmax;


case 74
global	KM9	;
KM9= Vmax;


case 75
global	KI9	;
KI9= Vmax;


case 76
global	KE9	;
KE9= Vmax;


case 77
global	KM10	;
KM10= Vmax;


case 78
global	KM101	;
KM101= Vmax;


case 79
global	KM102	;
KM102= Vmax;


case 80
global	KM103	;
KM103= Vmax;

case 81
global	KE10	;
KE10= Vmax;


case 82
global	KE11	;
KE11= Vmax;


case 83
global	KE12	;
KE12= Vmax;


case 84
global	KM131	;
KM131= Vmax;


case 85
global	KM132	;
KM132= Vmax;


case 86
global	KI131	;
KI131= Vmax;


case 87
global	KI132	;
KI132= Vmax;


case 88
global	KI133	;
KI133= Vmax;


case 89
global	KI134	;
KI134= Vmax;


case 90
global	KI135	;
KI135= Vmax;


case 91
global	KE13	;
KE13= Vmax;


case 92
global	KM161	;
KM161= Vmax;


case 93
global	KM162	;
KM162= Vmax;


case 94
global	KM163	;
KM163= Vmax;


case 95
global	KM16	;
KM16= Vmax;


case 96
global	KE21	;
KE21= Vmax;


case 97
global	KE22	;
KE22= Vmax;


case 98
global	KM231	;
KM231= Vmax;


case 99
global	KM232	;
KM232= Vmax;


case 100
global	KA231	;
KA231= Vmax;


case 101
global	KA232	;
KA232= Vmax;


case 102
global	KA233	;
KA233= Vmax;


case 103
global	KI23	;
KI23= Vmax;


case 104
global	KM311	;
KM311= Vmax;


case 105
global	KM312	;
KM312= Vmax;


case 106
global	KM313	;
KM313= Vmax;


case 107
global	KM32	;
KM32= Vmax;


case 108
global	KM33	;
KM33= Vmax;


case 121
global	KO	;
KO= Vmax;


case 122
global	KC	;
KC= Vmax;


case 123
global	KR	;
KR= Vmax;


case 124
global	KM112	;
KM112= Vmax;


case 125
global	KI1122	;
KI1122= Vmax;


case 126
global	KI1121	;
KI1121= Vmax;


case 127
global	KM1131	;
KM1131= Vmax;


case 128
global	KM1132	;
KM1132= Vmax;


case 129
global	KI113	;
KI113= Vmax;


case 130
global	KE113	;
KE113= Vmax;


case 131
global	KM121	;
KM121= Vmax;


case 132
global	KM1221	;
KM1221= Vmax;


case 133
global	KM1222	;
KM1222= Vmax;


case 134
global	KI1221	;
KI1221= Vmax;


case 135
global	KE122	;
KE122= Vmax;


case 136
global	KM123	;
KM123= Vmax;


case 137
global	KI123	;
KI123= Vmax;


case 138
global	KE123	;
KE123= Vmax;


case 139
global	KM1241	;
KM1241= Vmax;


case 140
global	KM1242	;
KM1242= Vmax;


case 141
global	KI124	;
KI124= Vmax;


case 142
global	KE124	;
KE124= Vmax;


case 143
global	KM1311	;
KM1311= Vmax;


case 144
global	KI1311	;
KI1311= Vmax;


case 145
global	KI1312	;
KI1312= Vmax;


case 146
global	KM1312	;
KM1312= Vmax;


case 147
global	KM1011	;
KM1011= Vmax;


case 148
global	KI1011	;
KI1011= Vmax;


case 149
global	KM1012	;
KM1012= Vmax;


case 150
global	KI1012	;
KI1012= Vmax;


case 161
global KE501
KE501= Vmax;

case 162
global Km511
Km511= Vmax;

case 163
global Km512
Km512= Vmax;


case 164
global Km513
Km513= Vmax;


case 165
global KE51
KE51= Vmax;


case 166
global Km514
Km514= Vmax;


case 167
global Km521
Km521= Vmax;


case 168
global KI521
KI521= Vmax;


case 169
global KI522
KI522= Vmax;


case 170
global KI523
KI523= Vmax;


case 171
global KE52
KE52= Vmax;


case 172
global KE531
KE531= Vmax;


case 173
global KE541
KE541= Vmax;


case 174
global Km551
Km551= Vmax;


case 175
global Km552
Km552= Vmax;


case 176
global Km553
Km553= Vmax;


case 177
global Km554
Km554= Vmax;


case 178
global KE55
KE55= Vmax;


case 179
global Km561
Km561= Vmax;


case 180
global Km562
Km562= Vmax;


case 181
global KI561
KI561= Vmax;


case 182
global KI562
KI562= Vmax;


case 183
global KI563
KI563= Vmax;


case 184
global KI564
KI564= Vmax;

case 185
global KI565
KI565= Vmax;


case 186
global KE56
KE56= Vmax;

case 187
global Km571
Km571= Vmax;

case 188
global Ki572
Ki572= Vmax;


case 189
global KE57
KE57= Vmax;


case 190
global Km581
Km581= Vmax;


case 191
global KI581
KI581= Vmax;


case 192
global KI582
KI582= Vmax;


case 193
global Km591
Km591= Vmax;


case 194
global Km592
Km592= Vmax;


case 195
global Km593
Km593= Vmax;


case 196
global KI591
KI591= Vmax;


case 197
global KI592
KI592= Vmax;


case 198
global KE59
KE59= Vmax;


case 199
global Km601
Km601= Vmax;


case 200
global Km602
Km602= Vmax;


case 201
global Km603
Km603= Vmax;


case 202
global Km604
Km604= Vmax;


case 203
global KE60
KE60= Vmax;


case 204
global KE61
KE61= Vmax;


case 205
global Km621
Km621= Vmax;


case 206
global V51
V51= Vmax;


case 207
global V52
V52= Vmax;


case 208
global V55
V55= Vmax;


case 209
global V56
V56= Vmax;


case 210
global V57
V57= Vmax;


case 211
global V58
V58= Vmax;


case 212
global V59
V59= Vmax;


case 213
global V60
V60= Vmax;


case 214
global V61
V61= Vmax;


case 215
global V62
V62= Vmax;


case 216
global Vdhap_in
Vdhap_in= Vmax;


case 217
global Vgap_in 
Vgap_in= Vmax;


case 218
global Vpga_in  
Vpga_in= Vmax;


end
