function Arate = ePhotosynthesis(driver, EnvFile, GRNFile, EnzymeFile, ATPCostFile, ...
                                 OutputFile, OutputParam, WeatherTemp, ...
                                 Radiation_PAR, Air_CO2, O2i, ...
                                 ProteinTotalRatio, GRNCin, GRNTin, ...
                                 RUBISCOMETHODin, RUBISCOTOTALin)
arguments
    driver (1,1) string
    EnvFile (1,1) string = ""
    GRNFile (1,1) string = ""
    EnzymeFile (1,1) string = ""
    ATPCostFile (1,1) string = ""
    OutputFile (1,1) string = "output.data"
    OutputParam (1,1) double = 0
    WeatherTemp (1,1) double = 25
    Radiation_PAR (1,1) double = 470*0.85*0.85
    Air_CO2 (1,1) double = 400
    O2i (1,1) double = 0.21  % i.e. 21%
    ProteinTotalRatio (1,1) double = 0.973
    GRNCin (1,1) double = 0
    GRNTin (1,1) double = 0
    RUBISCOMETHODin (1,1) double = 2
    RUBISCOTOTALin (1,1) double = 3
end
driver = str2double(driver);

global GRNC;
global GRNT;
global pcfactor;
global cATPsyn;
global CPSi;
global cNADPHsyn;
global cpsii;
global VfactorC;
global CO2_Env;
global CO2_cond;
global O2_cond;
global LI;
global Tp;
global RUBISCOMETHOD;
global RUBISCOTOTAL;

SucPath = 0;
ATPCost = 0.0;
GRNC = GRNCin;
GRNT = GRNTin;
RUBISCOMETHOD = RUBISCOMETHODin;
RUBISCOTOTAL = RUBISCOTOTALin;

if (strlength(EnvFile) > 0)
    Env_data = ReadParam(EnvFile);
    if isKey(Env_data, "WeatherTemp")
        WeatherTemp = Env_data("WeatherTemp");
    end
    if isKey(Env_data, "Radiation_PAR")
        Radiation_PAR = Env_data("Radiation_PAR");
    end
    if isKey(Env_data, "Air_CO2")
        Air_CO2 = Env_data("Air_CO2");
    end
    if isKey(Env_data, "ProteinTotalRatio")
        ProteinTotalRatio = Env_data("ProteinTotalRatio");
    end
    if isKey(Env_data, "GRNC")
        GRNC = Env_data("GRNC");
    end
    if isKey(Env_data, "GRNT")
        GRNT = Env_data("GRNT");
    end
    if isKey(Env_data, "SucPath")
        SucPath = Env_data("SucPath");
    end
end
if (strlength(ATPCostFile) > 0)
    ATPCost_data = ReadParam(ATPCostFile);
    if isKey(ATPCost_data, "ATPCost")
        ATPCost = ATPCost_data("ATPCost");
    end
end

CO2i = Air_CO2 * 0.7; % intercellular CO2
PPFDi = Radiation_PAR * 1E6/(2.35E5); % Convert W m^{-2} to u moles m^{-2} s^{-1}
CO2_Env = CO2i;
CO2_cond = CO2_Env/(3 * 10^4);
O2_cond = O2i * 1.26;
LI = PPFDi/30;
Tp = WeatherTemp;
pcfactor = 1/ProteinTotalRatio;

if (strlength(GRNFile) > 0)
    GRN_data = pcfactor * ReadGRN(GRNFile);
    VfactorC = GRN_data(1:33);
    if GRNC == 1
        cATPsyn = GRN_data(34);
        CPSi = GRN_data(35);
        cNADPHsyn = GRN_data(37);
        cpsii = GRN_data(36);
    else
        cATPsyn = 1;
        CPSi = 1;
        cNADPHsyn = 1;
        cpsii = 1;
    end
else
    GRN_data = ones(36, 1);
end
if (strlength(EnzymeFile) > 0)
    Enzyme_data = ReadEnzymes(EnzymeFile);
    global EnzymeAct;
    EnzymeAct = Enzyme_data(1:27) / 30; % unit change
    global Jmax;
    Jmax = EnzymeAct(27);
    global BFVmax;
    BFVmax = Enzyme_data(28:45);
    global FIVmax;
    FIVmax = Enzyme_data(46:66);
else
    Enzyme_data = ones(66, 1);
end

if OutputParam
    global export_mod_enabled;
    export_mod_enabled = 1;
end

fileID = 0;
if strlength(OutputFile) > 0
    fileID = fopen(OutputFile, 'w');
end

if (driver == 1)
    ResultRate = trDynaPS_Drive(CO2i, PPFDi, ATPCost, SucPath, 1, 1);
    if (fileID)
        Output(1)="Light intensity,Vc,Vo,VPGA,VT3P,Vstarch,Vt_glycerate,Vt_glycolate,CO2AR";
        Output(2)=Lii+","+ResultRate(1)+","+ResultRate(2)+","+ResultRate(3)+","+ResultRate(4)+","+ResultRate(5)+","+ResultRate(6)+","+ResultRate(7)+","+ResultRate(8);
        fprintf(fileID,Output(1));
        fprintf(fileID,"\n");
        fprintf(fileID,Output(2));
    end
elseif (driver == 2)
elseif (driver == 3)
    % CO2AR = CM_Driver(pop, currentPop)
elseif (driver == 4)
    % if (strlength(GRNFile) > 0 || strlength(EnzymeFile) > 0)
    %     Arate = EPS_Drive_GRNs(GRN_data, CO2i, PPFDi, WeatherTemp, ...
    %                            GRNC, GRNT, Enzyme_data);
    % else
    Arate = EPS_Drive(PPFDi, CO2i, WeatherTemp);
    % end
    if (fileID)
        fprintf(fileID, "%f\n", Arate);
    end
end
if fileID
    fclose(fileID);
end
