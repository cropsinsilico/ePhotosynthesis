function values = ePhotosynthesis(driver, EnvFile, GRNFile, EnzymeFile, ATPCostFile, ...
                                 OutputFile, OutputParamBase, WeatherTemp, ...
                                 PARi, CO2i, O2i, ...
                                 ProteinTotalRatio, GRNCin, GRNTin, ...
                                 RUBISCOMETHODin, RUBISCOTOTALin, ...
                                 PAR_in_Wpm2)
arguments
    driver (1,1) string
    EnvFile (1,1) string = ""
    GRNFile (1,1) string = ""
    EnzymeFile (1,1) string = ""
    ATPCostFile (1,1) string = ""
    OutputFile (1,1) string = "output.data"
    OutputParamBase (1,1) string = ""
    WeatherTemp (1,1) double = 25
    PARi (1,1) double = 470
    CO2i (1,1) double = 400
    O2i (1,1) double = 0.21  % i.e. 21%
    ProteinTotalRatio (1,1) double = 0.973
    GRNCin (1,1) double = 0
    GRNTin (1,1) double = 0
    RUBISCOMETHODin (1,1) double = 2
    RUBISCOTOTALin (1,1) double = 3
    PAR_in_Wpm2 (1, 1) double = 0
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
global Temp_cond;
global TestLi;
global TestCa;
global TestATPCost;
global TestSucPath;
global LI;
global Tp;
global RUBISCOMETHOD;
global RUBISCOTOTAL;
global GLight;
global VolRatioStCyto;
global dontCondition;
global dontGraph;
global dontClear;
% dontCondition = 1;
dontGraph = 1;
% dontClear = 1;

SucPath = 0;
ATPCost = 0.0;
GRNC = GRNCin;
GRNT = GRNTin;
RUBISCOMETHOD = RUBISCOMETHODin;
RUBISCOTOTAL = RUBISCOTOTALin;
pcfactor = 1.0 / 0.973;

if (strlength(EnvFile) > 0)
    fprintf("EVN FILE PROVIDED: %s\n", EnvFile);
    Env_data = ReadParam(EnvFile);
    if isKey(Env_data, "WeatherTemp")
        WeatherTemp = Env_data("WeatherTemp");
    end
    if isKey(Env_data, "PAR")
        PARi = Env_data("PAR");
    end
    % TODO: Propery unit parsing?
    if isKey(Env_data, "PAR_in_Wpm2")
        PAR_in_Wpm2 = Env_data("PAR_in_Wpm2");
    end
    if isKey(Env_data, "CO2")
        CO2i = Env_data("CO2");
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
    fprintf("ATPCOST FILE PROVIDED: %s\n", ATPCostFile);
    ATPCost_data = ReadParam(ATPCostFile);
    if isKey(ATPCost_data, "ATPCost")
        ATPCost = ATPCost_data("ATPCost");
    end
end

TestCa = CO2i;
TestLi = PARi;
if (PAR_in_Wpm2)
    TestLi = TestLi * 1E6/(2.35E5); % Convert W m^{-2} to u moles m^{-2} s^{-1}
end
TestATPCost = ATPCost;
TestSucPath = SucPath;

if dontCondition
    CO2_Env = 0.7 * TestCa; % intercellular CO2
    CO2_cond = CO2_Env/(3 * 10^4);
    % TestLi = TestLi / 30;
    O2_cond = O2i * 1.26;
    Tp = WeatherTemp;
    pcfactor = 1/ProteinTotalRatio;
    Temp_cond = Tp;
    GLight = TestLi*0.85*0.85;  % light umol m-2 s-1
    VolRatioStCyto = 1;
end

if (strlength(GRNFile) > 0)
    fprintf("GRN FILE PROVIDED: %s\n", GRNFile);
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
    fprintf("ENZYME FILE PROVIDED: %s\n", EnzymeFile);
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

if (strlength(OutputParamBase) > 0)
    global export_mod_file_base;
    global export_mod_enabled;
    export_mod_enabled = 1;
    export_mod_file_base = OutputParamBase;
end

if (driver == 1)
    error("trDynaPS driver is disabled on this branch");
    % ResultRate = trDynaPS_Drive(TestCa, TestLi, TestATPCost, TestSucPath, 1, 1);
    % Arate = TargetFunVal;
    % names = ["Light intensity", "Vc", "Vo", "VPGA", "VT3P", ...
    %          "Vstarch", "Vt_glycerate", "Vt_glycolate", "CO2AR"];
    % nres = size(ResultRate, 1);
    % values = zeros(nres+2, 1);
    % values(1, 1) = TestLi;
    % for i=1:nres
    %     values(i+1, 1) = ResultRate(i);
    % end
    % values(nres+2, 1) = Arate;
elseif (driver == 2)
    error("DynaPS driver is disabled on this branch");
    % ResultRate = DynaPS_Drive(1, 1);
    % Arate = TargetFunVal;
    % names = ["Light intensity", "PSIIabs", "PSIabs", "Vc", "Vo", ...
    %          "VPGA", "Vstarch", "Vsucrose", "CO2AR"];
    % nres = size(ResultRate, 1);
    % values = zeros(nres+2, 1);
    % values(1, 1) = TestLi;
    % for i=1:nres
    %     values(i+1, 1) = ResultRate(i);
    % end
    % values(nres+2, 1) = Arate;
elseif (driver == 3)
    error("CM driver is disabled on this branch");
    % CM_Drive;
    % Arate = TargetFunVal;
    % names = ["Light intensity", "CO2AR"];
    % values = [TestLi; Arate];
elseif (driver == 4)
    Arate = EPS_Drive(TestLi, TestCa, WeatherTemp);
    names = [];
    values = [Arate];
end
WriteTable(OutputFile, names, values);
end
