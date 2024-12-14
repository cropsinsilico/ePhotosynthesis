mfilePath = mfilename('fullpath');
if contains(mfilePath,'LiveEditorEvaluationHelper')
    mfilePath = matlab.desktop.editor.getActiveFilename;
end
mfileParts = fileparts(mfilePath);
topDir = [fileparts(mfileParts) filesep];
addpath(genpath(topDir));
test_dir = [mfileParts filesep];
data_dir = [test_dir 'data' filesep];
EnvFile = [data_dir 'InputEvn.txt'];
EnvFile_soyC3 = [data_dir 'InputEvn_soyC3.txt'];
GRNFile = [data_dir 'InputGRNC.txt'];
GRNFile_soyC3 = [data_dir 'InputGRNC_soyC3.txt'];
EnzymeFile = [data_dir 'InputEnzyme.txt'];
ATPCostFile = [data_dir 'InputATPCost.txt'];
OutputFile = [test_dir 'output.data'];
EnvFile0 = [topDir 'InputEvn.txt'];
GRNFile0 = "";
EnzymeFile0 = "";
ATPCostFile0 = [topDir 'InputATPCost.txt'];

%% Test EPS Driver with default input
expected = [23.4005893092884];
actual = ePhotosynthesis(4, EnvFile, GRNFile, EnzymeFile, ...
                         ATPCostFile, OutputFile);
assert(compare_floats(actual, expected));
expectedOutput = [data_dir 'EPS_default_output.data'];
assert(compare_files(OutputFile, expectedOutput));

%% Test EPS Driver with C3 input
% expected = [23.4005807699765];  % with PAR_in_Wpm2 = 0
expected = [25.8488652511769];  % with PAR_in_Wpm2 = 1
actual = ePhotosynthesis(4, EnvFile_soyC3, GRNFile_soyC3, EnzymeFile, ...
                         ATPCostFile, OutputFile);
assert(compare_floats(actual, expected));
expectedOutput = [data_dir 'EPS_soyC3_output.data'];
assert(compare_files(OutputFile, expectedOutput));

%% Test EPS Driver with output parameters
OutputParamBase = [data_dir 'actualParam_'];
expectedParamInit = [data_dir 'EPS_C3_init.txt'];
expectedParamLast = [data_dir 'EPS_C3_last.txt'];
actualParamInit = [OutputParamBase 'init.txt'];
actualParamLast = [OutputParamBase 'last.txt'];
% expected = [23.4005807699765];  % with PAR_in_Wpm2 = 0
expected = [25.8488652511769];  % with PAR_in_Wpm2 = 1
actual = ePhotosynthesis(4, EnvFile_soyC3, GRNFile_soyC3, EnzymeFile, ...
                         ATPCostFile, OutputFile, OutputParamBase);
assert(compare_floats(actual, expected));
assert(compare_files(actualParamInit, expectedParamInit));
assert(compare_files(actualParamLast, expectedParamLast));
expectedOutput = [data_dir 'EPS_soyC3_output.data'];
assert(compare_files(OutputFile, expectedOutput));
