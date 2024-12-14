abstol = 1.0e-02;
reltol = 1.0e-03;
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
GRNFile = [data_dir 'InputGRNC.txt'];
EnzymeFile = [data_dir 'InputEnzyme.txt'];
ATPCostFile = [data_dir 'InputATPCost.txt'];
OutputFile = [test_dir 'output.data'];
EnvFile0 = [topDir 'InputEvn.txt'];
GRNFile0 = "";
EnzymeFile0 = "";
ATPCostFile0 = [topDir 'InputATPCost.txt'];

%% Test EPS Driver with default input
expected = [29.4437170091883];
actual = ePhotosynthesis(4, EnvFile, GRNFile, EnzymeFile, ...
                         ATPCostFile, OutputFile);
reldiff = abs(actual - expected) ./ actual
absdiff = abs(actual - expected)
assert(all(reldiff <= reltol))
assert(all(absdiff <= abstol))
