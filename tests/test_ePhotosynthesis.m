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

%% Test trDynaPS Driver with default input
expected = [470.0; 21.6102; 7.2934; 0.025743; 4.1691; ...
            0.89649; 3.6467; 7.2934; 3.0177];
actual = ePhotosynthesis(1, EnvFile, GRNFile, EnzymeFile, ...
                         ATPCostFile, OutputFile);
reldiff = abs(actual - expected) ./ actual
absdiff = abs(actual - expected)
assert(all(reldiff <= reltol))
assert(all(absdiff <= abstol))

%% Test DynaPS Driver with default input
expected = [470.0; 124.067; 160.8063; 21.7496; 7.3405; ...
            0.025081; 3.3163; 0.94269; 3.1981];
actual = ePhotosynthesis(2, EnvFile, GRNFile, EnzymeFile, ...
                         ATPCostFile, OutputFile);
reldiff = abs(actual - expected) ./ actual
absdiff = abs(actual - expected)
assert(all(reldiff <= reltol))
assert(all(absdiff <= abstol))

%% Test CM Driver with default input
expected = [470.0; 17.6567];
actual = ePhotosynthesis(3, EnvFile, GRNFile, EnzymeFile, ...
                         ATPCostFile, OutputFile);
reldiff = abs(actual - expected) ./ actual
absdiff = abs(actual - expected)
assert(all(reldiff <= reltol))
assert(all(absdiff <= abstol))

%% Test EPS Driver with default input
expected = [14.6712];
actual = ePhotosynthesis(4, EnvFile, GRNFile, EnzymeFile, ...
                         ATPCostFile, OutputFile);
reldiff = abs(actual - expected) ./ actual
absdiff = abs(actual - expected)
assert(all(reldiff <= reltol))
assert(all(absdiff <= abstol))


%% Test trDynaPS Driver with Sim_Ephotosynthesis inputs
expected = [800.0; 24.1202; 8.1406; 0.0040041; 1.5878; ...
            2.5457; 4.0703; 8.1406; 5.9703];
actual = ePhotosynthesis(1, EnvFile0, GRNFile0, EnzymeFile0, ...
                         ATPCostFile0, OutputFile);
reldiff = abs(actual - expected) ./ actual
absdiff = abs(actual - expected)
assert(all(reldiff <= reltol))
assert(all(absdiff <= abstol))
