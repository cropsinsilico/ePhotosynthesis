function LeafAssim=LeafWrap(Env, Einput, Eio)

  global pcfactor;
  pcfactor=1.0/Env('ProteinTotalRatio');

  % Rearrange Einput to use the same order as sample
  inE = importdata('MeM_input5_0.txt');
  GeneOrder = inE.textdata(2:size(inE.textdata, 1), 1);
  GeneNo = size(GeneOrder, 1);
  GeneMap = containers.Map(Einput{:, 1}, Einput{:, 2});
  ExpValue = ones(GeneNo, 1);
  for i=1:GeneNo
    if isKey(GeneMap, GeneOrder{i, 1})
      ExpValue(i, 1) = GeneMap(GeneOrder{i, 1});
    else
      fprintf("Gene %s is missing.\n", GeneOrder{i, 1});
    end
  end

  LeafA = Leaf(Env('WeatherRH'), Env('WeatherTemperature'), Env('Air_CO2'), ...
	       Env('WeatherWind'), Env('Radiation_PAR'), Env('Radiation_NIR'), ...
	       Env('Radiation_LW'), Env('PhotosynthesisType'), Env('Vcmax25'), ...
	       Env('Jmax25'), Env('GRNC'), ExpValue, Eio{:, 2});

  LeafAssim = LeafA(2);
end
