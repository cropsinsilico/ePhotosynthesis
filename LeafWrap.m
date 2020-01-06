function LeafAssim=LeafWrap(Env, Einput, Eio)

  global pcfactor;
  pcfactor=1.0/Env('ProteinTotalRatio');

  LeafA = Leaf(Env('WeatherRH'), Env('WeatherTemperature'), Env('Air_CO2'), ...
	       Env('WeatherWind'), Env('Radiation_PAR'), Env('Radiation_NIR'), ...
	       Env('Radiation_LW'), Env('PhotosynthesisType'), Env('Vcmax25'), ...
	       Env('Jmax25'), Env('GRNC'), Einput{:, 2}, Eio{:, 2});

  LeafAssim = LeafA(2);
end
