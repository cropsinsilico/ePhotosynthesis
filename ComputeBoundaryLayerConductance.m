function BoundaryCound=ComputeBoundaryLayerConductance(WeatherTemperature,LeafTemperature,WeatherRH,WeatherWind,Pressure,LeafGs,NetAssimilation, Air_CO2)

%ComputeBoundaryLayerConductance(28,26,0.6,1,101325.0,0.01,20,380)
%Temporary variable
 ESaturation = 0.0; Ei = 0.0; Ea = 0.0; % Vapor pressure Pa
 LeafDimension = 0.06; %Leaf width/needle diameter m
 CForced = 4.322 / 1000.0; 
 CFree = 1.6361 / 1000.0;
 TemperatureKelvin = WeatherTemperature + 273.15;% Air temperature K
 LeafTemperatureKelvin = LeafTemperature + 273.15;%Leaf temperature K
 ESaturation = 0.611 * exp(17.502 * LeafTemperature / (LeafTemperature + 240.97)) * 1000.0;%Vapor pressure Pa
 Ei = ESaturation; %Vapor pressure Pa
 Ea = WeatherRH * ESaturation; %Vapor pressure Pa

 %Forced convection
 GbForced = CForced * TemperatureKelvin^0.56 * sqrt((TemperatureKelvin + 120.0)* (WeatherWind / LeafDimension / Pressure)); %// m/s

 %Free convection
 GbFree = GbForced;% m/s
 Eb = (LeafGs / 41.1 * Ei + GbFree * Ea) / (LeafGs / 41.1 + GbFree); % Stomatal conductance from moles/m2 leaf area/s to m/s
 TDifference = (LeafTemperatureKelvin / (1.0 - 0.378 * Eb /Pressure)) -(TemperatureKelvin / (1.0 - 0.378 * Ea / Pressure));
 GbFree = CFree * LeafTemperatureKelvin^0.56 * sqrt((LeafTemperatureKelvin + 120.0) /Pressure) * (abs(TDifference) / LeafDimension)^0.25;% m/s

% Maximum of two conductances
 if GbFree >= GbForced
     Gb = GbFree;
 else
     Gb = GbForced; % m/s
 end
 Gb = Gb * 41.4; %Conversion from m/s to moles/m2 leaf area/s
 Cb = Air_CO2 - 1.37 * NetAssimilation / Gb; % ppm
 BoundaryCound(1)=Gb;
 BoundaryCound(2)=Cb;
 BoundaryCound(3)=Eb;
end
