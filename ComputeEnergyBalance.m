%Compute energy balance
function CalLeafTemperature=ComputeEnergyBalance(Gb,Gs,NetAssimilation,LeafTemperature,RH,WeatherTemperature,Pressure,Boltzman,ConstantsCp,LatentHeatVaporization,Radiation_PAR,Radiation_NIR,Radiation_LW)

%ComputeEnergyBalance(10,10,25,26,0.6,28,101325.0,5.6697E-8,29.3,44000.0,200,0,0)
% Temporary variables
 MaximumLeafTemperature=80;
 Epsilon = 0.96;% Leaf thermal emissivity
 G = Gb * Gs / (Gb+Gs); % moles/m2 leaf area/s
 ESaturation = 0.0; Ea = 0.0; %Vapor pressure KPa
 LEFactor = 1.0; %%%%1.0%%%%%
 LWFactor = 2.0; 
 HFactor = 2.0;%%%%%GGGGGGGG
 EnergyBalanceResidual = 0.0;
 LeafEnergyFluxMe = 0.506 * NetAssimilation; %Energy in biochemical reactions W/m2
 ESaturation = 0.611 * exp(17.502 * LeafTemperature / (LeafTemperature + 240.97));%KPa
 Ea = RH * ESaturation;%KPa
 
     if Radiation_LW == 0
     Radiation_LW = LWFactor * Epsilon * Boltzman *(273.15 + WeatherTemperature) ^4;
     end


 CoefficientESaturation(1)=5.82436 / 10000.0 / 1000.0; % pow(LeafTemperature,4) KPa
 CoefficientESaturation(2)=1.5842 / 100.0 / 1000; %pow(LeafTemperature,3) KPa
 CoefficientESaturation(3)=1.55186 / 1000.0; %pow(LeafTemperature,2) KPa
 CoefficientESaturation(4)=44.513596 / 1000.0; %pow(LeafTemperature,1) KPa
 CoefficientESaturation(5)=607.919 / 1000.0; %pow(LeafTemperature,0) KPa

 CoefficientTemperature(1) = LEFactor * G * LatentHeatVaporization / (Pressure / 1000.0)*CoefficientESaturation(1) + LWFactor * Epsilon * Boltzman; %pow(LeafTemperature,4)
 CoefficientTemperature(2) = LEFactor * G * LatentHeatVaporization / (Pressure / 1000.0)*CoefficientESaturation(2) + 4 * LWFactor * Epsilon * Boltzman * 273.15; %pow(LeafTemperature,3)
 CoefficientTemperature(3) = LEFactor * G * LatentHeatVaporization / (Pressure / 1000.0)*CoefficientESaturation(3) + 6 * LWFactor * Epsilon * Boltzman * 273.15^2.0; %pow(LeafTemperature,2)
 CoefficientTemperature(4) = LEFactor * G * LatentHeatVaporization / (Pressure / 1000.0)*CoefficientESaturation(4) + 4 * LWFactor * Epsilon * Boltzman * 273.15^3.0 +HFactor * ConstantsCp * 0.924 * Gb; %pow(LeafTemperature,1)
 CoefficientTemperature(5) = LEFactor * G * LatentHeatVaporization / (Pressure / 1000.0)*(CoefficientESaturation(5) - Ea) + LWFactor * Epsilon * Boltzman * 273.15^4.0 -HFactor * ConstantsCp * 0.924 * Gb * WeatherTemperature - Radiation_PAR -Radiation_NIR - Radiation_LW + LeafEnergyFluxMe; %pow(LeafTemperature,0)

 %Check for NAN
    for Loop = 1:5
        if (isinf(CoefficientTemperature(Loop)) || isnan(CoefficientTemperature(Loop)))
            CoefficientTemperature(Loop)= 0;
        end
    end
    %Roots = rpoly(CoefficientTemperature, 4, RealRoots, ImaginaryRoots);
    global Roots;
    Roots=roots(CoefficientTemperature);
    %Select the lowest real root between -10 and 50
    LeafTemperature = MaximumLeafTemperature;
    for RootsLoop = 1:4
        if abs(abs(Roots(RootsLoop))-real(Roots(RootsLoop)))<0.0001&& Roots(RootsLoop) > -10.0 && LeafTemperature > Roots(RootsLoop)
            LeafTemperature = real(Roots(RootsLoop));
        end
    end
%     if LeafTemperature == MaximumLeafTemperature
% %        fprintf(LogOutputFile, "Error in ComputeEnergyBalance for Leaf ID %f: No real roots for leaf temperature=%f\n",
% %                 Photosynthesis->LeafID, LeafState->Temperature);
%        LeafTemperature = WeatherTemperature;
%     end

    %LeafState->Temperature = Weather->Temperature; %Uncomment to force leaf temperature to air temperature

    ESaturation = 0.611 * exp(17.502 * LeafTemperature / (LeafTemperature + 240.97)); %KPa
    Ea = RH * ESaturation; %KPa

    SensibleHeat = HFactor * ConstantsCp * 0.924 * Gb * (LeafTemperature - WeatherTemperature);
    LatentHeat = LEFactor * G * LatentHeatVaporization / (Pressure / 1000.0) * (ESaturation - Ea);
    Emission = LWFactor * Epsilon * Boltzman * (273.15 + LeafTemperature)^4.0;
    Transpiration = LatentHeat / LatentHeatVaporization * 10.0^6.0;

    EnergyBalanceResidual = Radiation_PAR + Radiation_NIR+Radiation_LW -Emission -SensibleHeat -LatentHeat - LeafEnergyFluxMe;

%     if (EnergyBalanceResidual/(Weather->Radiation.PAR + Weather->Radiation.NIR + Weather->Radiation.LW) > 0.01)
%     {
%         %fprintf(LogOutputFile, "Error in ComputeEnergyBalance for Leaf ID %f: Residual energy balance=%f W/m2\n",
%         %       Photosynthesis->LeafID, EnergyBalanceResidual);
%     }
    CalLeafTemperature(1)=LeafTemperature;
    CalLeafTemperature(2)=Transpiration;
    end
