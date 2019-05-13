function CalGs=ComputGsBallBerry(BallBerryIntercept,BallBerrySlope,WaterStressFactor,Eb,Cb,Gb,NetAssimilation,LeafTemp,Air_CO2)
%ComputGsBallBerry(0.008,10.5,1.0,1000,10,345,25,25,370)
%%Saturated water vapor pressure
ESaturation=0.611 * exp(17.502 * LeafTemp/ (LeafTemp+240.97)) * 1000.0;

%%Update stomatal conductance moles/m2 leaf area/s
Gs = BallBerryIntercept + BallBerrySlope*NetAssimilation * Eb / ESaturation/Cb;
Gs = Gs * WaterStressFactor; %Apply water stress factor

%     if (LeafState->Gs < Photosynthesis->BallBerryIntercept || Photosynthesis->CO2 < LeafState->GammaStar)
%     {
%         LeafState->Gs = Photosynthesis->BallBerryIntercept;
%     }

%% Update Ci
Ci = Air_CO2 - 1.37 * NetAssimilation / Gb - 1.6 * NetAssimilation / Gs; % Update leaf internal CO2 concentration ppm
%     if (LeafState->Ci <= LeafState->GammaStar)
%     {
%         LeafState->Ci = LeafState->GammaStar; // ppm
%     }
CalGs(1)=Gs;
CalGs(2)=Ci;
end