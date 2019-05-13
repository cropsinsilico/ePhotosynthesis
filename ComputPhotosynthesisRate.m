function PhotosynthesisRate=ComputPhotosynthesisRate(PhotosynthesisType,PhotosynthesQ10,Vcmax25,Jmax25,Rd25,R,LeafTemperature,Convert, Radiation_PAR,PhotosynthesisTheta,Ci,Air_O2,GRNC,Einput,Eio)
%Temporary variables double Q10Temperature;
if (PhotosynthesisType == 1.0 || PhotosynthesisType == 1.1)% C3 Farquhar or Metabolic
    Rd = Rd25 * exp(18.72 - 46.39 / (R * (LeafTemperature + 273.15)));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Compute C3 photosynthesis
%Constant
if (PhotosynthesisType == 1.0)% C3 Farquhar or Metabolic
Rate_TPu = 23.0; %u moles/m2 leaf area/s
%Temporary variables
LeafTemperatureKelvin = LeafTemperature + 273.15; %Leaf temperature in K
GammaStar = exp(19.02 - 37.83 / (R * LeafTemperatureKelvin));
Ko = exp(20.30 - 36.38 / (R * LeafTemperatureKelvin));
Kc = exp(38.05 - 79.43 / (R * LeafTemperatureKelvin));	
Vcmax = Vcmax25 * exp(26.35 - 65.33 / (R * LeafTemperatureKelvin));
PhiPS2 = 0.385 + 0.02166 * LeafTemperature - 3.37 * LeafTemperature^2.0 / 10000.0;% Match PS_FIT
I = Convert * Radiation_PAR * PhiPS2 * 0.5;
ThetaPS2 = PhotosynthesisTheta + 0.01713 * LeafTemperature - 3.75 * LeafTemperature^2.0 / 10000.0; % Match PS_FIT
Jmax = Jmax25 * exp(17.57 - 43.54 / (R * LeafTemperatureKelvin));
J = (I + Jmax - sqrt((I + Jmax)^2.0 - 4.0 * ThetaPS2 * I * Jmax)) / (2.0 * ThetaPS2);
LeafAc = (1.0 - GammaStar / Ci) * (Vcmax * Ci) /(Ci + Kc * (1.0 + Air_O2 / Ko));%Rubisco limited photosynthesis
LeafAj = (1.0 - GammaStar / Ci) * (J * Ci) /(4.5 * Ci + 10.5 * GammaStar); %Light limited photosynthesis
if (LeafAj < 0.0)
    LeafAj = 0.0;
end
LeafAp = (3.0 * Rate_TPu) / (1.0 - GammaStar /Ci); %TPU limited photosynthesis
if (LeafAp < 0.0)
    LeafAp=0.0;
end
GrossAssimilation=min(min(LeafAc,LeafAj),LeafAp);%Minimum of three limitations
NetAssimilation=GrossAssimilation-Rd;
  
    if (isinf(GrossAssimilation) || isnan(GrossAssimilation))
%         fprintf(LogOutputFile, "Error in ComputeC3Photosynthesis for Leaf ID: %f: GrossAssimilation=%f\n",
%                 Photosynthesis->LeafID, LeafMassFlux->GrossAssimilation);
        GrossAssimilation = 0.0;
        NetAssimilation =-Rd;

    end
end

if (PhotosynthesisType == 1.1)% C3 metabolic
GrossAssimilation=EPS_Drive_GRNs(Einput,Ci,Convert * Radiation_PAR,LeafTemperature,GRNC,0,Eio);
NetAssimilation=GrossAssimilation-Rd;
    if (isinf(GrossAssimilation) || isnan(GrossAssimilation))
%         fprintf(LogOutputFile, "Error in ComputeC3Photosynthesis for Leaf ID: %f: GrossAssimilation=%f\n",
%                 Photosynthesis->LeafID, LeafMassFlux->GrossAssimilation);
        GrossAssimilation = 0.0;
        NetAssimilation =-Rd;

    end
LeafTemperatureKelvin = LeafTemperature + 273.15; %Leaf temperature in K
GammaStar = exp(19.02 - 37.83 / (R * LeafTemperatureKelvin));
end


PhotosynthesisRate(1)=NetAssimilation;
PhotosynthesisRate(2)=GrossAssimilation;
PhotosynthesisRate(3)=Rd;
PhotosynthesisRate(4)=GammaStar;
end