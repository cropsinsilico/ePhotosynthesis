function InitRatios(ParaNum, Ratio)
if ParaNum <= 0
  return;
end
  
global PSRatio;
PSRatio=ones(103,1); 
if ParaNum <= 103
PSRatio(ParaNum)=Ratio;
end

global SUCRatio;
SUCRatio=ones(66,1);
if ParaNum > 103&&ParaNum<=169
SUCRatio(ParaNum-103)=Ratio;
end

global PRRatio;
PRRatio=ones(48,1);
if ParaNum > 169&&ParaNum<=217
PRRatio(ParaNum-169)=Ratio;
end

global RacRatio;
RacRatio=ones(16,1);
if ParaNum > 217&&ParaNum<=233
RacRatio(ParaNum-217)=Ratio;
end

global FIRatio;
FIRatio=ones(23,1);
if ParaNum > 233&&ParaNum<=256
FIRatio(ParaNum-233)=Ratio;
end

global BFRatio;
BFRatio=ones(49,1);
if ParaNum > 256&&ParaNum<=305
BFRatio(ParaNum-256)=Ratio;
end    

global XanRatio;
XanRatio=ones(4,1);
if ParaNum > 305&&ParaNum<=309
XanRatio(ParaNum-305)=Ratio;
end

end
