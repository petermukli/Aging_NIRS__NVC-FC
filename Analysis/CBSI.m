function [HbO2_true,Hb_true]= CBSI(HbO2,Hb,NCh)
%input: HbO2,Hb and NCh
%output: HbO2_true and HbO2true
%HbO2: oxihaemoglobin concentration (in column)
%Hb: deoxihaemoglobin concentration (in column)
%NCh: number of channels
for i=1:NCh
oxy= HbO2(:,i);
deoxy= Hb(:,i);
alpha = std(oxy)./std(deoxy);
oxy0 = oxy - alpha .* deoxy;
oxy0 = oxy0 / 2;
HbO2_true(:,i)= oxy0;
Hb_true=-HbO2_true/alpha;
end

