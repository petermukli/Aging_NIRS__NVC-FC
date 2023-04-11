function [HbR, HbO] = calc_hemo_760_850(bs760, bs850, int760, int850, opl760, opl850)

% calculation of oxy- and deoxyhemoglobin from intensity data recorded at 760 nm and 850 nm
% differential modified Beer-Lambert Law

% Scott, P. S. J. (2015). Optical absorption and emission data. 
% Retrieved from the website:http://omlc.ogi.edu

% Written by Frigyes Samuel Racz and Peter Mukli
aHbR_760 = 1.54852;
aHbO_760 = 0.586;
aHbR_850 = 0.69132;
aHbO_850 = 1.058;

od760 = log10(repmat(bs760, [length(int760) 1]) ./ abs(int760));
od850 = log10(repmat(bs850, [length(int850) 1]) ./ abs(int850));

HbO = ((aHbR_760 * od850 ./ opl850) - (aHbR_850 * od760 ./ opl760)) ./ ...
    (aHbR_760 * aHbO_850 - aHbR_850 * aHbO_760);
HbR = ((aHbO_760 * od850 ./ opl850) - (aHbO_850 * od760 ./ opl760)) ./ ...
    (aHbO_760*aHbR_850 - aHbO_850 * aHbR_760);

end

