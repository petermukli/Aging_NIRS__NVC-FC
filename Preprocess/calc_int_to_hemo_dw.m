function [ results, varargout ] = calc_int_to_hemo_dw(nirdata, Nch, Fs, ...
    lf, hf, opl_760, opl_850)
% preprocessing pipeline for intensitities measured at 760 nm and 850 nm (NIRX NirScout)
% this pipeline includes modified Beer-Lambert law too
% Written by Frigyes Samuel Racz and Peter Mukli

int760_all = nirdata.int760(:, 1 : end-1);
int850_all = nirdata.int850(:, 1 : end-1);

% detrending

int760_det_all = zeros(size(int760_all));
int850_det_all = zeros(size(int850_all));
mean760_all = zeros(1, Nch);
mean850_all = zeros(1, Nch);

for ch = 1 : Nch
    p760 = polyfit((1:length(int760_all(:, ch)))', int760_all(:, ch), 1);
    p850 = polyfit((1:length(int850_all(:, ch)))', int850_all(:, ch), 1);

    int760_det_all(:, ch) = detrend(int760_all(:, ch));
    int850_det_all(:, ch) = detrend(int850_all(:, ch));
end

% wavelet filtering



for ch = 1 : Nch
    mean760_all(:, ch) = mean(int760_det_all(:, ch));
    mean850_all(:, ch) = mean(int850_det_all(:, ch));
end

int760_dwfilt_all = zeros(size(int760_det_all));
int850_dwfilt_all = zeros(size(int760_det_all));

for ch = 1 : Nch
    int760_dwfilt_all(:, ch) = ...
        func_wfilter_nirs(int760_det_all(:, ch)) + mean760_all(:, ch);
    int850_dwfilt_all(:, ch) = ...
        func_wfilter_nirs(int850_det_all(:, ch)) + mean850_all(:, ch);
end
% 'Expected NumOctaves to be a scalar with value <= 1.'


% bandpass filtering

for ch = 1 : Nch

    int760_det = int760_det_all(:, ch);
    int850_det = int850_det_all(:, ch);    
    int760_dwfilt = int760_dwfilt_all(:, ch);
    int850_dwfilt = int850_dwfilt_all(:, ch);
    
    mean760 = mean(int760_dwfilt);
    mean850 = mean(int850_dwfilt);
    
    try
        int760_bpfilt = bp_filter((int760_dwfilt), Fs, lf, hf, 5) + mean760;
        int850_bpfilt = bp_filter((int850_dwfilt), Fs, lf, hf, 5) + mean850;
    catch
        int760_bpfilt = int760_dwfilt;
        int850_bpfilt = int850_dwfilt;
    end

    % intensity preprocessing
    if p760(2) < 0
        int760 = int760_det + abs(p760(2))+1000;
    else
        int760 = int760_det + p760(2)+1000;
    end    
    int760(int760 == 0) = 1;
    int760 = abs(int760 ./ mean(int760));
    bs760 = abs(mean(int760));

    if p850(2) < 0
        int850 = int850_det + abs(p850(2))+1000;
    else
        int850 = int850_det + p850(2)+1000;
    end    
    int850(int850 == 0) = 1;
    int850 = abs(int850 ./ mean(int850));
    bs850 = abs(mean(int850));


    if p760(2) < 0
        int760_dwfilt = int760_dwfilt + abs(p760(2))+1000;
    else
        int760_dwfilt = int760_dwfilt + p760(2)+1000;
    end    
    int760_dwfilt(int760_dwfilt == 0) = 1;
    int760_dwfilt = abs(int760_dwfilt ./ mean(int760_dwfilt));
    bs760_dwfilt = abs(mean(int760_dwfilt));

    if p850(2) < 0
        int850_dwfilt = int850_dwfilt + abs(p850(2))+1000;
    else
        int850_dwfilt = int850_dwfilt + p850(2)+1000;
    end    
    int850_dwfilt(int850_dwfilt == 0) = 1;
    int850_dwfilt = abs(int850_dwfilt ./ mean(int850_dwfilt));
    bs850_dwfilt = abs(mean(int850_dwfilt));

    if p760(2) < 0
        int760_bpfilt = int760_bpfilt + abs(p760(2))+1000;
    else
        int760_bpfilt = int760_bpfilt + p760(2)+1000;
    end    
    int760_bpfilt(int760_bpfilt == 0) = 1;
    int760_bpfilt = abs(int760_bpfilt ./ mean(int760_bpfilt));
    bs760_bpfilt = abs(mean(int760_bpfilt));


    if p850(2) < 0
        int850_bpfilt = int850_bpfilt + abs(p850(2))+1000;
    else
        int850_bpfilt = int850_bpfilt + p850(2)+1000;
    end    
    int850_bpfilt(int850_bpfilt == 0) = 1;
    int850_bpfilt = abs(int850_bpfilt ./ mean(int850_bpfilt));
    bs850_bpfilt = abs(mean(int850_bpfilt));

    % calculating hemodynamics using 760 & 850 nm - differential modified Beer-Lambert law


    [HbR_760_850, HbO_760_850] = calc_hemo_760_850(...
        bs760, bs850, int760, int850, ...
        opl_760, opl_850);

    [HbR_dwfilt_760_850, HbO_dwfilt_760_850] = calc_hemo_760_850(...
        bs760_dwfilt, bs850_dwfilt, int760_dwfilt, int850_dwfilt, ...
        opl_760, opl_850);

    [HbR_bpfilt_760_850, HbO_bpfilt_760_850] = calc_hemo_760_850(...
        bs760_bpfilt, bs850_bpfilt, int760_bpfilt, int850_bpfilt,...
        opl_760, opl_850);



    raw.int760_dwfilt(:, ch) = int760_dwfilt;
    raw.int760_bpfilt(:, ch) = int760_bpfilt;
    raw.int850_dwfilt(:, ch) = int850_dwfilt;
    raw.int850_bpfilt(:, ch) = int850_bpfilt;

    results.hbr(:, ch) = HbR_760_850;
    results.hbo(:, ch) = HbO_760_850;
    results.hbr_dwfilt(:, ch) = HbR_dwfilt_760_850;
    results.hbo_dwfilt(:, ch) = HbO_dwfilt_760_850;
    results.hbr_bpfilt(:, ch) = HbR_bpfilt_760_850;
    results.hbo_bpfilt(:, ch) = HbO_bpfilt_760_850;
    
    clear HbO_* HbR_* 
    clear int760_det int850_det int760_dwfilt int850_dwfilt 

end


varargout{1} = raw;


end





