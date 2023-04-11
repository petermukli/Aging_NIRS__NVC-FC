function [ cleandata ] = func_wfilter_nirs( data)

% Wavelet filtering of signals
% Written by Frigyes Samuel Racz

if size(data,1) > size(data,2)
    data = data';
    mk_inv = 1;
else
    mk_inv = 0;
end

L = 5; % wavelet decomposition level

% padding with zeros for discrete wavelet transform
modulus = mod(length(data), 2^L);
if modulus ~=0
    extra = zeros(1,(2^L)-modulus);
else
    extra = [];
end

if ~isempty(extra)
    data_tmp = [data,extra]; % pad with zeros
else
    data_tmp = data;
end

% wavelet filtering
[thr, sorh, ~] = ddencmp('den','wv',data);
SWC = swt(data_tmp, L, 'coif5');
SWC_th = wthresh(SWC, sorh, thr);
wdata_tmp = iswt(SWC_th, 'coif5');

% remove padding
if ~isempty(extra)
    wdata_tmp = wdata_tmp(:, 1:end-numel(extra));
end

cleandata = data - wdata_tmp;

if mk_inv
    cleandata = cleandata';
end

end

