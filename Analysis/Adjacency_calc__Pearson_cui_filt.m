function result = Adjacency_calc__Pearson_cui_filt(chromofor, id, vec, varargin)
%%% Adjacency Matrix calculation
% Pearson correlation coefficient
% pipeline included CBSI
% surrogate thresholding, elimination of negative or insignificant coefficients

if nargin > 3
    fig = varargin{1};
end

task_name = {'0-back', '1-back', '0-back', '2-back'};
 
 % set absolute threshold value
th = 0;

% cycle through sessions
for session = 1 : 4

	% preallocations
    temp_v1_o = zeros(size(chromofor(id).nback(session).hbo_bpfilt_reg));
    temp_v1_r = zeros(size(chromofor(id).nback(session).hbr_bpfilt_reg));

	% bridge detrending
    for ch = 1 : 48
        temp_v1_r(:, ch) = bridge(chromofor(id).nback(session).hbr_bpfilt_reg(:, ch));
        temp_v1_o(:, ch) = bridge(chromofor(id).nback(session).hbo_bpfilt_reg(:, ch));
    end
    
	% Correlation-based signal improvement (Cui. et al. 2010)
    [temp_v1_o, temp_v1_r] = CBSI(temp_v1_o, temp_v1_r, 48);
    temp_v1 = temp_v1_o + temp_v1_r;
    
    temp_v1 = temp_v1(:, vec); 
	
	% calculating adjacency matrix
    [temp_pv1, p_temp_pv1] = corrcoef(temp_v1);
	% surrogate thresholding
    temp_pv1 = temp_pv1 .* (p_temp_pv1 < 0.05);

    clear p_temp*
        
    result.nback(session).pearson = temp_pv1;
	
	% absolute thresholding (eliminate anticorrelations)
    for k = 1 : length(th)
		
		% this is a function from the Brain Connectivity Toolbox (Rubinov and Sporns, 2010)
        result.nback(session).pearson = ...
            threshold_absolute(result.nback(session).pearson, th(k));

    end
        
    clear temp*

end