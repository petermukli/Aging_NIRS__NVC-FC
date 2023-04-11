%%% MAIN SCRIPT for analysis of static functional connectivity %%%
% Written by Peter Mukli

%% Add path

addpath('Preprocess')
addpath('Analysis');
addpath('Analysis\2015_01_25 BCT')
load('channel_info.mat');

mypath = pwd;

%% Declarations

Fs=3.90625;  % sampling frequency

% bandpass filtering
l_p = 0.0045; % lower cutoff
h_p = 0.4; % upper cutoff
fltn = 'B4_-_B1'; 

% preprocess
mn_p = {'cui_filt'};
% mn_p = {'nocui_nofilt', 'cui_nofilt', 'nocui_filt', 'cui_filt'};
mn_m = {'pearson'};
% mn_m = {'pearson', 'spearman', 'hilbert_r', 'hilbert_pli' 'hilbert_c'};

igotthis = dir('.\Data_raw\CN*');
ID = zeros(1, length(igotthis));
for i = 1 : length(igotthis)
    ID(i) = str2double(igotthis(i).name(3:end));
end
clear i igotthis

%% Preprocess intensity data, convert into chromofor signals


% Preallocations
alany = struct('age', [], 'path_760', [], 'path_850', [], ...
    'nback', struct('raw_time', [], 'int760', [], 'int850', [], ...
    'int760_reg', [], 'int850_reg', []));

chromofor = struct('nback', struct('hbr', [], 'hbo', [], ...
    'hbr_cwfilt', [], 'hbo_cwfilt', [], 'hbr_bpfilt', [], 'hbo_bpfilt', [], ...
    'hbr_reg', [], 'hbo_reg', []));

% Cycle through subjects
for id = ID
    try
        disp([int2str(id), ' | ', int2str(length(ID))])
    
		% find and load intensity files
		
        mydir = ['CN', sprintf('%03d', id)];
        cd(['.\Data_raw\', mydir])
        fn.wl760 = dir([pwd, './*.wl1']);
        fn.wl850 = dir([pwd, './*.wl2']);
		
		% find and load marker files
        mrk_fn = dir('*.evt');
        mrk = load(mrk_fn(1).name);
        temp_wl760 = load(fn.wl760(1).name);
        temp_wl850 = load(fn.wl850(1).name);
        cd(mypath)

		% Correction for age-related change in DPF
        cella{1} = ['B', int2str(id + 2)];
        cella{2} = ['C', int2str(id + 2)];
        cella{3} = ['D', int2str(id + 2)];
		
		% load DPF values from a separate excel file
        alany(id).age = xlsread('.\DPF_calculator_Scholkmann.xlsx', 'Sheet1', cella{1});
        alany(id).path_760 = xlsread('.\DPF_calculator_Scholkmann.xlsx', 'Sheet1', cella{2});
        alany(id).path_850 = xlsread('.\DPF_calculator_Scholkmann.xlsx', 'Sheet1', cella{3});
        alany(id).NCh = 48;

        % cycle through n-back sessions
        for session =  1 : 4
        
            if session < 4
                alany(id).nback(session).raw_time = (mrk(session, 1) : mrk(session + 1, 1)) ./ Fs;
                % to work with intensity data, wavelength=760 nm
				alany(id).nback(session).int760_reg = ...
                    temp_wl760(mrk(session, 1) : mrk(session + 1, 1), relevant_48channels(:, 3));
                % to work with intensity data, wavelength=850 nm
                alany(id).nback(session).int850_reg = ...
                    temp_wl850(mrk(session, 1) : mrk(session + 1, 1), relevant_48channels(:, 3));                
            else
                alany(id).nback(session).raw_time = (mrk(session, 1) : size(temp_wl760, 1)) ./ Fs;
                % to work with intensity data, wavelength=760 nm
                alany(id).nback(session).int760 = ...
                    temp_wl760(mrk(session, 1) : end, relevant_48channels(:, 3));
                % to work with intensity data, wavelength=850 nm
                alany(id).nback(session).int850 = ...
                    temp_wl850(mrk(session, 1) : end, relevant_48channels(:, 3));
            end
            
            
            %%% Quality control based on coefficient of variation 
            %%% (>7.5% rejected)
            for c = 1 : 48
                qc_760 = std(alany(id).nback(session).int760(:, c)) / ...
                    mean(alany(id).nback(session).int760(:, c));
                qc_850 = std(alany(id).nback(session).int850(:, c)) / ...
                    mean(alany(id).nback(session).int850(:, c));
                alany(id).nback(session).QC(c) = (qc_760<0.075) && (qc_850<0.075);
                clear qc*
            end
          
		% Preprocessing employing discrete wavelete transform, modified Beer-Lambert Law
          [results, raw_output] = ...
            calc_int_to_hemo_dw(alany(id).nback(session), alany(id).NCh, Fs, ...
                l_p, h_p, alany(id).path_760, alany(id).path_850);

            alany(id).nback(session).int760_dwfilt = raw_output.int760_dwfilt;
            alany(id).nback(session).int760_bpfilt = raw_output.int760_bpfilt;
            alany(id).nback(session).int850_dwfilt = raw_output.int850_dwfilt;
            alany(id).nback(session).int850_bpfilt = raw_output.int760_bpfilt;
            
            clear raw_output

			% Assigning chromofores as outputs
            chromofor(id).nback(session).hbr = results.hbr; % deoxyhemoglobin (raw)
            chromofor(id).nback(session).hbo = results.hbo; % oxyhemoglobin (raw)
            chromofor(id).nback(session).hbt = ...
                results.hbr + results.hbo; % total hemoglobin (raw)
            chromofor(id).nback(session).hbr_dwfilt = results.hbr_dwfilt; % HbR after discret wavelet filtering
            chromofor(id).nback(session).hbo_dwfilt = results.hbo_dwfilt; % HbO after discret wavelet filtering
            chromofor(id).nback(session).hbt_dwfilt = ...
                results.hbr_dwfilt + results.hbo_dwfilt; % HbT after discret wavelet filtering
            chromofor(id).nback(session).hbr_bpfilt = results.hbr_bpfilt; % HbR after bandpass filtering
            chromofor(id).nback(session).hbo_bpfilt = results.hbo_bpfilt; % HbO after bandpass filtering
            chromofor(id).nback(session).hbt_bpfilt = ...
                results.hbr_bpfilt + results.hbo_bpfilt; % HbT after bandpass filtering

            clear results raw_output

    
        end
		

		% Repeat preprocessing pipeline, as above, but here after removing short channel signal
        % cycle through n-back sessions

		for session = 1 : 4

            if session < 4
                alany(id).nback(session).raw_time = (mrk(session, 1) : mrk(session + 1, 1)) ./ Fs;
                % to work with short-channel corrected intensity data, wavelength=760 nm
				alany(id).nback(session).int760_reg = ...
                    temp_wl760(mrk(session, 1) : mrk(session + 1, 1), relevant_48channels(:, 3)) - ...
					alany(id).nback(session).int760(mrk(session, 1) : mrk(session + 1, 1), alany(id).NCh + 1);
                % to work with short-channel corrected intensity data, wavelength=850 nm
                alany(id).nback(session).int850_reg = ...
                    temp_wl850(mrk(session, 1) : mrk(session + 1, 1), relevant_48channels(:, 3)) - ...
					alany(id).nback(session).int850(mrk(session, 1) : mrk(session + 1, 1), alany(id).NCh + 1);					;                
            else
                alany(id).nback(session).raw_time = (mrk(session, 1) : size(temp_wl760, 1)) ./ Fs;
                % to work with short-channel corrected intensity data, wavelength=760 nm
                alany(id).nback(session).int760 = ...
                    temp_wl760(mrk(session, 1) : end, relevant_48channels(:, 3)) - ...
					alany(id).nback(session).int760(mrk(session, 1) : end, alany(id).NCh + 1);
				% to work with short-channel corrected intensity data, wavelength=850 nm
                alany(id).nback(session).int850 = ...
                    temp_wl850(mrk(session, 1) : end, relevant_48channels(:, 3)) - ...
					alany(id).nback(session).int850(mrk(session, 1) : end, alany(id).NCh + 1);
            end
					
					 					 
			% Preprocessing employing discrete wavelete transform, modified Beer-Lambert Law
			[results, raw_output] = ...
				calc_int_to_hemo_dw(alany(id).nback(session), alany(id).NCh, ...
					Fs, l_p, h_p, alany(id).path_760, alany(id).path_850);            
				

			alany(id).nback(session).int760_dwfilt_reg = raw_output.int760_dwfilt_reg;
			alany(id).nback(session).int760_bpfilt_reg = raw_output.int760_bpfilt_reg;
			alany(id).nback(session).int850_dwfilt_reg = raw_output.int850_dwfilt_reg;
			alany(id).nback(session).int850_bpfilt_reg = raw_output.int760_bpfilt_reg;
				
			clear raw_output            
			
			% Assigning chromofores as outputs (short-channel corrected)			
			chromofor(id).nback(session).hbr_reg = results.hbr_reg; % deoxyhemoglobin (raw)
			chromofor(id).nback(session).hbo_reg = results.hbo_reg; % oxyhemoglobin (raw)
			chromofor(id).nback(session).hbt_reg = ...
				results.hbr_reg + results.hbo_reg;  % total oxyhemoglobin (raw)
			chromofor(id).nback(session).hbr_dwfilt_reg = results.hbr_dwfilt_reg;  % HbR after discret wavelet filtering
			chromofor(id).nback(session).hbo_dwfilt_reg = results.hbo_dwfilt_reg; % HbO after discret wavelet filtering
			chromofor(id).nback(session).hbt_dwfilt_reg = ...
				results.hbr_dwfilt_reg + results.hbo_dwfilt_reg; % HbT after discret wavelet filtering
			chromofor(id).nback(session).hbr_bpfilt_reg = results.hbr_bpfilt; % HbR after bandpass filtering
			chromofor(id).nback(session).hbo_bpfilt_reg = results.hbo_bpfilt; % HbO after bandpass filtering
			chromofor(id).nback(session).hbt_bpfilt_reg = ...
				results.hbr_bpfilt_reg + results.hbo_bpfilt_reg; % HbT after bandpass filtering
				
			clear results raw_output
		
		
		end
		
		
		% Identifying channels that pass quality control tests - each channel needs to be "good" in all n-back sessions
		for c = 1 : 48    
			alany(id).QC(c) = alany(id).nback(1).QC(c) & alany(id).nback(2).QC(c)  & ...
				alany(id).nback(3).QC(c) & alany(id).nback(4).QC(c); 
		end

    % error control
    catch ME
        continue
    end
end


%%%

ID2 = [];
j=0;
for i = 1 : length(chromofor)
    if size(chromofor(i).nback, 2)==4
        j = j + 1;
        ID2(j)=i;
    end
end


%% Preallocations %%


% fields indicate binary and weighted network metrics: density, clustering coefficient and efficiency
konnekt_param = struct('density_bu', [], 'density_wu', [], ...
    'clustering_bu', [], 'clustering_wu', [], ...
    'efficiency_bu', [], 'efficiency_wu', []); %cell(2,1)

metrika2 = struct(mn_m{1}, konnekt_param);

% metrika2 = struct('pearson', konnekt_param, 'spearman', konnekt_param, ...
%    'hilbert_r', konnekt_param, ...
%    'hilbert_pli', konnekt_param, 'hilbert_c', konnekt_param);

nback = repmat(repmat(struct('connection_matrix', []), 1, 4), j, 1);

mynetwork = repmat(struct('nback', repmat(konnekt_param, 1, 4)), j, 1);

% mynetwork = struct('nback', struct('nocui_nofilt', metrika2, ...
%     'cui_nofilt', metrika2, 'nocui_filt', metrika2, 'cui_filt', metrika2));

clear j


for j = 1 : length(ID2)
    
    try

        mydir = ['.\Data_raw\CN', sprintf('%03d', j)];
        % determining channels that pass quality control
		vec = 1 : 48;
        vec_temp = find(alany(ID2(j)).QC == 0);
        for ii = 1 : length(vec_temp)
            vec(vec==vec_temp(ii))=[];
        end
		
		% call function: analysis of static functional connectivity
        myresult = Adjacency_calc__Pearson_cui_filt(chromofor, ID2(j), vec, 1);
		
		
		% assign adjacency matrices
        nback(j, 1).connection_matrix = myresult.nback(1); 
        nback(j, 2).connection_matrix = myresult.nback(2);  
        nback(j, 3).connection_matrix = myresult.nback(3);  
        nback(j, 4).connection_matrix = myresult.nback(4);  

        close all
        clear myresult

        mynetwork(j).age = alany(ID2(j)).age;
        mynetwork(j).vec = vec;

		% cycle through sessions for network analysis
		% calling functions of Brain Connectivity Toolbox
		
		% Rubinov, Mikail, and Olaf Sporns. 
		% “Complex network measures of brain connectivity: uses and interpretations.” 
		% NeuroImage vol. 52,3 (2010): 1059-69. doi:10.1016/j.neuroimage.2009.10.003
		
        for session = 1 : 4

            disp(['id: ', int2str(j), '; session: ',...
                int2str(session)])
                        
			% Normalized local node degree ~ network density
            mynetwork(j).nback(session).density_bu{1} = ...
                degrees_und(nback(j, session).connection_matrix.pearson) ./ length(vec);
            
			% Normalized global node degree ~ network density
			[mynetwork(j).nback(session).density_bu{2}, ~, ~] = ...
                density_und(nback(j, session).connection_matrix.pearson); 

			% Normalized local connection strength ~ network density
            [mynetwork(j).nback(session).density_wu{1}] = ...
                strengths_und(nback(j, session).connection_matrix.pearson) ./ length(vec);

			% Normalized global connection strength ~ network density
            mynetwork(j).nback(session).density_wu{2} = ...
                mean(mynetwork(j).nback(session).density_wu{1});

            temp_bin = double(nback(j, session).connection_matrix.pearson~=0);

			% Local clustering coefficient (binary)
            mynetwork(j).nback(session).clustering_bu{1} = ...
                clustering_coef_bu(temp_bin);
			% Global clustering coefficient (binary)
            mynetwork(j).nback(session).clustering_bu{2} = ...
                mean(mynetwork(j).nback(session).clustering_bu{1});
			% Local clustering coefficient (weighted)
            mynetwork(j).nback(session).clustering_wu{1} = ...
                clustering_coef_wu(nback(j, session).connection_matrix.pearson);
			% Global clustering coefficient (weighted)
            mynetwork(j).nback(session).clustering_wu{2} = ...
                mean(mynetwork(j).nback(session).clustering_wu{1});

			% Local efficiency (binary)
            mynetwork(j).nback(session).efficiency_bu{1} = ...
                efficiency_bin(temp_bin, 1);
			% Global efficiency (weighted)
            mynetwork(j).nback(session).efficiency_bu{2} = ...
                efficiency_bin(temp_bin);

			% Local efficiency (weighted)
            mynetwork(j).nback(session).efficiency_wu{1} = ...
                efficiency_wei(temp_bin, 1);
			% Global efficiency (weighted)
            mynetwork(j).nback(session).efficiency_wu{2} = ...
                efficiency_wei(temp_bin);

            clear temp_bin

        end
    catch ME
        alany(ID2(j)).hiba = ME;
        continue
    end
end

save([pwd, '.\NIRS_preprocessed_individual\', 'NIRS_data__'preprocessed'], ...
    'alany',  'chromofor')
	
%
save([pwd, '.\', 'CN_study_StaticConnectivity__Pearson__', fltn], ...
    'mynetwork', 'nback')