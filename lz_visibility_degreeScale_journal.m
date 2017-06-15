save('/home/lz206/Dropbox/projects/VG/VG_OneDrive/results_journal/power_sf_mat', 'power_sf_mat');
%%%% 05/22/2017 modified for the filepath

%%%% 05/25/2016 modified for journal

%%%% 01/02/2016
%%%% compute the degree distribution adn scale-freeness

function [kpk_log_tot, kpk_tot, power_sf] = lz_visibility_degreeScale_journal(subject) %, dropThreshold)

%%%% Input:
%%%% subject:       'm1n03, i6a08,i3i03,n8n10,a2u04,i3h06,i4n06,u8n09,l9e11'.
%%%% dropThreshold: the threshold set smaller than which the degree data was dropped from fitting. %%%% deleted on 02/10/2016

%%%% Output:
%%%% kpk_log_tot:   Cell array, Channel X Condition. The resultant power of scale-freeness.
%%%%                In each cell, the 1st column is the degree value, the 2nd
%%%%                column is the corresponding possibility of the degree value.
%%%% kpk_tot:       same format as kpk_log_tot, but the columns are k and p(k), for all k (without detecting the maximum point).
%%%% power_sf:      power of scale-freeness. Channel X Condition.

%==========================================================================
% load the visibility matrix
if ispc
    eval(['VG_ec  = load(''C:\Users\Li_Lab537\Dropbox\projects\VG\VG_OneDrive\results_journal\VG_',subject,'_ec_05222017_noDownSamp'');'])
else
    eval(['VG_ec  = load(''/home/lz206/Dropbox/projects/VG/VG_OneDrive/results_journal/VG_',subject,'_ec_05222017_noDownSamp'');'])
end
tmp = VG_ec; clear VG_ec;
eval(['VG_ec = tmp.VG_',subject,'_ec;'])

if ispc
    eval(['VG_eo  = load(''C:\Users\Li_Lab537\Dropbox\projects\VG\VG_OneDrive\results_journal\VG_',subject,'_eo_05222017_noDownSamp'');'])
else
    eval(['VG_eo  = load(''/home/lz206/Dropbox/projects/VG/VG_OneDrive/results_journal/VG_',subject,'_eo_05222017_noDownSamp'');'])
end
tmp = VG_eo; clear VG_eo;
eval(['VG_eo = tmp.VG_',subject,'_eo;'])

if ispc
    eval(['VG_rt  = load(''C:\Users\Li_Lab537\Dropbox\projects\VG\VG_OneDrive\results_journal\VG_',subject,'_rt_05222017_noDownSamp'');'])
else
    eval(['VG_rt  = load(''/home/lz206/Dropbox/projects/VG/VG_OneDrive/results_journal/VG_',subject,'_rt_05222017_noDownSamp'');'])
end
tmp = VG_rt; clear VG_rt;
eval(['VG_rt = tmp.VG_',subject,'_rt;'])

if ispc
    eval(['VG_mix = load(''C:\Users\Li_Lab537\Dropbox\projects\VG\VG_OneDrive\results_journal\VG_',subject,'_mix_05222017_noDownSamp'');'])
else
    eval(['VG_mix = load(''/home/lz206/Dropbox/projects/VG/VG_OneDrive/results_journal/VG_',subject,'_mix_05222017_noDownSamp'');'])
end
tmp = VG_mix; clear VG_mix;
eval(['VG_mix = tmp.VG_',subject,'_mix;'])

m = 14; % # of channels
n = size(VG_ec,1); % # of nodes

% initiate power of scale-freeness
degree_powerScale_ec  = zeros(m,3); % channel X block
degree_powerScale_eo  = zeros(m,3); % channel X block
degree_powerScale_rt  = zeros(m,3); % channel X block
degree_powerScale_mix = zeros(m,3); % channel X block

% compute degree for each node
% initiate degree
degree_eo  = zeros(n,m,3);
degree_ec  = zeros(n,m,3);
degree_rt  = zeros(n,m,3);
degree_mix = zeros(n,m,3);

for chani = 1: m % loop for channel
    for blocki =1: 3 % loop for block
        degree_eo(:,chani,blocki)  = squeeze(sum(VG_eo(:,:,chani,blocki)));
        degree_ec(:,chani,blocki)  = squeeze(sum(VG_ec(:,:,chani,blocki)));
        degree_rt(:,chani,blocki)  = squeeze(sum(VG_rt(:,:,chani,blocki)));
        degree_mix(:,chani,blocki) = squeeze(sum(VG_mix(:,:,chani,blocki)));
    end
end

% compute p(k) for each condition
% initiate k and p(k) matrix for each condition
kpk_ec  = cell(m,3);
kpk_eo  = cell(m,3);
kpk_rt  = cell(m,3);
kpk_mix = cell(m,3);
for chani = 1: m
    for blocki = 1: 3
        seq_bins = unique(degree_ec(:,chani,blocki));
        histcount = histc(degree_ec(:,chani,blocki), seq_bins);
%         kpk_ec  = cell(length(seq_bins),2,m,3); % node X 2 column X channel X blocks. Save values of k vs. p(k). First column: the degree value, second column: the associated p(k) value
        kpk_ec(chani,blocki) = {[seq_bins, histcount/n]};
%         kpk_ec(chani,blocki) = {[log2(1./seq_bins), log2(histcount/n)]};
            %%%% Save values of k vs. p(k). First column: log2(1/k), second column:  log2[p(k)]
%         kpk_ec{chani,blocki}(kpk_ec{chani,blocki}(:,1) > dropThreshold,:) = [];    
%             %%%% drop the entries that larger than the threshold
        
        seq_bins = unique(degree_eo(:,chani,blocki));
        histcount = histc(degree_eo(:,chani,blocki), seq_bins);
        kpk_eo(chani,blocki) = {[seq_bins, histcount/n]};
%         kpk_eo(chani,blocki) = {[log2(1./seq_bins), log2(histcount/n)]};
%         kpk_eo{chani,blocki}(kpk_eo{chani,blocki}(:,1) > dropThreshold,:) = [];     
        
        seq_bins = unique(degree_rt(:,chani,blocki));
        histcount = histc(degree_rt(:,chani,blocki), seq_bins);
        kpk_rt(chani,blocki) = {[seq_bins, histcount/n]};
%         kpk_rt(chani,blocki) = {[log2(1./seq_bins), log2(histcount/n)]};
%         kpk_rt{chani,blocki}(kpk_rt{chani,blocki}(:,1) > dropThreshold,:) = [];     
        
        seq_bins = unique(degree_mix(:,chani,blocki));
        histcount = histc(degree_mix(:,chani,blocki), seq_bins);
        kpk_mix(chani,blocki) = {[seq_bins, histcount/n]};
%         kpk_mix(chani,blocki) = {[log2(1./seq_bins), log2(histcount/n)]};
%         kpk_mix{chani,blocki}(kpk_mix{chani,blocki}(:,1) > dropThreshold,:) = [];
    end
end

% combine degree values across blocks for each channel and condition
% the following is the code developed based on http://www.mathworks.com/matlabcentral/answers/125787-merging-two-matrices-by-first-column-values
% initiate k vs. p(k) cell array
kpk_ec_tot  = cell(m,1);
kpk_eo_tot  = cell(m,1);
kpk_rt_tot  = cell(m,1);
kpk_mix_tot = cell(m,1);
% initiate log(1/k) vs. log[p(k)] cell array
kpk_log_ec_tot  = cell(m,1);
kpk_log_eo_tot  = cell(m,1);
kpk_log_rt_tot  = cell(m,1);
kpk_log_mix_tot = cell(m,1);


parfor chani = 1:m
    % ec
    q = unique([kpk_ec{chani,1}(:,1);kpk_ec{chani,2}(:,1);kpk_ec{chani,3}(:,1)]); % unique sorted catenate.
    Ind_blk1 = ismember(q, kpk_ec{chani,1}(:,1)); % All dataes in block 1
    Ind_blk2 = ismember(q, kpk_ec{chani,2}(:,1)); % All dataes in block 2
    Ind_blk3 = ismember(q, kpk_ec{chani,3}(:,1)); % All dataes in block 3
    x = zeros(length(q),4);
    x(:,1) = q; % value of k
    x(Ind_blk1,2) = kpk_ec{chani,1}(:,2); % p(k) in block 1
    x(not(Ind_blk1),2) = 0;
    x(Ind_blk2,3) = kpk_ec{chani,2}(:,2); % p(k) in block 2
    x(not(Ind_blk2),3) = 0;
    x(Ind_blk3,4) = kpk_ec{chani,3}(:,2); % p(k) in block 3
    x(not(Ind_blk3),4) = 0;
    kpk_ec_tot{chani} = [x(:,1),x(:,2)+x(:,3)+x(:,4)];
    kpk_ec_tot{chani}(:,2) = kpk_ec_tot{chani}(:,2)/3;
    % option: identify the k associated with top p(k), and drop all items that smaller than that k value
    [~,I] = max(kpk_ec_tot{chani}(:,2));
    
    % define linear zone. 05/26/2016
    I1 = 100;
    I2 = 250; % not used. practically: I1:end.
    
    tmp = kpk_ec_tot{chani};
    Ind1 = find(tmp(:,1) == I1);
    Ind2 = find(tmp(:,1) == I2);
    kpk_ec_tot_thre = tmp(Ind1:end,:); 
    % convert to log style
    kpk_log_ec_tot{chani}(:,1) = log2(kpk_ec_tot_thre(:,1));
    kpk_log_ec_tot{chani}(:,2) = log2(kpk_ec_tot_thre(:,2));
    
    % eo
    q = unique([kpk_eo{chani,1}(:,1);kpk_eo{chani,2}(:,1);kpk_eo{chani,3}(:,1)]); % unique sorted catenate.
    Ind_blk1 = ismember(q, kpk_eo{chani,1}(:,1)); % All dataes in block 1
    Ind_blk2 = ismember(q, kpk_eo{chani,2}(:,1)); % All dataes in block 2
    Ind_blk3 = ismember(q, kpk_eo{chani,3}(:,1)); % All dataes in block 3
    x = zeros(length(q),4);
    x(:,1) = q;
    x(Ind_blk1,2) = kpk_eo{chani,1}(:,2);
    x(not(Ind_blk1),2) = 0;
    x(Ind_blk2,3) = kpk_eo{chani,2}(:,2);
    x(not(Ind_blk2),3) = 0;
    x(Ind_blk3,4) = kpk_eo{chani,3}(:,2);
    x(not(Ind_blk3),4) = 0;
    kpk_eo_tot{chani} = [x(:,1),x(:,2)+x(:,3)+x(:,4)];
    kpk_eo_tot{chani}(:,2) = kpk_eo_tot{chani}(:,2)/3;
    % option: identify the k associated with top p(k), and drop all items that smaller than that k value
    [~,I] = max(kpk_eo_tot{chani}(:,2));
    
    % define linear zone. 05/26/2016
    tmp = kpk_eo_tot{chani};
    Ind1 = find(tmp(:,1) == I1);
    Ind2 = find(tmp(:,1) == I2);
    kpk_eo_tot_thre = tmp(Ind1:end,:);
    % convert to log style
    kpk_log_eo_tot{chani}(:,1) = log2(kpk_eo_tot_thre(:,1));
    kpk_log_eo_tot{chani}(:,2) = log2(kpk_eo_tot_thre(:,2));
    
    % rt
    q = unique([kpk_rt{chani,1}(:,1);kpk_rt{chani,2}(:,1);kpk_rt{chani,3}(:,1)]); % unique sorted catenate.
    Ind_blk1 = ismember(q, kpk_rt{chani,1}(:,1)); % All dataes in block 1
    Ind_blk2 = ismember(q, kpk_rt{chani,2}(:,1)); % All dataes in block 2
    Ind_blk3 = ismember(q, kpk_rt{chani,3}(:,1)); % All dataes in block 3
    x = zeros(length(q),4);
    x(:,1) = q;
    x(Ind_blk1,2) = kpk_rt{chani,1}(:,2);
    x(not(Ind_blk1),2) = 0;
    x(Ind_blk2,3) = kpk_rt{chani,2}(:,2);
    x(not(Ind_blk2),3) = 0;
    x(Ind_blk3,4) = kpk_rt{chani,3}(:,2);
    x(not(Ind_blk3),4) = 0;
    kpk_rt_tot{chani} = [x(:,1),x(:,2)+x(:,3)+x(:,4)];
    kpk_rt_tot{chani}(:,2) = kpk_rt_tot{chani}(:,2)/3;
    % option: identify the k associated with top p(k), and drop all items that smaller than that k value
    [~,I] = max(kpk_rt_tot{chani}(:,2));
    
    % define linear zone. 05/26/2016
    tmp = kpk_rt_tot{chani};
    Ind1 = find(tmp(:,1) == I1);
    Ind2 = find(tmp(:,1) == I2);
    kpk_rt_tot_thre = tmp(Ind1:end,:);
    % convert to log style
    kpk_log_rt_tot{chani}(:,1) = log2(kpk_rt_tot_thre(:,1));
    kpk_log_rt_tot{chani}(:,2) = log2(kpk_rt_tot_thre(:,2));
    
    % gng
    q = unique([kpk_mix{chani,1}(:,1);kpk_mix{chani,2}(:,1);kpk_mix{chani,3}(:,1)]); % unique sorted catenate.
    Ind_blk1 = ismember(q, kpk_mix{chani,1}(:,1)); % All dataes in block 1
    Ind_blk2 = ismember(q, kpk_mix{chani,2}(:,1)); % All dataes in block 2
    Ind_blk3 = ismember(q, kpk_mix{chani,3}(:,1)); % All dataes in block 3
    x = zeros(length(q),4);
    x(:,1) = q;
    x(Ind_blk1,2) = kpk_mix{chani,1}(:,2);
    x(not(Ind_blk1),2) = 0;
    x(Ind_blk2,3) = kpk_mix{chani,2}(:,2);
    x(not(Ind_blk2),3) = 0;
    x(Ind_blk3,4) = kpk_mix{chani,3}(:,2);
    x(not(Ind_blk3),4) = 0;
    kpk_mix_tot{chani} = [x(:,1),x(:,2)+x(:,3)+x(:,4)];
    kpk_mix_tot{chani}(:,2) = kpk_mix_tot{chani}(:,2)/3;
    % option: identify the k associated with top p(k), and drop all items that smaller than that k value
    [~,I] = max(kpk_mix_tot{chani}(:,2));
        
    % define linear zone. 05/26/2016    
    tmp = kpk_mix_tot{chani};
    Ind1 = find(tmp(:,1) == I1);
    Ind2 = find(tmp(:,1) == I2);
    kpk_mix_tot_thre = tmp(Ind1:end,:);
    % convert to log style
    kpk_log_mix_tot{chani}(:,1) = log2(kpk_mix_tot_thre(:,1));
    kpk_log_mix_tot{chani}(:,2) = log2(kpk_mix_tot_thre(:,2));
end

kpk_tot     = cat(2,kpk_ec_tot,kpk_eo_tot,kpk_rt_tot,kpk_mix_tot);
kpk_log_tot = cat(2,kpk_log_ec_tot,kpk_log_eo_tot,kpk_log_rt_tot,kpk_log_mix_tot);

% fitting, extract slope for each channel, each condition
% initiate power of scale-freeness (slope) matrix
power_sf = zeros(m,4); % channel X condition. Sequence of condition: ec, eo, rt, irt, mix.

for chani = 1: m
    for conditioni = 1: 4
        tmp = polyfit(kpk_log_tot{chani,conditioni}(:,1), kpk_log_tot{chani,conditioni}(:,2), 1);
        power_sf(chani,conditioni) = -tmp(1);
    end
end


        
        
