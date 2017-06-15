% 05/24/2016 plot PSD for wide band filtered (0.001~6.25 Hz) data
% ---------------- Wide band
hbo = load('C:\Users\Li_Lab537\Documents\Projects\GoNoGo_NIRx_EEG\a2u04\NIRs\2015-07-08_001_EO\exportData\wideBand\NIRS-2015-07-08_001_oxyhb_T1to7503_C1to14.txt');
hbR = load('C:\Users\Li_Lab537\Documents\Projects\GoNoGo_NIRx_EEG\a2u04\NIRs\2015-07-08_001_EO\exportData\wideBand\NIRS-2015-07-08_001_deoxyhb_T1to7503_C1to14.txt');
srate = 12.5;

chan2plot = 14;

figure;
subplot(211); 
[Pxx,F] = periodogram(hbo(:,chan2plot),[],length(hbo(:,chan2plot)),srate);
loglog(F,10*log10(Pxx));xlim([0.001 6.25]);

subplot(212); 
[Pxx,F] = periodogram(hbR(:,chan2plot),[],length(hbR(:,chan2plot)),srate);
loglog(F,10*log10(Pxx));xlim([0.001 6.25]);

% ---------------- Narrow band
hbo = load('C:\Users\Li_Lab537\Documents\Projects\GoNoGo_NIRx_EEG\a2u04\NIRs\2015-07-08_001_EO\exportData\narrowBand\NIRS-2015-07-08_001_oxyhb_T1to7503_C1to14.txt');
hbR = load('C:\Users\Li_Lab537\Documents\Projects\GoNoGo_NIRx_EEG\a2u04\NIRs\2015-07-08_001_EO\exportData\narrowBand\NIRS-2015-07-08_001_deoxyhb_T1to7503_C1to14.txt');
srate = 12.5;

chan2plot = 14;

figure;
subplot(211); 
[Pxx,F] = periodogram(hbo(:,chan2plot),[],length(hbo(:,chan2plot)),srate);
loglog(F,10*log10(Pxx));xlim([0.001 6.25]);

subplot(212); 
[Pxx,F] = periodogram(hbR(:,chan2plot),[],length(hbR(:,chan2plot)),srate);
loglog(F,10*log10(Pxx));xlim([0.001 6.25]);

%% re-construct VG, based on the cleaned data through nirsLAB
subjectArray = {'a2u04','i3h06','i3i03','i4n06','i6a08','l9e11','m1n03','n8n10','u8n09'};

for subji = 1: length(subjectArray)   
    
    subject = cell2mat(subjectArray(subji));
    [~,~] = lz_visibility_cvt_journal_eceo(subject);  
    
end
%%
subjectArray = {'a2u04','i3h06','i3i03','i4n06','i6a08','l9e11','m1n03','n8n10','u8n09'};
for subji = 1: length(subjectArray)   
    
    subject = cell2mat(subjectArray(subji));
    [~] = lz_visibility_cvt_journal_mix(subject);  
    
end
%%
subjectArray = {'a2u04','i3h06','i3i03','i4n06','i6a08','l9e11','m1n03','n8n10','u8n09'};
for subji = 1: length(subjectArray)   
    
    subject = cell2mat(subjectArray(subji));
    [~] = lz_visibility_cvt_journal_rt(subject);  
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  power of scale-freeness and plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ========================================================================
% STEP 4: COMPUTE DEGREE DISTRIBUTION
[~,kpk_tot(:,:,1), power_sf_mat(:,:,1)] = lz_visibility_degreeScale_journal('l9e11');% ,-.1);
[~,kpk_tot(:,:,2), power_sf_mat(:,:,2)] = lz_visibility_degreeScale_journal('i6a08');% ,-.1);
[~,kpk_tot(:,:,3), power_sf_mat(:,:,3)] = lz_visibility_degreeScale_journal('i3i03');% ,-.1);
[~,kpk_tot(:,:,4), power_sf_mat(:,:,4)] = lz_visibility_degreeScale_journal('n8n10');% ,-.1);
[~,kpk_tot(:,:,5), power_sf_mat(:,:,5)] = lz_visibility_degreeScale_journal('a2u04');% ,-.1);
[~,kpk_tot(:,:,6), power_sf_mat(:,:,6)] = lz_visibility_degreeScale_journal('i3h06');% ,-.1);
[~,kpk_tot(:,:,7), power_sf_mat(:,:,7)] = lz_visibility_degreeScale_journal('u8n09');% ,-.1);
[~,kpk_tot(:,:,8), power_sf_mat(:,:,8)] = lz_visibility_degreeScale_journal('i4n06');% ,-.1);
[~,kpk_tot(:,:,9), power_sf_mat(:,:,9)] = lz_visibility_degreeScale_journal('m1n03');% ,-.1);

% in order to present the subject-averaged degree distribution, combine the degree distribution across that from all subjects.
% refer to the correspondingly section in `lz_visibility_degreeScale.m'.
for chani = 1:14
    for condi = 1:4
        
        q = unique([kpk_tot{chani,condi,1}(:,1); kpk_tot{chani,condi,2}(:,1); kpk_tot{chani,condi,3}(:,1); ...
            kpk_tot{chani,condi,4}(:,1); kpk_tot{chani,condi,5}(:,1); kpk_tot{chani,condi,6}(:,1); ...
            kpk_tot{chani,condi,7}(:,1); kpk_tot{chani,condi,8}(:,1); kpk_tot{chani,condi,9}(:,1)]);
        
        x = zeros(length(q), 10);
        x(:,1) = q;
        
        for subi = 1:9
            Ind_blk = ismember(q,kpk_tot{chani,condi,subi}(:,1));
            x(Ind_blk,subi+1) = kpk_tot{chani,condi,subi}(:,2);
            x(not(Ind_blk),subi+1) = 0;
        end
        
        kpk_tot_mnCrossSub{chani,condi} = [ x(:,1), sum(x(:,2:end),2) ];
        kpk_tot_mnCrossSub{chani,condi}(:,2) = kpk_tot_mnCrossSub{chani,condi}(:,2)/9;
        
    end
end

% plot log-log distribution across all subjects
figure(13);clf;
for chani  = 1: 14; 
    subplot(4,4,chani)
    loglog(kpk_tot_mnCrossSub{chani,1}(:,1),(kpk_tot_mnCrossSub{chani,1}(:,2)),'b+')
    hold on; loglog(kpk_tot_mnCrossSub{chani,2}(:,1),(kpk_tot_mnCrossSub{chani,2}(:,2)),'g+')
    hold on; loglog(kpk_tot_mnCrossSub{chani,3}(:,1),(kpk_tot_mnCrossSub{chani,3}(:,2)),'k.')
    hold on; loglog(kpk_tot_mnCrossSub{chani,4}(:,1),(kpk_tot_mnCrossSub{chani,4}(:,2)),'m.')
%     legend('ec','eo','rt','irt','mix')
    title(['[log-log] k vs. p(k), Channel ',num2str(chani), ', Mean Across Subjects'])
    xlabel('k');
    ylabel('p(k)');
%     xlim([0, 1200])
%     ylim([7e-5, 5e-2])
end

% plot zoom-in plot to highlight the ``straight-line'' effects
figure(14);clf;
for chani  = 1: 14; 
    subplot(4,4,chani)
    loglog(kpk_tot_mnCrossSub{chani,1}(:,1),(kpk_tot_mnCrossSub{chani,1}(:,2)),'b+')
    hold on; loglog(kpk_tot_mnCrossSub{chani,2}(:,1),(kpk_tot_mnCrossSub{chani,2}(:,2)),'g+')
    hold on; loglog(kpk_tot_mnCrossSub{chani,3}(:,1),(kpk_tot_mnCrossSub{chani,3}(:,2)),'k.')
    hold on; loglog(kpk_tot_mnCrossSub{chani,4}(:,1),(kpk_tot_mnCrossSub{chani,4}(:,2)),'m.')
    legend('ec','eo','rt','mix')
    title(['[log-log] k vs. p(k), Channel ',num2str(chani), ', Mean Across Subjects'])
    xlabel('k');
    ylabel('p(k)');
    xlim([100, 500])
    ylim([7e-5, 2e-2])
end

% number of channel
m = 14;

% ANOVA for power of scalef-reeness
anova_p_powerScale = zeros(m,1);
group = {'EC','EC','RT','MIX'};
for chani = 1: m
    anova_p_powerScale(chani) = anova1(squeeze(power_sf_mat(chani,:,:))', [], 'off');
end

% average power of scale-freeness and bar plot across subjects
power_sf_mnCrossSub = mean(power_sf_mat,3);
power_sf_stdCrossSub = std(power_sf_mat,0,3);
figure(15);clf
for chani = 1: m
    subplot(2,7,chani)
    barwitherr(power_sf_stdCrossSub(chani,:,:), power_sf_mnCrossSub(chani,:,:))
    set(gca,'XTickLabel',{'EC','EO','RT','MIX'})
    title(['Channel ', num2str(chani)])
    xlabel('Condition')
    ylabel('r')
    ylim([2,5])
    box on; 
    grid on;
    set(gca, 'fontsize',12, 'linew',2)
end
%%
%%%% 06/15/2016 integrate above three figures into one
figure(19); clf;
for chani = 1: 14
    subplot(14,3, (chani-1) * 3 + 1)
    loglog(kpk_tot_mnCrossSub{chani,1}(:,1),(kpk_tot_mnCrossSub{chani,1}(:,2)),'b+')
    hold on; loglog(kpk_tot_mnCrossSub{chani,2}(:,1),(kpk_tot_mnCrossSub{chani,2}(:,2)),'g+','linew',.1)
    hold on; loglog(kpk_tot_mnCrossSub{chani,3}(:,1),(kpk_tot_mnCrossSub{chani,3}(:,2)),'k.','linew',.1)
    hold on; loglog(kpk_tot_mnCrossSub{chani,4}(:,1),(kpk_tot_mnCrossSub{chani,4}(:,2)),'m.','linew',.1)
    grid on; box on;
    if chani ~= 14
        set(gca,'xticklabel',[],'yticklabel',[])
    end
%     legend('ec','eo','rt','irt','mix')
%     title(['[log-log] k vs. p(k), Channel ',num2str(chani), ', Mean Across Subjects'])
    if chani ==14
        xlabel('k');
        ylabel('p(k)');
    end
%     xlabel('k');
%     ylabel('p(k)');
    xlim([0, 1200])
    ylim([7e-5, 5e-2])
    
    subplot(14,3, (chani-1) * 3 + 2)
    loglog(kpk_tot_mnCrossSub{chani,1}(:,1),(kpk_tot_mnCrossSub{chani,1}(:,2)),'b+')
    hold on; loglog(kpk_tot_mnCrossSub{chani,2}(:,1),(kpk_tot_mnCrossSub{chani,2}(:,2)),'g+')
    hold on; loglog(kpk_tot_mnCrossSub{chani,3}(:,1),(kpk_tot_mnCrossSub{chani,3}(:,2)),'k.')
    hold on; loglog(kpk_tot_mnCrossSub{chani,4}(:,1),(kpk_tot_mnCrossSub{chani,4}(:,2)),'m.')
    %     legend('ec','eo','rt','mix')
    %     title(['[log-log] k vs. p(k), Channel ',num2str(chani), ', Mean Across Subjects'])
    %     xlabel('k');
    %     ylabel('p(k)');
    grid on; box on;
    if chani ~= 14
        set(gca,'xticklabel',[],'yticklabel',[])
    end
    if chani ==14
        xlabel('k');
        ylabel('p(k)');
    end
    xlim([100, 500])
    ylim([7e-5, 2e-2])
    
    subplot(14,3, (chani-1) * 3 + 3)
    barwitherr(power_sf_stdCrossSub(chani,:,:), power_sf_mnCrossSub(chani,:,:))
    set(gca,'XTickLabel',{'EC','EO','RT','MIX'})
%     title(['Channel ', num2str(chani)])
    xlim([0,5])
    ylim([2,5])
    box on; 
    grid on;
    if chani ~= 14
        set(gca,'xticklabel',[],'yticklabel',[])
    end
    if chani == 14
        xlabel('Condition')
        ylabel('\gamma','Interpreter','Tex')        
    end
    set(gca, 'fontsize',12, 'linew',1)
end

%% 06/30/2016 re-plot figure(19), but using 2 as base (not use loglog function)
dotSize = 30;

figure(20); clf;
for chani = 1: 14
    subplot(14,3, (chani-1) * 3 + 1)
    scatter(log2(kpk_tot_mnCrossSub{chani,1}(:,1)),log2((kpk_tot_mnCrossSub{chani,1}(:,2))), dotSize, 'b','+')
    hold on; scatter(log2(kpk_tot_mnCrossSub{chani,2}(:,1)),log2((kpk_tot_mnCrossSub{chani,2}(:,2))), dotSize, 'g','+')
    hold on; scatter(log2(kpk_tot_mnCrossSub{chani,3}(:,1)),log2((kpk_tot_mnCrossSub{chani,3}(:,2))), dotSize, 'k','.')
    hold on; scatter(log2(kpk_tot_mnCrossSub{chani,4}(:,1)),log2((kpk_tot_mnCrossSub{chani,4}(:,2))), dotSize, 'm','.')
    ylim([log2(7e-5), log2(5e-2)])
    xlim([log2(0), log2(1200)])

    set(gca, 'XTickLabel',[])                      %# suppress current x-labels
    xt = get(gca, 'XTick');
    yl = get(gca, 'YLim');
    str = cellstr( num2str(xt(:),'2^{%d}') );      %# format x-ticks as 2^{xx}
    hTxt = text(xt, yl(ones(size(xt))), str, ...   %# create text at same locations
        'Interpreter','tex', ...                   %# specify tex interpreter
        'VerticalAlignment','top', ...             %# v-align to be underneath
        'HorizontalAlignment','center');           %# h-aligh to be centered

    set(gca, 'YTickLabel',[])                      %# suppress current x-labels
    yt = get(gca, 'YTick');
    xl = get(gca, 'XLim');
    str = cellstr( num2str(yt(:),'2^{%d}') );      %# format x-ticks as 2^{xx}
    hTxt = text(yt, xl(ones(size(yt))), str, ...   %# create text at same locations
        'Interpreter','tex', ...                   %# specify tex interpreter
        'VerticalAlignment','middle', ...             %# v-align to be underneath
        'HorizontalAlignment','center');           %# h-aligh to be centered
    
    box on;
    if chani ~= 14
        set(gca,'xticklabel',[],'yticklabel',[])
    end
%     legend('ec','eo','rt','irt','mix')
%     title(['[log-log] k vs. p(k), Channel ',num2str(chani), ', Mean Across Subjects'])
    if chani ==14
        xlabel('k');
        ylabel('p(k)');
    end
%     xlabel('k');
%     ylabel('p(k)');

legend('EC', 'EO', 'RT', 'GNG')
    
    subplot(14,3, (chani-1) * 3 + 2)
    scatter(log2(kpk_tot_mnCrossSub{chani,1}(:,1)),(log2(kpk_tot_mnCrossSub{chani,1}(:,2))),dotSize, 'b','+')
    hold on; scatter(log2(kpk_tot_mnCrossSub{chani,2}(:,1)),(log2(kpk_tot_mnCrossSub{chani,2}(:,2))),dotSize, 'g','+')
    hold on; scatter(log2(kpk_tot_mnCrossSub{chani,3}(:,1)),(log2(kpk_tot_mnCrossSub{chani,3}(:,2))),dotSize, 'k','.')
    hold on; scatter(log2(kpk_tot_mnCrossSub{chani,4}(:,1)),(log2(kpk_tot_mnCrossSub{chani,4}(:,2))),dotSize, 'm','.')
    xlim([log2(100), log2(500)])
    ylim([log2(7e-5), log2(2e-2)])    %     legend('ec','eo','rt','mix')
    set(gca, 'XTickLabel',[])                      %# suppress current x-labels
    xt = get(gca, 'XTick');
    yl = get(gca, 'YLim');
    str = cellstr( num2str(xt(:),'2^{%d}') );      %# format x-ticks as 2^{xx}
    hTxt = text(xt, yl(ones(size(xt))), str, ...   %# create text at same locations
        'Interpreter','tex', ...                   %# specify tex interpreter
        'VerticalAlignment','top', ...             %# v-align to be underneath
        'HorizontalAlignment','center');           %# h-aligh to be centered
    
%     set(gca, 'YTickLabel',[])                      %# suppress current x-labels
    yt = get(gca, 'YTick');
    xl = get(gca, 'XLim');
    str = cellstr( num2str(yt(:),'2^{%d}') );      %# format x-ticks as 2^{xx}
    hTxt = text(yt, xl(ones(size(yt))), str, ...   %# create text at same locations
        'Interpreter','tex', ...                   %# specify tex interpreter
        'VerticalAlignment','middle', ...             %# v-align to be underneath
        'HorizontalAlignment','left');           %# h-aligh to be centered
    
    box on;
    if chani ~= 14
        set(gca,'xticklabel',[],'yticklabel',[])
    end
    if chani ==14
        xlabel('k');
        ylabel('p(k)');
    end

    
    subplot(14,3, (chani-1) * 3 + 3)
    barwitherr(power_sf_stdCrossSub(chani,:,:), power_sf_mnCrossSub(chani,:,:))
    set(gca,'XTickLabel',{'EC','EO','RT','MIX'})
%     title(['Channel ', num2str(chani)])
    xlim([0,5])
    ylim([2,5])
    box on; 
    grid on;
    if chani ~= 14
        set(gca,'xticklabel',[],'yticklabel',[])
    end
    if chani == 14
        xlabel('Condition')
        ylabel('\gamma','Interpreter','Tex')        
    end
    set(gca, 'fontsize',12, 'linew',1)
end
%%
% t-test for power of scale-freeness, between rest and task
t_h_powerScale = zeros(m,1);
t_p_powerScale = zeros(m,1);
group = {'EC','EO','RT','MIX'};
for chani = 1: m
   [t_h_powerScale(chani), t_p_powerScale(chani)] = ttest(reshape(cat(2,power_sf_mat(chani,1:2,:)),1,18), reshape(cat(2,power_sf_mat(chani,3:4,:)),1,18),'Tail','right');
end
%%%% t-test EC vs EO
for chani = 1: m
   [t_h_powerScale_ECvsEO(chani), t_p_powerScale_ECvsEO(chani)] = ttest(power_sf_mat(chani,1,:), power_sf_mat(chani,2,:),'Tail','right');
end
%%%% t-test EC vs RT
for chani = 1: m
   [t_h_powerScale_ECvsRT(chani), t_p_powerScale_ECvsRT(chani)] = ttest(power_sf_mat(chani,1,:), power_sf_mat(chani,3,:),'Tail','right');
end
%%%% t-test EC vs MIX
for chani = 1: m
   [t_h_powerScale_ECvsMIX(chani), t_p_powerScale_ECvsMIX(chani)] = ttest(power_sf_mat(chani,1,:), power_sf_mat(chani,4,:),'Tail','right');
end
%%%% t-test EO vs RT
for chani = 1: m
   [t_h_powerScale_EOvsRT(chani), t_p_powerScale_EOvsRT(chani)] = ttest(power_sf_mat(chani,2,:), power_sf_mat(chani,3,:),'Tail','right');
end
%%%% t-test EO vs MIT
for chani = 1: m
   [t_h_powerScale_EOvsMIT(chani), t_p_powerScale_EOvsMIT(chani)] = ttest(power_sf_mat(chani,2,:), power_sf_mat(chani,4,:),'Tail','right');
end
%%%% t-test RT vs MIT
for chani = 1: m
   [t_h_powerScale_RTvsMIT(chani), t_p_powerScale_RTvsMIT(chani)] = ttest(power_sf_mat(chani,3,:), power_sf_mat(chani,4,:),'Tail','right');
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%% COMPUTE MEAN PATHLENGTH %%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pathlength_l9e11_mn = lz_visibility_pathlength_journal('l9e11');
pathlength_i6a08_mn = lz_visibility_pathlength_journal('i6a08');
pathlength_i3i03_mn = lz_visibility_pathlength_journal('i3i03');
pathlength_n8n10_mn = lz_visibility_pathlength_journal('n8n10');
pathlength_a2u04_mn = lz_visibility_pathlength_journal('a2u04');
pathlength_i3h06_mn = lz_visibility_pathlength_journal('i3h06');
pathlength_u8n09_mn = lz_visibility_pathlength_journal('u8n09');
pathlength_i4n06_mn = lz_visibility_pathlength_journal('i4n06');
pathlength_m1n03_mn = lz_visibility_pathlength_journal('m1n03');

pathlength_mn = cat(3,pathlength_l9e11_mn,pathlength_i6a08_mn,pathlength_i3i03_mn,...
    pathlength_n8n10_mn, pathlength_a2u04_mn, pathlength_i3h06_mn, pathlength_u8n09_mn,...
    pathlength_i4n06_mn,pathlength_m1n03_mn);

% ANOVA for mean path length
anova_p_pathlength = zeros(m,1);
group = {'EC','EO','RT','MIX'};
for chani = 1:m
    anova_p_pathlength(chani) = anova1(squeeze(pathlength_mn(chani,:,:))', [], 'off');
end

% average pathlength across subjects and bar plot across subjects
pathlength_mnCrossSub = mean(pathlength_mn,3);
pathlength_stdCrossSub = std(pathlength_mn,0,3);
figure(6);clf
for chani = 1: m
    subplot(4,4,chani)
    h = barwitherr(pathlength_stdCrossSub(chani,:,:), pathlength_mnCrossSub(chani,:,:));
    set(gca,'XTickLabel',{'EC','EO','RT','MIX'})
    title(['Average Pathlength, Channel ', num2str(chani)])
    xlabel('Condition')
    ylabel('L')
    ylim([0.05, 0.1])
    grid on;
end

% t-test for power of scale-freeness, between rest and task
t_h_pathlength = zeros(m,1);
t_p_pathlength = zeros(m,1);
group = {'EC','EO','RT','MIX'};
for chani = 1: m
   [t_h_pathlength(chani), t_p_pathlength(chani)] = ttest(squeeze(mean(pathlength_mn(chani,1:2,:),2)), squeeze(mean(pathlength_mn(chani,3:4,:),2)),'Tail','right');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% compute fractal dimension for each channel, condition, block, and subject.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subArray = {'m1n03', 'i6a08','i3i03','n8n10','a2u04','i3h06','i4n06','u8n09','l9e11'};
% initiate the D and H matrix
D_rd = zeros(14, 3, 4, 9); % subject X channel X block X condition

H_rd = zeros(14, 3, 4, 9); % subject X channel X block X condition

for subi = 1: 9
    sub2comp = cell2mat(subArray(subi));
    
    % using relative dispersion
    [D_rd(:,:,:,subi), H_rd(:,:,:,subi)] = lz_fractality_journal_swv(sub2comp, 1, 0);
        % subject X Channel X block X condition
%     % using generalized Hurst Exponont approach
%     H_genhurst(:,:,:,subi) = lz_fractality_genhurst_journal(sub2comp, 1, 0);
 
end

%%% plot D values in hist across all time series
figure(20);clf;
subplot(121)
hist(D_rd(:))
set(gca,'fontsize',22,'fontweight','bold')
xlabel('Fractal Dimension')
ylabel('Frequency')
box on;
subplot(122)
boxplot(D_rd(:));
set(gca,'fontsize',30,'fontweight','bold','linew',2)

%%% plot least square for an examplary time series in scatter plot
chan2plot = 1; %% channel ID to be plotted
blk2plot  = 3; %% block ID to be plotted
% subject ID is l9e11

minPow2plot = 2;
maxPow2plot = 7;

figure(11);clf; 
plot(minPow2plot:maxPow2plot, log2(rd_eo_orig(minPow2plot:maxPow2plot,chan2plot,blk2plot)),'ro','linew',4)
    %%%% "rd_eo_orig" can be computed from "everyday_biocas2016.m".
% hold on
% plot(minPow2plot:maxPow2plot, log2(rd_eo_shuf(minPow2plot:maxPow2plot,chan2plot,1)),'bv','linew',1)
% legend1 = legend('Unshuffled Signal', 'Shuffled Signal');
% set(legend1, 'Location', 'southeast')
set(gca,'FontSize',24,'FontWeight','bold','linewidth',2)
ls = lsline;                    
set(ls(1),'color','r','linew',3)
% set(ls(2),'color','r','linew',3)
grid on
box on
ylabel('log_2 (\tilde{\sigma})');
xlabel('log_2(n)');

% plot raw signals and 

%%
% plot result of Hurst to compare across conditions
figure(21); clf;
for chani = 1:14
    H2plot = mean(squeeze(mean(squeeze(H_rd(chani,:,:,:)),1)),2); % average across subjects and blocks
    idx2plot = repelem([1 2 3 4],1);
    subplot(4,4,chani)
    scatter(idx2plot,H2plot,'ro','linew',4);
    ax = gca;
    ax.XTickLabel = {'','EC','EO','RT', 'GNG'};
    ls = lsline;
    set(ls,'color','k','linew',4)
    title(['Channel ', num2str(chani)])
    grid on;
    box on;
    set(gca,'FontSize',26,'FontWeight','bold','linewidth',2)
end

%%
% plot scatter between power_sf and Hurst (different method)
H_rd_mnAcross_blk = squeeze(mean(H_rd, 2));
% H_genhurst_mnAcross_blk = squeeze(mean(H_genhurst, 2));

H_rd2Comp = reshape(H_rd_mnAcross_blk,14*4*9,1);
power_sf2Comp = reshape(power_sf_mat,14*4*9,1);
% H_genhurst2Comp = reshape(H_genhurst_mnAcross_blk,14*4*9,1);

fig22 = figure(22);clf
scatter(H_rd2Comp,power_sf2Comp,'dia');
xlabel('Hurst Exponent')
ylabel('Power of Scale-freeness')
hold on
my_poly=polyfit(H_rd2Comp,power_sf2Comp,1);
% X2= 0.88:0.01:1; % X data range 
X2=0.75:0.01:0.95
Y2=polyval(my_poly,X2);

plot(X2,Y2,'linew',3);
box on;
grid on;
set(gca,'fontsize',28,'linew',5,'fontweight','bold')

H_rd2Comp(204) = NaN; % exclude outlier (the min of H_rd2Comp).

[R,P] = corrcoef(H_rd2Comp,power_sf2Comp,'rows','complete')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% use logfit.m function to re-extract power of scale-freeness for each channel, block, and subject
subject = 'm1n03';

idx = find(kpk_tot{2,4}(:,1)==57);
[kpk_tot] = lz_visibility_degreeScale_journal_logfit(subject,'4p16') % format of kpk_tot: cell array, channel X condition
figure;
subplot(211)
[slope, intercept] = logfit(kpk_tot{2,4}(:,1),kpk_tot{2,4}(:,2),'loglog','skip',idx)
title('mix')
subplot(212)
[slope, intercept] = logfit(kpk_tot{2,2}(:,1),kpk_tot{2,2}(:,2),'loglog','skip',idx)
title('eo')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% use logfit.m function to compare the degree distribution between pre-filtered and post-filtered time series
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subject = 'm1n03';

%%%% 1. no filtered:
ec_nofilt = load('C:\Users\Li_Lab537\Documents\Projects\GoNoGo_NIRx_EEG\m1n03\NIRS\2015-07-14_002_EC\exportData\wideBand\NIRS-2015-07-14_002_oxyhb_T1to7503_C1to14.txt');
tmp = ec_nofilt; 
clear ec_nofilt;
ec_nofilt.blk1 = tmp(564 : 564+2124,:);
ec_nofilt.blk2 = tmp(564+2124 : 564+2124+2124,:);
ec_nofilt.blk3 = tmp(564+2124+2124 : 564+2124+2124+2124,:);
ec_nofilt = cat(3,ec_nofilt.blk1,ec_nofilt.blk2,ec_nofilt.blk3);
clear tmp;
% convert time-seris to VG
n = size(ec_nofilt,1); % length of time-series
m = size(ec_nofilt,2); % # of channel
VG_ec_nofilt  = zeros(n,n,m,3);
parfor blocki = 1:3
    VG_ec_nofilt(:,:,:,blocki)  = lz_VG_build(ec_nofilt(:,:,blocki));
end
clear blocki 

%%%% 2. band-filtered:
ec = load('C:\Users\Li_Lab537\Documents\Projects\GoNoGo_NIRx_EEG\m1n03\NIRS\2015-07-14_002_EC\exportData\wideBand\NIRS-2015-07-14_002_oxyhb_T1to7503_C1to14.txt');
tmp = ec; 
clear ec;
ec.blk1 = tmp(564 : 564+2124,:);
ec.blk2 = tmp(564+2124 : 564+2124+2124,:);
ec.blk3 = tmp(564+2124+2124 : 564+2124+2124+2124,:);
ec = cat(3,ec.blk1,ec.blk2,ec.blk3);
clear tmp;
% filter each signal with [0.01-0.2] Hz if needed
srate = 12.5;
order = 4;
lobound = 0.01;
hibound = 0.2;
type = 'bandpass';
for blocki = 1:3
    ec_bandFilt(:,:,blocki)  = lz_butterworth(ec(:,:,blocki),srate,order,lobound,hibound,type);
end
clear srate order lobound hibound type blocki eo ec rt irt mix;
% convert time-seris to VG
n = size(ec_bandFilt,1); % length of time-series
m = size(ec_bandFilt,2); % # of channel
VG_ec_bandFilt  = zeros(n,n,m,3);
parfor blocki = 1:3
    VG_ec_bandFilt(:,:,:,blocki)  = lz_VG_build(ec_bandFilt(:,:,blocki));
end
clear blocki 

%%%% 3. low-pass filtered:
ec = load('C:\Users\Li_Lab537\Documents\Projects\GoNoGo_NIRx_EEG\m1n03\NIRS\2015-07-14_002_EC\exportData\wideBand\NIRS-2015-07-14_002_oxyhb_T1to7503_C1to14.txt');
tmp = ec; 
clear ec;
ec.blk1 = tmp(564 : 564+2124,:);
ec.blk2 = tmp(564+2124 : 564+2124+2124,:);
ec.blk3 = tmp(564+2124+2124 : 564+2124+2124+2124,:);
ec = cat(3,ec.blk1,ec.blk2,ec.blk3);
clear tmp;
% low-pass filtered each signal with [0.2] Hz if needed
srate = 12.5;
order = 4;
lobound = 0.01;
hibound = 0.2;
type = 'low';
for blocki = 1:3
    ec_lowFilt(:,:,blocki)  = lz_butterworth(ec(:,:,blocki),srate,order,lobound,hibound,type);
end
clear srate order lobound hibound type blocki eo rt irt mix;
% convert time-seris to VG
n = size(ec,1); % length of time-series
m = size(ec,2); % # of channel
VG_ec_lowFilt  = zeros(n,n,m,3);
parfor blocki = 1:3
    VG_ec_lowFilt(:,:,:,blocki)  = lz_VG_build(ec_lowFilt(:,:,blocki));
end
clear blocki 

%%%% 4. high-pass filtered:
ec = load('C:\Users\Li_Lab537\Documents\Projects\GoNoGo_NIRx_EEG\m1n03\NIRS\2015-07-14_002_EC\exportData\wideBand\NIRS-2015-07-14_002_oxyhb_T1to7503_C1to14.txt');
tmp = ec; 
clear ec;
ec.blk1 = tmp(564 : 564+2124,:);
ec.blk2 = tmp(564+2124 : 564+2124+2124,:);
ec.blk3 = tmp(564+2124+2124 : 564+2124+2124+2124,:);
ec = cat(3,ec.blk1,ec.blk2,ec.blk3);
clear tmp;
% high-pass filtered each signal with [0.01] Hz if needed
srate = 12.5;
order = 4;
lobound = 0.01;
hibound = 0.2;
type = 'high';
for blocki = 1:3
    ec_highFilt(:,:,blocki)  = lz_butterworth(ec(:,:,blocki),srate,order,lobound,hibound,type);
end
clear srate order lobound hibound type blocki eo rt irt mix;
% convert time-seris to VG
n = size(ec,1); % length of time-series
m = size(ec,2); % # of channel
VG_ec_highFilt  = zeros(n,n,m,3);
parfor blocki = 1:3
    VG_ec_highFilt(:,:,:,blocki)  = lz_VG_build(ec_highFilt(:,:,blocki));
end
clear blocki 

%%%% 5. no filtered but shuffle:
for chani = 1:m
    for blocki = 1:3
        ec_shuffle(:,chani,blocki) = shuffle(ec_nofilt(:,chani,blocki),1);
    end
end
% convert time-seris to VG
n = size(ec_shuffle,1); % length of time-series
m = size(ec_shuffle,2); % # of channel
VG_ec_shuffle  = zeros(n,n,m,3);
parfor blocki = 1:3
    VG_ec_shuffle(:,:,:,blocki)  = lz_VG_build(ec_shuffle(:,:,blocki));
end
clear blocki 
%%
%%%% compute degree distribution
%%%% 1. nofiltered.
m = 14; % # of channels
n = size(VG_ec_nofilt,1); % # of nodes
% initiate power of scale-freeness
degree_powerScale_ec_nofilt  = zeros(m,3); % channel X block
% compute degree for each node
% initiate degree
degree_ec  = zeros(n,m,3);
for chani = 1: m % loop for channel
    for blocki =1: 3 % loop for block
        degree_ec_nofilt(:,chani,blocki)  = squeeze(sum(VG_ec_nofilt(:,:,chani,blocki)));
    end
end
% compute p(k) for each condition
% initiate k and p(k) matrix for each condition
kpk_ec_nofilt  = cell(m,3);
for chani = 1: m
    for blocki = 1: 3
        seq_bins = unique(degree_ec_nofilt(:,chani,blocki));
        histcount = histc(degree_ec_nofilt(:,chani,blocki), seq_bins);
        kpk_ec_nofilt(chani,blocki) = {[seq_bins, histcount/n]};
    end
end
% combine degree values across blocks for each channel and condition
% the following is the code developed based on http://www.mathworks.com/matlabcentral/answers/125787-merging-two-matrices-by-first-column-values
% initiate k vs. p(k) cell array
kpk_ec_tot_nofilt  = cell(m,1);
for chani = 1:m
    % ec
    q = unique([kpk_ec_nofilt{chani,1}(:,1);kpk_ec_nofilt{chani,2}(:,1);kpk_ec_nofilt{chani,3}(:,1)]); % unique sorted catenate.
    Ind_blk1 = ismember(q, kpk_ec_nofilt{chani,1}(:,1)); % All dataes in block 1
    Ind_blk2 = ismember(q, kpk_ec_nofilt{chani,2}(:,1)); % All dataes in block 2
    Ind_blk3 = ismember(q, kpk_ec_nofilt{chani,3}(:,1)); % All dataes in block 3
    x = zeros(length(q),4);
    x(:,1) = q; % value of k
    x(Ind_blk1,2) = kpk_ec_nofilt{chani,1}(:,2); % p(k) in block 1
    x(not(Ind_blk1),2) = 0;
    x(Ind_blk2,3) = kpk_ec_nofilt{chani,2}(:,2); % p(k) in block 2
    x(not(Ind_blk2),3) = 0;
    x(Ind_blk3,4) = kpk_ec_nofilt{chani,3}(:,2); % p(k) in block 3
    x(not(Ind_blk3),4) = 0;
    kpk_ec_tot_nofilt{chani} = [x(:,1),x(:,2)+x(:,3)+x(:,4)];
    kpk_ec_tot_nofilt{chani}(:,2) = kpk_ec_tot_nofilt{chani}(:,2)/3;
end

%%%% 2. band_filtered.
m = 14; % # of channels
n = size(VG_ec_bandFilt,1); % # of nodes
% initiate power of scale-freeness
degree_powerScale_ec_bandFilt  = zeros(m,3); % channel X block
% compute degree for each node
% initiate degree
degree_ec  = zeros(n,m,3);
for chani = 1: m % loop for channel
    for blocki =1: 3 % loop for block
        degree_ec_bandFilt(:,chani,blocki)  = squeeze(sum(VG_ec_bandFilt(:,:,chani,blocki)));
    end
end
% compute p(k) for each condition
% initiate k and p(k) matrix for each condition
kpk_ec_bandFilt  = cell(m,3);
for chani = 1: m
    for blocki = 1: 3
        seq_bins = unique(degree_ec_bandFilt(:,chani,blocki));
        histcount = histc(degree_ec_bandFilt(:,chani,blocki), seq_bins);
        kpk_ec_bandFilt(chani,blocki) = {[seq_bins, histcount/n]};
    end
end
% combine degree values across blocks for each channel and condition
% the following is the code developed based on http://www.mathworks.com/matlabcentral/answers/125787-merging-two-matrices-by-first-column-values
% initiate k vs. p(k) cell array
kpk_ec_tot_bandFilt  = cell(m,1);
for chani = 1:m
    % ec
    q = unique([kpk_ec_bandFilt{chani,1}(:,1);kpk_ec_bandFilt{chani,2}(:,1);kpk_ec_bandFilt{chani,3}(:,1)]); % unique sorted catenate.
    Ind_blk1 = ismember(q, kpk_ec_bandFilt{chani,1}(:,1)); % All dataes in block 1
    Ind_blk2 = ismember(q, kpk_ec_bandFilt{chani,2}(:,1)); % All dataes in block 2
    Ind_blk3 = ismember(q, kpk_ec_bandFilt{chani,3}(:,1)); % All dataes in block 3
    x = zeros(length(q),4);
    x(:,1) = q; % value of k
    x(Ind_blk1,2) = kpk_ec_bandFilt{chani,1}(:,2); % p(k) in block 1
    x(not(Ind_blk1),2) = 0;
    x(Ind_blk2,3) = kpk_ec_bandFilt{chani,2}(:,2); % p(k) in block 2
    x(not(Ind_blk2),3) = 0;
    x(Ind_blk3,4) = kpk_ec_bandFilt{chani,3}(:,2); % p(k) in block 3
    x(not(Ind_blk3),4) = 0;
    kpk_ec_tot_bandFilt{chani} = [x(:,1),x(:,2)+x(:,3)+x(:,4)];
    kpk_ec_tot_bandFilt{chani}(:,2) = kpk_ec_tot_bandFilt{chani}(:,2)/3;
end

%%%% 3. low-pass filtered.
m = 14; % # of channels
n = size(VG_ec_lowFilt,1); % # of nodes
% initiate power of scale-freeness
degree_powerScale_ec_lowFilt  = zeros(m,3); % channel X block
% compute degree for each node
% initiate degree
degree_ec  = zeros(n,m,3);
for chani = 1: m % loop for channel
    for blocki =1: 3 % loop for block
        degree_ec_lowFilt(:,chani,blocki)  = squeeze(sum(VG_ec_lowFilt(:,:,chani,blocki)));
    end
end
% compute p(k) for each condition
% initiate k and p(k) matrix for each condition
kpk_ec_lowFilt  = cell(m,3);
for chani = 1: m
    for blocki = 1: 3
        seq_bins = unique(degree_ec_lowFilt(:,chani,blocki));
        histcount = histc(degree_ec_lowFilt(:,chani,blocki), seq_bins);
        kpk_ec_lowFilt(chani,blocki) = {[seq_bins, histcount/n]};
    end
end
% combine degree values across blocks for each channel and condition
% the following is the code developed based on http://www.mathworks.com/matlabcentral/answers/125787-merging-two-matrices-by-first-column-values
% initiate k vs. p(k) cell array
kpk_ec_tot_lowFilt  = cell(m,1);
for chani = 1:m
    % ec
    q = unique([kpk_ec_lowFilt{chani,1}(:,1);kpk_ec_lowFilt{chani,2}(:,1);kpk_ec_lowFilt{chani,3}(:,1)]); % unique sorted catenate.
    Ind_blk1 = ismember(q, kpk_ec_lowFilt{chani,1}(:,1)); % All dataes in block 1
    Ind_blk2 = ismember(q, kpk_ec_lowFilt{chani,2}(:,1)); % All dataes in block 2
    Ind_blk3 = ismember(q, kpk_ec_lowFilt{chani,3}(:,1)); % All dataes in block 3
    x = zeros(length(q),4);
    x(:,1) = q; % value of k
    x(Ind_blk1,2) = kpk_ec_lowFilt{chani,1}(:,2); % p(k) in block 1
    x(not(Ind_blk1),2) = 0;
    x(Ind_blk2,3) = kpk_ec_lowFilt{chani,2}(:,2); % p(k) in block 2
    x(not(Ind_blk2),3) = 0;
    x(Ind_blk3,4) = kpk_ec_lowFilt{chani,3}(:,2); % p(k) in block 3
    x(not(Ind_blk3),4) = 0;
    kpk_ec_tot_lowFilt{chani} = [x(:,1),x(:,2)+x(:,3)+x(:,4)];
    kpk_ec_tot_lowFilt{chani}(:,2) = kpk_ec_tot_lowFilt{chani}(:,2)/3;
end


%%%% 4. high-pass filtered.
m = 14; % # of channels
n = size(VG_ec_highFilt,1); % # of nodes
% initiate power of scale-freeness
degree_powerScale_ec_highFilt  = zeros(m,3); % channel X block
% compute degree for each node
% initiate degree
degree_ec  = zeros(n,m,3);
for chani = 1: m % loop for channel
    for blocki =1: 3 % loop for block
        degree_ec_highFilt(:,chani,blocki)  = squeeze(sum(VG_ec_highFilt(:,:,chani,blocki)));
    end
end
% compute p(k) for each condition
% initiate k and p(k) matrix for each condition
kpk_ec_highFilt  = cell(m,3);
for chani = 1: m
    for blocki = 1: 3
        seq_bins = unique(degree_ec_highFilt(:,chani,blocki));
        histcount = histc(degree_ec_highFilt(:,chani,blocki), seq_bins);
        kpk_ec_highFilt(chani,blocki) = {[seq_bins, histcount/n]};
    end
end
% combine degree values across blocks for each channel and condition
% the following is the code developed based on http://www.mathworks.com/matlabcentral/answers/125787-merging-two-matrices-by-first-column-values
% initiate k vs. p(k) cell array
kpk_ec_tot_highFilt  = cell(m,1);
for chani = 1:m
    % ec
    q = unique([kpk_ec_highFilt{chani,1}(:,1);kpk_ec_highFilt{chani,2}(:,1);kpk_ec_highFilt{chani,3}(:,1)]); % unique sorted catenate.
    Ind_blk1 = ismember(q, kpk_ec_highFilt{chani,1}(:,1)); % All dataes in block 1
    Ind_blk2 = ismember(q, kpk_ec_highFilt{chani,2}(:,1)); % All dataes in block 2
    Ind_blk3 = ismember(q, kpk_ec_highFilt{chani,3}(:,1)); % All dataes in block 3
    x = zeros(length(q),4);
    x(:,1) = q; % value of k
    x(Ind_blk1,2) = kpk_ec_highFilt{chani,1}(:,2); % p(k) in block 1
    x(not(Ind_blk1),2) = 0;
    x(Ind_blk2,3) = kpk_ec_highFilt{chani,2}(:,2); % p(k) in block 2
    x(not(Ind_blk2),3) = 0;
    x(Ind_blk3,4) = kpk_ec_highFilt{chani,3}(:,2); % p(k) in block 3
    x(not(Ind_blk3),4) = 0;
    kpk_ec_tot_highFilt{chani} = [x(:,1),x(:,2)+x(:,3)+x(:,4)];
    kpk_ec_tot_highFilt{chani}(:,2) = kpk_ec_tot_highFilt{chani}(:,2)/3;
end

%%%% 5. no-filtered but shuffled.
m = 14; % # of channels
n = size(VG_ec_shuffle,1); % # of nodes
% initiate power of scale-freeness
degree_powerScale_ec_shuffle  = zeros(m,3); % channel X block
% compute degree for each node
% initiate degree
degree_ec_shuffle  = zeros(n,m,3);
for chani = 1: m % loop for channel
    for blocki =1: 3 % loop for block
        degree_ec_shuffle(:,chani,blocki)  = squeeze(sum(VG_ec_shuffle(:,:,chani,blocki)));
    end
end
% compute p(k) for each condition
% initiate k and p(k) matrix for each condition
kpk_ec_shuffle  = cell(m,3);
for chani = 1: m
    for blocki = 1: 3
        seq_bins = unique(degree_ec_shuffle(:,chani,blocki));
        histcount = histc(degree_ec_shuffle(:,chani,blocki), seq_bins);
        kpk_ec_shuffle(chani,blocki) = {[seq_bins, histcount/n]};
    end
end
% combine degree values across blocks for each channel and condition
% the following is the code developed based on http://www.mathworks.com/matlabcentral/answers/125787-merging-two-matrices-by-first-column-values
% initiate k vs. p(k) cell array
kpk_ec_tot_shuffle  = cell(m,1);
for chani = 1:m
    % ec
    q = unique([kpk_ec_shuffle{chani,1}(:,1);kpk_ec_shuffle{chani,2}(:,1);kpk_ec_shuffle{chani,3}(:,1)]); % unique sorted catenate.
    Ind_blk1 = ismember(q, kpk_ec_shuffle{chani,1}(:,1)); % All dataes in block 1
    Ind_blk2 = ismember(q, kpk_ec_shuffle{chani,2}(:,1)); % All dataes in block 2
    Ind_blk3 = ismember(q, kpk_ec_shuffle{chani,3}(:,1)); % All dataes in block 3
    x = zeros(length(q),4);
    x(:,1) = q; % value of k
    x(Ind_blk1,2) = kpk_ec_shuffle{chani,1}(:,2); % p(k) in block 1
    x(not(Ind_blk1),2) = 0;
    x(Ind_blk2,3) = kpk_ec_shuffle{chani,2}(:,2); % p(k) in block 2
    x(not(Ind_blk2),3) = 0;
    x(Ind_blk3,4) = kpk_ec_shuffle{chani,3}(:,2); % p(k) in block 3
    x(not(Ind_blk3),4) = 0;
    kpk_ec_tot_shuffle{chani} = [x(:,1),x(:,2)+x(:,3)+x(:,4)];
    kpk_ec_tot_shuffle{chani}(:,2) = kpk_ec_tot_shuffle{chani}(:,2)/3;
end
%% plot and fit by logfit.m
chan2plot = 2;
srate = 12.5;

figure(2); clf

subplot(521)
k_min = 37; 
k_max = 115; 
idx1 = find(kpk_ec_tot_nofilt{chan2plot}(:,1) == k_min);
idx2 = find(kpk_ec_tot_nofilt{chan2plot}(:,1) == k_max);
xlab  = linspace(1,ntime/srate,ntime);
plot(xlab,ec_nofilt(:,chan2plot,1),'linew',1);
xlabel('time (s)')
ylab = ylabel('Amplitue (\mu M)');
set(ylab,'Interpreter','tex');
title('ec\_non-filter')
box on;grid on;
set(gca, 'fontw','bold','fontsize',16,'xlim',[1,ntime/srate],'linew',2)
subplot(522)

[slope_nofilt, intercept] = logfit(kpk_ec_tot_nofilt{chan2plot}(:,1),kpk_ec_tot_nofilt{chan2plot}(:,2),'loglog','linewidth',3,'skipBegin',idx1,'skipEnd',length(kpk_ec_tot_nofilt{chan2plot}(:,1))-idx2)
title('ec\_non-filter')
box on;grid on;
set(gca, 'fontw','bold','fontsize',16,'xlim',[1e0 1e3],'linew',2)
text2print = ['PS = ',num2str(roundn(-slope_nofilt,-3))];
text(180,0.01,text2print,'color','r', 'HorizontalAlignment','left')
xlabel('log_2(k)'); ylabel('log_2[p(k)]');

    

subplot(523)
k_min = 130; 
k_max = 350; 
idx1 = find(kpk_ec_tot_lowFilt{chan2plot}(:,1) == k_min);
idx2 = find(kpk_ec_tot_lowFilt{chan2plot}(:,1) == k_max);
xlab  = linspace(1,ntime/srate,ntime);
plot(xlab,ec_lowFilt(:,chan2plot,1),'linew',1);
xlabel('time (s)')
ylab = ylabel('Amplitue (\mu M)');
set(ylab,'Interpreter','tex');
title('ec\_low-pass filtered')
box on;grid on;
set(gca, 'fontw','bold','fontsize',16,'xlim',[1,ntime/srate],'linew',2)
subplot(524)
[slope_lowFilt, intercept] = logfit(kpk_ec_tot_lowFilt{chan2plot}(:,1),kpk_ec_tot_lowFilt{chan2plot}(:,2),'loglog','linewidth',3,'skipBegin',idx1,'skipEnd',length(kpk_ec_tot_lowFilt{chan2plot}(:,1))-idx2)
title('ec\_low-pass filtered')
box on;grid on;
set(gca, 'fontw','bold','fontsize',16,'xlim',[1e0 1e3],'linew',2)
% text2print = ['PS = ',num2str(roundn(-slope_lowFilt,-3))];
% text(180,0.01,text2print,'color','r', 'HorizontalAlignment','left')
xlabel('log_2(k)'); ylabel('log_2[p(k)]');

subplot(525)
k_min = 37; 
k_max = 130; 
idx1 = find(kpk_ec_tot_highFilt{chan2plot}(:,1)==k_min);
idx2 = find(kpk_ec_tot_highFilt{chan2plot}(:,1) == k_max);
xlab  = linspace(1,ntime/srate,ntime);
plot(xlab,ec_highFilt(:,chan2plot,1),'linew',1);
xlabel('time (s)')
ylab = ylabel('Amplitue (\mu M)');
set(ylab,'Interpreter','tex');
title('ec\_high-pass filtered')
box on;grid on;
set(gca, 'fontw','bold','fontsize',16,'xlim',[1,ntime/srate],'linew',2)
subplot(526)
[slope_highFilt, intercept] = logfit(kpk_ec_tot_highFilt{chan2plot}(:,1),kpk_ec_tot_highFilt{chan2plot}(:,2),'loglog','linewidth',3,'skipBegin',idx1,'skipEnd',length(kpk_ec_tot_highFilt{chan2plot}(:,1))-idx2)
title('ec\_high-pass filtered')
box on;grid on;
set(gca, 'fontw','bold','fontsize',16,'xlim',[1e0 1e3],'linew',2)
% text2print = ['PS = ',num2str(roundn(-slope_highFilt,-3))];
% text(180,0.01,text2print,'color','r', 'HorizontalAlignment','left')
xlabel('log_2(k)'); ylabel('log_2[p(k)]');

subplot(527)
k_min = 95; 
k_max = 325; 
idx1 = find(kpk_ec_tot_bandFilt{chan2plot}(:,1) == k_min);
idx2 = find(kpk_ec_tot_bandFilt{chan2plot}(:,1) == k_max);
ntime = size(ec_bandFilt,1);
xlab  = linspace(1,ntime/srate,ntime);
plot(xlab, ec_bandFilt(:,chan2plot,1),'linew',1);
title('ec\_band-pass filtered')
xlabel('time (s)')
ylab = ylabel('Amplitue (\mu M)');
set(ylab,'Interpreter','tex');
box on;grid on;
set(gca, 'fontw','bold','fontsize',16,'xlim',[1,ntime/srate],'linew',2)
subplot(528)
[slope_bandFilt, intercept] = logfit(kpk_ec_tot_bandFilt{chan2plot}(:,1),kpk_ec_tot_bandFilt{chan2plot}(:,2),'loglog','linewidth',3,'skipBegin',idx1,'skipEnd',length(kpk_ec_tot_bandFilt{chan2plot}(:,1))-idx2)
title('ec\_band-pass filted')
box on;grid on;
set(gca, 'fontw','bold','fontsize',16,'xlim',[1e0 1e3],'linew',2)
% text2print = ['PS = ',num2str(roundn(-slope_bandFilt, -3))];
% text(180,0.01,text2print,'color','r', 'HorizontalAlignment','left')
xlabel('log_2(k)'); ylabel('log_2[p(k)]');

chan22plot = 1;
subplot(529)
k_min = 10; 
k_max = 50;
tmp = abs(kpk_ec_tot_shuffle{chan22plot}(:,1) - k_max);
[~, idx] = min(tmp);
k_max = kpk_ec_tot_shuffle{chan22plot}(idx,1);
idx1 = find(kpk_ec_tot_shuffle{chan22plot}(:,1) == k_min);
idx2 = find(kpk_ec_tot_shuffle{chan22plot}(:,1) == k_max);
ntime = size(ec_shuffle,1);
xlab  = linspace(1,ntime/srate,ntime);
plot(xlab, ec_shuffle(:,chan22plot,1),'linew',1);
title('shuffled')
xlabel('time (s)')
ylab = ylabel('Amplitue (\mu M)');
set(ylab,'Interpreter','tex');
box on;grid on;
set(gca, 'fontw','bold','fontsize',16,'xlim',[1,ntime/srate],'linew',2)
subplot(5,2,10)
[slope_shuffle, intercept] = logfit(kpk_ec_tot_shuffle{chan22plot}(:,1),kpk_ec_tot_shuffle{chan22plot}(:,2),'loglog','linewidth',3,'skipBegin',idx1,'skipEnd',length(kpk_ec_tot_shuffle{chan22plot}(:,1))-idx2)
title('shuffled')
box on;grid on;
set(gca, 'fontw','bold','fontsize',16,'xlim',[1e0 1e3],'linew',2)
text2print = ['PS = ',num2str(roundn(-slope_shuffle, -3))];
text(180,0.01,text2print,'color','r', 'HorizontalAlignment','left')
xlabel('log_2(k)'); ylabel('log_2[p(k)]');
%% re-plot above figure based on log2
chan2plot = 2;
srate = 12.5;

figure(3); clf

subplot(221)
k_min = 37; 
k_max = 115; 
idx1 = find(kpk_ec_tot_nofilt{chan2plot}(:,1) == k_min);
idx2 = find(kpk_ec_tot_nofilt{chan2plot}(:,1) == k_max);
xlab  = linspace(1,ntime/srate,ntime);
plot(xlab,ec_nofilt(:,chan2plot,1),'linew',1);
xlabel('time (s)')
ylab = ylabel('Amplitue (\mu M)');
set(ylab,'Interpreter','tex');
title('ec\_non-filter')
box on;grid off;
set(gca, 'fontw','bold','fontsize',16,'xlim',[1,ntime/srate],'linew',2)
subplot(222)
scatter(log2(kpk_ec_tot_nofilt{chani,1}(:,1)),log2((kpk_ec_tot_nofilt{chani,1}(:,2))), dotSize, 'b','+')
ylim([log2(7e-5), log2(5e-2)])
box on
set(gca, 'fontw','bold','fontsize',16,'xlim',[2, 7],'linew',2)

% set(gca, 'XTickLabel',[])                      %# suppress current x-labels
% xt = get(gca, 'XTick');
% yl = get(gca, 'YLim');
% str = cellstr( num2str(xt(:),'2^{%d}') );      %# format x-ticks as 2^{xx}
% hTxt = text(xt, yl(ones(size(xt))), str, ...   %# create text at same locations
%     'Interpreter','tex', ...                   %# specify tex interpreter
%     'VerticalAlignment','top', ...             %# v-align to be underneath
%     'HorizontalAlignment','center');           %# h-aligh to be centered
% box on;grid on;
% % set(gca, 'fontw','bold','fontsize',16,'xlim',[log2(1e0) log2(1e3)],'linew',2)


chan22plot = 1;
subplot(223)
k_min = 10; 
k_max = 50;
tmp = abs(kpk_ec_tot_shuffle{chan22plot}(:,1) - k_max);
[~, idx] = min(tmp);
k_max = kpk_ec_tot_shuffle{chan22plot}(idx,1);
idx1 = find(kpk_ec_tot_shuffle{chan22plot}(:,1) == k_min);
idx2 = find(kpk_ec_tot_shuffle{chan22plot}(:,1) == k_max);
ntime = size(ec_shuffle,1);
xlab  = linspace(1,ntime/srate,ntime);
plot(xlab, ec_shuffle(:,chan22plot,1),'linew',1);
title('shuffled')
xlabel('time (s)')
ylab = ylabel('Amplitue (\mu M)');
set(ylab,'Interpreter','tex');
box on;grid off;
set(gca, 'fontw','bold','fontsize',16,'xlim',[1,ntime/srate],'linew',2)
subplot(224)
scatter(log2(kpk_ec_tot_shuffle{chani,1}(:,1)),log2((kpk_ec_tot_shuffle{chani,1}(:,2))), dotSize, 'b','+')
ylim([log2(7e-5), log2(5e-2)])
box on
set(gca, 'fontw','bold','fontsize',16,'xlim',[2, 7],'linew',2)

% set(gca, 'XTickLabel',[])                      %# suppress current x-labels
% xt = get(gca, 'XTick');
% yl = get(gca, 'YLim');
% str = cellstr( num2str(xt(:),'2^{%d}') );      %# format x-ticks as 2^{xx}
% hTxt = text(xt, yl(ones(size(xt))), str, ...   %# create text at same locations
%     'Interpreter','tex', ...                   %# specify tex interpreter
%     'VerticalAlignment','top', ...             %# v-align to be underneath
%     'HorizontalAlignment','center');           %# h-aligh to be centered
% box on;grid on;
% set(gca, 'fontw','bold','fontsize',16,'xlim',[log2(1e0) log2(1e3)],'linew',2)
%%
figure(3);clf
subplot(121)
plot(ec_shuffle(:,chan2plot,1),'linew',1);
title('ec\_shuffle')
subplot(1,2,2)
% [slope_shuffle, intercept] = logfit(kpk_ec_tot_shuffle{chan2plot}(:,1),kpk_ec_tot_shuffle{chan2plot}(:,2),'loglog','skip',idx)
title('ec\_shuffle')

%%%%%%%%%
%% My own way to plot
% chan2plot = 2;
% 
% figure(2); clf
% subplot(421)
% plot(ec_nofilt(:,chan2plot,1),'linew',1);
% title('ec\_non-filter')
% subplot(422)
% loglog(kpk_ec_tot_nofilt{chan2plot}(:,1),kpk_ec_tot_nofilt{chan2plot}(:,2),'b+')
% set(gca,'xlim',[0 1e3])
% title('ec\_non-filter')
% % power of scale-freeness
% k_min = 37; 
% k_max = 115; 
% idx1 = find(kpk_ec_tot_nofilt{chan2plot}(:,1) == k_min);
% idx2 = find(kpk_ec_tot_nofilt{chan2plot}(:,1) == k_max);
% tmp = polyfit(log(kpk_ec_tot_nofilt{chan2plot}(idx1:idx2,1)), log(kpk_ec_tot_nofilt{chan2plot}(idx1:idx2,2)), 1);
% power_sf = -tmp(1);
% % text2print = ['PS = ',num2str(power_sf)];
% % text(50,0.1,text2print,'color','r')
% 
% subplot(423)
% plot(ec_lowFilt(:,chan2plot,1),'linew',1);
% title('ec\_low-pass filtered')
% subplot(424)
% loglog(kpk_ec_tot_lowFilt{chan2plot}(:,1),kpk_ec_tot_lowFilt{chan2plot}(:,2),'b+')
% set(gca,'xlim',[0 1e3])
% title('ec\_low-pass filtered')
% % power of scale-freeness
% k_min = 130; 
% k_max = 350; 
% idx1 = find(kpk_ec_tot_lowFilt{chan2plot}(:,1) == k_min);
% idx2 = find(kpk_ec_tot_lowFilt{chan2plot}(:,1) == k_max);
% tmp = polyfit(log(kpk_ec_tot_lowFilt{chan2plot}(idx1:idx2,1)), log(kpk_ec_tot_lowFilt{chan2plot}(idx1:idx2,2)), 1);
% power_sf = -tmp(1);
% % text2print = ['PS = ',num2str(power_sf)];
% % text(50,0.1,text2print,'color','r')
% 
% subplot(425)
% plot(ec_highFilt(:,chan2plot,1),'linew',1);
% title('ec\_high-pass filtered')
% subplot(426)
% loglog(kpk_ec_tot_highFilt{chan2plot}(:,1),kpk_ec_tot_highFilt{chan2plot}(:,2),'b+')
% set(gca,'xlim',[0 1e3])
% title('ec\_high-pass filtered')
% % power of scale-freeness
% k_min = 37; 
% k_max = 130; 
% idx1 = find(kpk_ec_tot_highFilt{chan2plot}(:,1)==k_min);
% idx2 = find(kpk_ec_tot_highFilt{chan2plot}(:,1) == k_max);
% tmp = polyfit(log(kpk_ec_tot_highFilt{chan2plot}(idx1:idx2,1)), log(kpk_ec_tot_highFilt{chan2plot}(idx1:idx2,2)), 1);
% power_sf = -tmp(1);
% % text2print = ['PS = ',num2str(power_sf)];
% % text(50,0.1,text2print,'color','r')
% 
% subplot(427)
% plot(ec_bandFilt(:,chan2plot,1),'linew',1);
% title('ec\_band-pass filtered')
% subplot(428)
% loglog(kpk_ec_tot_bandFilt{chan2plot}(:,1),kpk_ec_tot_bandFilt{chan2plot}(:,2),'b+')
% set(gca,'xlim',[0 1e3])
% title('ec\_band-pass filted')
% % power of scale-freeness
% k_min = 95; 
% k_max = 325; 
% idx1 = find(kpk_ec_tot_bandFilt{chan2plot}(:,1) == k_min);
% idx2 = find(kpk_ec_tot_bandFilt{chan2plot}(:,1) == k_max);
% tmp = polyfit(log(kpk_ec_tot_bandFilt{chan2plot}(idx1:idx2,1)), log(kpk_ec_tot_bandFilt{chan2plot}(idx1:idx2,2)), 1);
% power_sf = -tmp(1);
% text2print = ['PS = ',num2str(power_sf)];
% text(50,0.1,text2print,'color','r')

% subplot(529)
% plot(ec_shuffle(:,chan2plot,1),'linew',1);
% title('ec\_shuffle')
% subplot(5,2,10)
% loglog(kpk_ec_tot_shuffle{chan2plot}(:,1),kpk_ec_tot_shuffle{chan2plot}(:,2),'b+')
% set(gca,'xlim',[0 1e3])
% title('ec\_shuffle')
% % power of scale-freeness
% k_min = 10; 
% k_max = 50; 
% idx1 = find(kpk_ec_tot_shuffle{chan2plot}(:,1) == k_min);
% idx2 = find(kpk_ec_tot_shuffle{chan2plot}(:,1) == k_max);
% tmp = polyfit(log(kpk_ec_tot_shuffle{chan2plot}(idx1:idx2,1)), log(kpk_ec_tot_shuffle{chan2plot}(idx1:idx2,2)), 1);
% power_sf = -tmp(1);
% text2print = ['PS = ',num2str(power_sf)];
% text(50,0.1,text2print,'color','r')

%%%%%%%%%%%%%%%%%%%
% visually compare non-filtered and high-pass filtered:
figure(3);clf
title('compare non-filtered and high-pass filtered:')
loglog(kpk_ec_tot_nofilt{chan2plot}(:,1),kpk_ec_tot_nofilt{chan2plot}(:,2),'b+');
hold on
loglog(kpk_ec_tot_highFilt{chan2plot}(:,1),kpk_ec_tot_highFilt{chan2plot}(:,2),'r+')
set(gca,'xlim',[0 1e3])
legend('original (non-filtered)','high-pass filted (at 0.01 Hz)')

% visually compare non-filtered and high-pass filtered:
figure(3);clf
title('compare non-filtered and low-pass filtered:')
loglog(kpk_ec_tot_nofilt{chan2plot}(:,1),kpk_ec_tot_nofilt{chan2plot}(:,2),'b+');
hold on
loglog(kpk_ec_tot_lowFilt{chan2plot}(:,1),kpk_ec_tot_lowFilt{chan2plot}(:,2),'r+')
set(gca,'xlim',[0 1e3])
legend('original (non-filtered)','low-pass filted (at 0.2 Hz)')

%% visually compare non-filtered and random shuffled:
% open parallel
% repeat the simulation (shuffling time series) for 100 times
% initial
num_rep = 100; % number of simulation
power_sf_shuffle_rep = zeros(num_rep,1);

parfor irep = 1: num_rep
    n = size(ec_nofilt,1); % length of time-series
    ec_shuffle_tmp = zeros(n,3);
    for blocki = 1:3
        ec_shuffle_tmp(:,blocki) = shuffle(ec_nofilt(:,chan2plot,blocki),1); % time X block
    end
    % convert time-seris to VG
    VG_ec_shuffle_tmp  = zeros(n,n,3); % time X time X block .... (only one channel!)
    for blocki = 1:3
        VG_ec_shuffle_tmp(:,:,blocki)  = lz_VG_build(ec_shuffle_tmp(:,blocki));
    end

    n = size(VG_ec_shuffle_tmp,1); % # of nodes
    % compute degree for each node
    % initiate degree
    degree_ec_shuffle_tmp  = zeros(n,3);
    for blocki =1: 3 % loop for block
        degree_ec_shuffle_tmp(:,blocki)  = squeeze(sum(VG_ec_shuffle_tmp(:,:,blocki)));
    end
    % compute p(k) 
    % initiate k and p(k) matrix for each condition
    kpk_ec_shuffle_tmp  = cell(1,3); 
    for blocki = 1: 3
        seq_bins = unique(degree_ec_shuffle_tmp(:,blocki));
        histcount = histc(degree_ec_shuffle_tmp(:,blocki), seq_bins);
        kpk_ec_shuffle_tmp(blocki) = {[seq_bins, histcount/n]};
    end
    % combine degree values across blocks
    % the following is the code developed based on http://www.mathworks.com/matlabcentral/answers/125787-merging-two-matrices-by-first-column-values
    % initiate k vs. p(k) cell array
    kpk_ec_tot_shuffle_tmp  = cell(1,1);
    q = unique([kpk_ec_shuffle_tmp{1}(:,1);kpk_ec_shuffle_tmp{2}(:,1);kpk_ec_shuffle_tmp{3}(:,1)]); % unique sorted catenate.
    Ind_blk1 = ismember(q, kpk_ec_shuffle_tmp{1}(:,1)); % All dataes in block 1
    Ind_blk2 = ismember(q, kpk_ec_shuffle_tmp{2}(:,1)); % All dataes in block 2
    Ind_blk3 = ismember(q, kpk_ec_shuffle_tmp{3}(:,1)); % All dataes in block 3
    x = zeros(length(q),4);
    x(:,1) = q; % value of k
    x(Ind_blk1,2) = kpk_ec_shuffle_tmp{1}(:,2); % p(k) in block 1
    x(not(Ind_blk1),2) = 0;
    x(Ind_blk2,3) = kpk_ec_shuffle_tmp{2}(:,2); % p(k) in block 2
    x(not(Ind_blk2),3) = 0;
    x(Ind_blk3,4) = kpk_ec_shuffle_tmp{3}(:,2); % p(k) in block 3
    x(not(Ind_blk3),4) = 0;
    kpk_ec_tot_shuffle_tmp{1} = [x(:,1),x(:,2)+x(:,3)+x(:,4)];
    kpk_ec_tot_shuffle_tmp{1}(:,2) = kpk_ec_tot_shuffle_tmp{1}(:,2)/3;
    
    % power of scale-freeness
    k_min = 10;
    k_max = 50;
    idx1 = find(kpk_ec_tot_shuffle_tmp{1}(:,1) == k_min);
    idx2 = find(kpk_ec_tot_shuffle_tmp{1}(:,1) == k_max);
%     tmp = polyfit(log(kpk_ec_tot_shuffle_tmp{1}(idx1:idx2,1)), log(kpk_ec_tot_shuffle_tmp{1}(idx1:idx2,2)), 1);
    tmp = polyfit(log(kpk_ec_tot_shuffle_tmp{1}(idx1:end,1)), log(kpk_ec_tot_shuffle_tmp{1}(idx1:end,2)), 1);
    
    % plot
%     figure; loglog((kpk_ec_tot_shuffle_tmp{1}(:,1)), log(kpk_ec_tot_shuffle_tmp{1}(:,2)))
    power_sf_shuffle_rep(irep) = -tmp(1);
    
end
clear blocki
power_sf_shuffle_mn  = mean(power_sf_shuffle_rep)
power_sf_shuffle_std = std(power_sf_shuffle_rep)

figure(4); boxplot(power_sf_shuffle_rep);
figure(5); clf; histfit(power_sf_shuffle_rep); 
set(gca,'fontsize',52,'fontweight','bold','linew',2)
xlabel('Power of Scale-freeness');
ylabel('Frequency');
box on,grid on;
xlim([3 4.3]); 
vline(3.097,'g');
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot4paper
subject = 'l9e11'
%%%% 1. no filtered:
eo_nofilt = load('C:\Users\Li_Lab537\Documents\Projects\GoNoGo_NIRx_EEG\l9e11\NIRS\2015-06-26_003_EO\exportData\wideBand\NIRS-2015-06-26_003_oxyhb_T1to7503_C1to14.txt');
tmp = eo_nofilt; 
clear eo_nofilt;
eo_nofilt.blk1 = tmp(564 : 564+2124,:);
eo_nofilt.blk2 = tmp(564+2124 : 564+2124+2124,:);
eo_nofilt.blk3 = tmp(564+2124+2124 : 564+2124+2124+2124,:);
eo_nofilt = cat(3,eo_nofilt.blk1,eo_nofilt.blk2,eo_nofilt.blk3);
clear tmp;
figure(1);clf
subplot(211)
plot(eo_nofilt(:,chan2plot,blk2plot))
% bandpass filter
srate = 12.5;
order = 4;
lobound = 0.01;
hibound = 0.2;
type = 'bandpass';
for blocki = 1:3
    eo_bandFilt(:,:,blocki)  = lz_butterworth(eo_nofilt(:,:,blocki),srate,order,lobound,hibound,type);
end
figure(1);clf
subplot(211)
ntime = size(eo_bandFilt,1);
xlab  = linspace(1,ntime/srate,ntime);
plot(xlab,eo_bandFilt(:,chan2plot,blk2plot),'linew',3)
xlabel('time (s)')
ylab = ylabel('Amplitue ($\mu$ M)');
% set(ylab,'Interpreter','latex')
set(gca,'xlim',[1,ntime/srate],'linew',2, 'fontsize',22,'fontweight','bold')
subplot(212)
autocorr(eo_bandFilt(:,chan2plot,blk2plot),40)
set(gca,'linew',2, 'fontsize',22,'fontweight','bold')
clear srate order lobound hibound type blocki eo ec rt irt mix;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% replot fig. 4 for journal paper
%% with base of 2
chan2plot = 2;
srate = 12.5;

figure(12); clf

subplot(221)
k_min = 37; 
k_max = 115; 
idx1 = find(kpk_ec_tot_nofilt{chan2plot}(:,1) == k_min);
idx2 = find(kpk_ec_tot_nofilt{chan2plot}(:,1) == k_max);
xlab  = linspace(1,ntime/srate,ntime);
plot(xlab,ec_nofilt(:,chan2plot,1),'linew',1);
xlabel('time (s)')
ylab = ylabel('Amplitue (\mu M)');
set(ylab,'Interpreter','tex');
title('ec\_non-filter')
box on;grid on;
set(gca, 'fontw','bold','fontsize',16,'xlim',[1,ntime/srate],'linew',2)
subplot(222)

[slope_nofilt, intercept] = logfit_log2(kpk_ec_tot_nofilt{chan2plot}(:,1),kpk_ec_tot_nofilt{chan2plot}(:,2),'loglog','linewidth',3,'skipBegin',idx1,'skipEnd',length(kpk_ec_tot_nofilt{chan2plot}(:,1))-idx2)
title('ec\_non-filter')
box on;grid off;
set(gca, 'fontw','bold','fontsize',16,'xlim',[1e0 1e3],'linew',2)
% text2print = ['PS = ',num2str(roundn(-slope_nofilt,-3))];
% text(180,0.01,text2print,'color','r', 'HorizontalAlignment','left')
xlabel('log_2(k)'); ylabel('log_2[p(k)]');


chan22plot = 2;
subplot(223)
k_min = 10; 
k_max = 50;
tmp = abs(kpk_ec_tot_shuffle{chan22plot}(:,1) - k_max);
[~, idx] = min(tmp);
k_max = kpk_ec_tot_shuffle{chan22plot}(idx,1);
idx1 = find(kpk_ec_tot_shuffle{chan22plot}(:,1) == k_min);
idx2 = find(kpk_ec_tot_shuffle{chan22plot}(:,1) == k_max);
ntime = size(ec_shuffle,1);
xlab  = linspace(1,ntime/srate,ntime);
plot(xlab, ec_shuffle(:,chan22plot,1),'linew',1);
title('shuffled')
xlabel('time (s)')
ylab = ylabel('Amplitue (\mu M)');
set(ylab,'Interpreter','tex');
box on;grid on;
set(gca, 'fontw','bold','fontsize',16,'xlim',[1,ntime/srate],'linew',2)
subplot(2,2,4)
[slope_shuffle, intercept] = logfit_log2(kpk_ec_tot_shuffle{chan22plot}(:,1),kpk_ec_tot_shuffle{chan22plot}(:,2),'loglog','linewidth',3,'skipBegin',idx1,'skipEnd',length(kpk_ec_tot_shuffle{chan22plot}(:,1))-idx2)
title('shuffled')
box on;grid off;
set(gca, 'fontw','bold','fontsize',16,'xlim',[1e0 1e3],'linew',2)
% text2print = ['PS = ',num2str(roundn(-slope_shuffle, -3))];
% text(180,0.01,text2print,'color','r', 'HorizontalAlignment','left')
xlabel('log_2(k)'); ylabel('log_2[p(k)]');
% 
%     xlim([log2(100), log2(500)])
%     ylim([log2(7e-5), log2(2e-2)])    %     legend('ec','eo','rt','mix')
    set(gca, 'XTickLabel',[])                      %# suppress current x-labels
    xt = get(gca, 'XTick');
    xt = log2(xt)
    yl = get(gca, 'YLim');
    str = cellstr( num2str(xt(:),'2^{%d}') );      %# format x-ticks as 2^{xx}
    hTxt = text(xt, yl(ones(size(xt))), str, ...   %# create text at same locations
        'Interpreter','tex', ...                   %# specify tex interpreter
        'VerticalAlignment','top', ...             %# v-align to be underneath
        'HorizontalAlignment','center');           %# h-aligh to be centered
    
%     set(gca, 'YTickLabel',[])                      %# suppress current x-labels
    yt = get(gca, 'YTick');
    xl = get(gca, 'XLim');
    str = cellstr( num2str(yt(:),'2^{%d}') );      %# format x-ticks as 2^{xx}
    hTxt = text(yt, xl(ones(size(yt))), str, ...   %# create text at same locations
        'Interpreter','tex', ...                   %# specify tex interpreter
        'VerticalAlignment','middle', ...             %# v-align to be underneath
        'HorizontalAlignment','left');           %# h-aligh to be centered

