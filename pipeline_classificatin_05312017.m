%%%% 05/24/2017
%%%% performing classification related study using features extracted
%%%% through script "pipeline_feature_05222017"

%% load power of scale-freeness
%%
if ismac
    load('/Users/lizhu/Dropbox/projects/VG/VG_OneDrive/results_journal/power_sf_mat.mat')
elseif ~ispc
    load('/home/lz206/Dropbox/projects/VG/VG_OneDrive/results_journal/power_sf_mat.mat')
end
%% normalize within each channel and each subject
%% 
for iSu = 1:9
    power_sf_norm(:,:,iSu) = normr(power_sf_mat(:,:,iSu));
end
%% check point: plot histogram for each subject, based on the normalized data
%% average power of scale-freeness and bar plot across subjects
power_sf_mnCrossSub = mean(power_sf_norm,3);
power_sf_stdCrossSub = std(power_sf_norm,0,3);
figure(15);clf
for chani = 1: 14
    subplot(2,7,chani)
    barwitherr(power_sf_stdCrossSub(chani,:,:), power_sf_mnCrossSub(chani,:,:))
    set(gca,'XTickLabel',{'EC','EO','RT','MIX'})
    title(['Channel ', num2str(chani)])
    xlabel('Condition')
    ylabel('r')
    ylim([0,1])
    box on; 
    grid on;
    set(gca, 'fontsize',12, 'linew',2)
end
%% pool resting-state and task-based data for each subject
psf_rs = power_sf_norm(:,1:2,:);
psf_rs = reshape(psf_rs, 14, [], 1);
psf_tb = power_sf_norm(:,3:4,:);
psf_tb = reshape(psf_tb, 14, [], 1);
%% t-test on PS for across subjects: resting-state vs. task-based
%%
for iCh = 1:14
    t_Ch(iCh) = lz_ttest2(psf_rs(iCh,:), psf_tb(iCh,:));
end
[~, Ch_order] = sort(t_Ch, 'descend');
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% organize figure for paper: kNN
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sampData = cat( 1, psf_rs', psf_tb' );
sampLabl = [cellstr(repmat('Resting-State', 18, 1)); cellstr(repmat('Task-Based', 18, 1))];
%%%% initiate
ac_kNN = nan(6,14);
se_kNN = nan(6,14);
sp_kNN = nan(6,14);
for iNumNeighbor = 1: 6
    for num_Ch_for_feature = 1: 14
        ind_Ch_for_feature = Ch_order(1:num_Ch_for_feature);
        [ac_kNN(iNumNeighbor,num_Ch_for_feature), se_kNN(iNumNeighbor,num_Ch_for_feature),...
            sp_kNN(iNumNeighbor,num_Ch_for_feature)] = ...
            lz_knn_leave_two_out_PSF( sampData(:,ind_Ch_for_feature), sampLabl, iNumNeighbor );
    end
end
%%%% results
figure(11); clf;
subplot(131)
plot(ac_kNN', 'linew', 2)
set(gca, 'linew', 3, 'fontweight', 'bold', 'fontsize', 16)
xlabel('#Channel for Classification')
ylabel('Accuracy (%)')
xlim([1 14])
ylim([0 1])
legend('k = 1', 'k = 2', 'k = 3', 'k = 4', 'k = 5', 'k = 6', 'location', 'southeast')
title('kNN Accuracy')

subplot(132)
plot(se_kNN', 'linew', 2)
set(gca, 'linew', 3, 'fontweight', 'bold', 'fontsize', 16)
xlabel('#Channel for Classification')
ylabel('Sensitivity (%)')
xlim([1 14])
ylim([0 1])
legend('k = 1', 'k = 2', 'k = 3', 'k = 4', 'k = 5', 'k = 6', 'location', 'southeast')
title('kNN Sensitivity')

subplot(133)
plot(sp_kNN', 'linew', 2)
set(gca, 'linew', 3, 'fontweight', 'bold', 'fontsize', 16)
xlabel('#Channel for Classification')
ylabel('Specificity (%)')
xlim([1 14])
ylim([0 1])
legend('k = 1', 'k = 2', 'k = 3', 'k = 4', 'k = 5', 'k = 6', 'location', 'southeast')
title('kNN Specificity')
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% organize figure for paper: SVM
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
boxconstraint= 30;

tic
% kernel function
kernelArray = {'linear', 'rbf', 'quadratic', 'polynomial', 'polynomial-4'};
%%%% initiate
ac_SVM = nan(5,14);
se_SVM = nan(5,14);
sp_SVM = nan(5,14);

for iKernel = 1: 5
    
    kernel = kernelArray{iKernel};
    for num_Ch_for_feature = 1: 14
        ind_Ch_for_feature = Ch_order(1:num_Ch_for_feature);
        [ac_SVM(iKernel,num_Ch_for_feature), se_SVM(iKernel,num_Ch_for_feature), ...
            sp_SVM(iKernel,num_Ch_for_feature)] = ...
            lz_svm_leave_two_out_PSF_2( sampData(:,ind_Ch_for_feature), ...
            sampLabl, kernel, boxconstraint );
    end
end
toc

%%%% results
%%%% results
figure(12); clf;
% subplot(131)
plot(ac_SVM', 'linew', 2)
set(gca, 'linew', 3, 'fontweight', 'bold', 'fontsize', 16)
xlabel('#Channel for Classification')
ylabel('Accuracy (%)')
xlim([1 14])
ylim([0 1])
legend('SVM, linear kernel', 'SVM, RBF kernel', 'SVM, quadratic kernel', ...
    'SVM, 3^{rd} order polynomial kernel', 'SVM, 4^{th} order polynomial kernel', ...
    'location', 'southeast')
title('SVM Accuracy')


% subplot(132)
% plot(se_SVM', 'linew', 2)
% set(gca, 'linew', 3, 'fontweight', 'bold', 'fontsize', 16)
% xlabel('#Channel for Classification')
% ylabel('Sensitivity (%)')
% xlim([1 14])
% ylim([0 1])
% legend('SVM, linear kernel', 'SVM, RBF kernel', 'SVM, quadratic kernel', ...
%     'SVM, 3^{rd} order polynomial kernel', 'SVM, 4^{th} order polynomial kernel', ...
%     'location', 'southeast')
% title('SVM Sensitivity')
% 
% subplot(133)
% plot(sp_SVM', 'linew', 2)
% set(gca, 'linew', 3, 'fontweight', 'bold', 'fontsize', 16)
% xlabel('#Channel for Classification')
% ylabel('Specificity (%)')
% xlim([1 14])
% ylim([0 1])
% legend('SVM, linear kernel', 'SVM, RBF kernel', 'SVM, quadratic kernel', ...
%     'SVM, 3^{rd} order polynomial kernel', 'SVM, 4^{th} order polynomial kernel', ...
%     'location', 'southeast')
% title('SVM Specificity')
%%
% visualize SVM using low-dimensional features
kernel = 'linear';
num_Ch_for_feature = 3;
ind_Ch_for_feature = Ch_order(1:num_Ch_for_feature);

%%%%% Train the SVM Classifier
if num_Ch_for_feature == 2
    ifplot = 'true';
elseif num_Ch_for_feature == 3
    ifplot = 'False';
end

if strcmp(kernel, '@lz_svm_sigmoid')
    Md1 = svmtrain(sampData(:,ind_Ch_for_feature), sampLabl, ...
        'kernel_function', @lz_svm_sigmoid, 'showplot', ifplot, 'boxconstraint', boxconstraint);
else
    Md1 = svmtrain(sampData(:,ind_Ch_for_feature), sampLabl, ...
        'kernel_function', kernel, 'rbf_sigma', .85, 'showplot', ifplot, 'boxconstraint', boxconstraint);
end

if num_Ch_for_feature == 3
    figure(13); clf;
    lz_svm_3d_matlab_vis(Md1, sampData(:,ind_Ch_for_feature), sampLabl);
end

legend('location', 'northwest')
xlabel('Channel 2'); ylabel('Channel 13'); zlabel('Channel 4');
set(gca, 'linew', 1.5, 'fontweight', 'bold', 'fontsize', 16)
grid off
%%
%% 06/04/2017 
%% use LDA
tic
% kernel function
DiscrimType = {'linear', 'quadratic'};
%%%% initiate
ac_lda = nan(2,14);
se_lda = nan(2,14);
sp_lda = nan(2,14);

for iType = 1: 2
    
    thisType = DiscrimType{iType};
    for num_Ch_for_feature = 1: 14
        ind_Ch_for_feature = Ch_order(1:num_Ch_for_feature);
        [ac_lda(iType,num_Ch_for_feature), se_lda(iType,num_Ch_for_feature), ...
            sp_lda(iType,num_Ch_for_feature)] = ...
            lz_lda_leave_two_out_PSF( sampData(:,ind_Ch_for_feature), ...
            sampLabl, thisType );
    end
end
toc

%%%% results
figure(14); clf;
% subplot(131)
plot(ac_lda', 'linew', 2)
set(gca, 'linew', 3, 'fontweight', 'bold', 'fontsize', 16)
xlabel('#Channel for Classification')
ylabel('Accuracy (%)')
xlim([1 14])
ylim([0 1])
legend('linear', 'quadratic', 'location', 'southwest')
title('LDA Accuracy')

subplot(132)
plot(se_lda', 'linew', 2)
set(gca, 'linew', 3, 'fontweight', 'bold', 'fontsize', 16)
xlabel('#Channel for Classification')
ylabel('Sensitivity (%)')
xlim([1 14])
ylim([0 1])
legend('linear', 'quadratic', 'location', 'southwest')
title('LDA Sensitivity')

subplot(133)
plot(sp_lda', 'linew', 2)
set(gca, 'linew', 3, 'fontweight', 'bold', 'fontsize', 16)
xlabel('#Channel for Classification')
ylabel('Specificity (%)')
xlim([1 14])
ylim([0 1])
legend('linear', 'quadratic', 'location', 'southwest')
title('LDA Specificity')


















