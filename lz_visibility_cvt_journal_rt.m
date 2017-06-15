%%%% 05/24/2016 
%%%% compute visibility graph for individual subjects for MIX condition.
%%%% filter 0.01~0.2 Hz. no downsampling


function [VG_rt] = lz_visibility_cvt_journal_rt(subject)

%%%% Input:
%%%% subject: 'm1n03, i6a08,i3i03,n8n10,a2u04,i3h06,i4n06,u8n09,l9e11'.
%%%% filtOrNot: 2: filter and no downsampling, 1: filter and downsampling, 0: no filter and no downsampling.

%%%% Output:
%%%% VG_XX: visibility graphs. Time X Time X Channel X Block

%=====================================================================
% prepare for the data
%%%% load stim and data
srate = 12.5;
nbchan = 14;

if strcmp(subject, 'l9e11')
    if ~ispc
        rtirt_stim        = load('/home/lz206/Dropbox/projects/GNG_LI_LAB537/l9e11/NIRS/2015-06-26_005_RTIRT/NIRS-2015-06-26_005.evt');
       rt                = load('/home/lz206/Dropbox/projects/GNG_LI_LAB537/l9e11/NIRS/2015-06-26_005_RTIRT/exportData/wideBandRT/NIRS-2015-06-26_005_oxyhb_T1to16055_C1to14.txt');
    else
        rtirt_stim        = load('C:\Users\Li_Lab537\Documents\Projects\GoNoGo_NIRx_EEG\l9e11\NIRS\2015-06-26_005_RTIRT\NIRS-2015-06-26_005.evt');
        rt = load('C:\Users\Li_Lab537\Documents\Projects\GoNoGo_NIRx_EEG\l9e11\NIRS\2015-06-26_005_RTIRT\exportData\wideBandRT\NIRS-2015-06-26_005_oxyhb_T1to16055_C1to14.txt');
    end
end

if strcmp(subject, 'i4n06')
    if ~ispc
        rtirt_stim        = load('/home/lz206/Dropbox/projects/GNG_LI_LAB537/i4n06/NIRS/2015-06-30_005_RTIRT/NIRS-2015-06-30_005.evt');
        rt = load('/home/lz206/Dropbox/projects/GNG_LI_LAB537/i4n06/NIRS/2015-06-30_005_RTIRT/exportData/wideBandRT/NIRS-2015-06-30_005_oxyhb_T1to12675_C1to14.txt');
    else
        rtirt_stim        = load('C:\Users\Li_Lab537\Documents\Projects\GoNoGo_NIRx_EEG\i4n06\NIRS\2015-06-30_005_RTIRT\NIRS-2015-06-30_005.evt');
        rt = load('C:\Users\Li_Lab537\Documents\Projects\GoNoGo_NIRx_EEG\i4n06\NIRS\2015-06-30_005_RTIRT\exportData\wideBandRT\NIRS-2015-06-30_005_oxyhb_T1to12675_C1to14.txt');
    end
end

if strcmp(subject, 'u8n09')
    if ~ispc
        rtirt_stim        = load('/home/lz206/Dropbox/projects/GNG_LI_LAB537/u8n09/NIRS/2015-06-25_003_RTIRT/NIRS-2015-06-25_003.evt');
        rt = load('/home/lz206/Dropbox/projects/GNG_LI_LAB537/u8n09/NIRS/2015-06-25_003_RTIRT/exportData/wideBandRT/NIRS-2015-06-25_003_oxyhb_T1to12485_C1to14.txt');
    else
        rtirt_stim        = load('C:\Users\Li_Lab537\Documents\Projects\GoNoGo_NIRx_EEG\u8n09\NIRS\2015-06-25_003_RTIRT\NIRS-2015-06-25_003.evt');
        rt = load('C:\Users\Li_Lab537\Documents\Projects\GoNoGo_NIRx_EEG\u8n09\NIRS\2015-06-25_003_RTIRT\exportData\wideBandRT\NIRS-2015-06-25_003_oxyhb_T1to12485_C1to14.txt');
    end
end

if strcmp(subject, 'i3h06')
    if ~ispc
        rtirt_stim        = load('/home/lz206/Dropbox/projects/GNG_LI_LAB537/i3h06/NIRs/2015-07-07_003_RTIRT/NIRS-2015-07-07_003.evt');
        rt = load('/home/lz206/Dropbox/projects/GNG_LI_LAB537/i3h06/NIRs/2015-07-07_003_RTIRT/exportData/wideBandRT/NIRS-2015-07-07_003_oxyhb_T1to13240_C1to14.txt');
    else
        rtirt_stim        = load('C:\Users\Li_Lab537\Documents\Projects\GoNoGo_NIRx_EEG\i3h06\NIRs\2015-07-07_003_RTIRT\NIRS-2015-07-07_003.evt');
        rt = load('C:\Users\Li_Lab537\Documents\Projects\GoNoGo_NIRx_EEG\i3h06\NIRs\2015-07-07_003_RTIRT\exportData\wideBandRT\NIRS-2015-07-07_003_oxyhb_T1to13240_C1to14.txt');
    end
end

if strcmp(subject, 'a2u04')
    if ~ispc
        rtirt_stim        = load('/home/lz206/Dropbox/projects/GNG_LI_LAB537/a2u04/NIRs/2015-07-08_005_RTIRT/NIRS-2015-07-08_005.evt');
        rt = load('/home/lz206/Dropbox/projects/GNG_LI_LAB537/a2u04/NIRs/2015-07-08_005_RTIRT/exportData/wideBandRT/NIRS-2015-07-08_005_oxyhb_T1to13030_C1to14.txt');
    else
        rtirt_stim        = load('C:\Users\Li_Lab537\Documents\Projects\GoNoGo_NIRx_EEG\a2u04\NIRs\2015-07-08_005_RTIRT\NIRS-2015-07-08_005.evt');
        rt = load('C:\Users\Li_Lab537\Documents\Projects\GoNoGo_NIRx_EEG\a2u04\NIRs\2015-07-08_005_RTIRT\exportData\wideBandRT\NIRS-2015-07-08_005_oxyhb_T1to13030_C1to14.txt');
    end
end

if strcmp(subject, 'n8n10')
    if ~ispc
        rtirt_stim        = load('/home/lz206/Dropbox/projects/GNG_LI_LAB537/n8n10/NIRS/2015-07-10_004_RTIRT/NIRS-2015-07-10_004.evt');
        rt = load('/home/lz206/Dropbox/projects/GNG_LI_LAB537/n8n10/NIRS/2015-07-10_004_RTIRT/exportData/wideBandRT/NIRS-2015-07-10_004_oxyhb_T1to12810_C1to14.txt');
    else
        rtirt_stim        = load('C:\Users\Li_Lab537\Documents\Projects\GoNoGo_NIRx_EEG\n8n10\NIRS\2015-07-10_004_RTIRT\NIRS-2015-07-10_004.evt');
        rt = load('C:\Users\Li_Lab537\Documents\Projects\GoNoGo_NIRx_EEG\n8n10\NIRS\2015-07-10_004_RTIRT\exportData\wideBandRT\NIRS-2015-07-10_004_oxyhb_T1to12810_C1to14.txt');
    end
end

if strcmp(subject, 'i3i03')
    if ~ispc
        rtirt_stim        = load('/home/lz206/Dropbox/projects/GNG_LI_LAB537/i3i03/NIRS/2015-07-10_009_RTIRT/NIRS-2015-07-10_009.evt');
        rt = load('/home/lz206/Dropbox/projects/GNG_LI_LAB537/i3i03/NIRS/2015-07-10_009_RTIRT/exportData/wideBandRT/NIRS-2015-07-10_009_oxyhb_T1to12650_C1to14.txt');
    else
        rtirt_stim        = load('C:\Users\Li_Lab537\Documents\Projects\GoNoGo_NIRx_EEG\i3i03\NIRS\2015-07-10_009_RTIRT\NIRS-2015-07-10_009.evt');
        rt = load('C:\Users\Li_Lab537\Documents\Projects\GoNoGo_NIRx_EEG\i3i03\NIRS\2015-07-10_009_RTIRT\exportData\wideBandRT\NIRS-2015-07-10_009_oxyhb_T1to12650_C1to14.txt');
    end
end

if strcmp(subject, 'i6a08')
    if ~ispc
        rtirt_stim        = load('/home/lz206/Dropbox/projects/GNG_LI_LAB537/i6a08/NIRS/2015-07-11_003_RTIRT/NIRS-2015-07-11_003.evt');
        rt = load('/home/lz206/Dropbox/projects/GNG_LI_LAB537/i6a08/NIRS/2015-07-11_003_RTIRT/exportData/wideBandRT/NIRS-2015-07-11_003_oxyhb_T1to12455_C1to14.txt');
    else
        rtirt_stim        = load('C:\Users\Li_Lab537\Documents\Projects\GoNoGo_NIRx_EEG\i6a08\NIRS\2015-07-11_003_RTIRT\NIRS-2015-07-11_003.evt');
        rt = load('C:\Users\Li_Lab537\Documents\Projects\GoNoGo_NIRx_EEG\i6a08\NIRS\2015-07-11_003_RTIRT\exportData\wideBandRT\NIRS-2015-07-11_003_oxyhb_T1to12455_C1to14.txt');
    end
end

if strcmp(subject, 'm1n03')
    if ~ispc
        rtirt_stim        = load('/home/lz206/Dropbox/projects/GNG_LI_LAB537/m1n03/NIRS/2015-07-14_003_RTIRT/NIRS-2015-07-14_003.evt');
        rt = load('/home/lz206/Dropbox/projects/GNG_LI_LAB537/m1n03/NIRS/2015-07-14_003_RTIRT/exportData/wideBandRT/NIRS-2015-07-14_003_oxyhb_T1to12320_C1to14.txt');
    else
        rtirt_stim        = load('C:\Users\Li_Lab537\Documents\Projects\GoNoGo_NIRx_EEG\i6a08\NIRS\2015-07-11_003_RTIRT\NIRS-2015-07-11_003.evt');
        rt = load('C:\Users\Li_Lab537\Documents\Projects\GoNoGo_NIRx_EEG\i6a08\NIRS\2015-07-11_003_RTIRT\exportData\wideBandRT\NIRS-2015-07-11_003_oxyhb_T1to12455_C1to14.txt');
    end
end

%%%% stim
% convert the binary .evt file to decible
stim_deci = zeros(size(rtirt_stim,1),1);

for stimi = 1:length(rtirt_stim)
    stim_deci(stimi) = rtirt_stim(stimi,2)*1 + rtirt_stim(stimi,3)*2 + rtirt_stim(stimi,4)*4 + rtirt_stim(stimi,5)*8;
end
rtirt_stim = [rtirt_stim, stim_deci];

stim_rt = rtirt_stim(find(rtirt_stim(:,10)==7) );

%%%% initial blocked data
tmp = rt;
rt = zeros(170*srate, 14, 3);

for blki = 1:3
    for chani = 1:nbchan
        rt(:,chani,blki) = tmp(stim_rt(blki)-14*srate+1:stim_rt(blki)+156*srate,chani);
        rt(:,chani,blki) = rt(:,chani,blki) - mean(rt(4*srate:14*srate,chani,blki));
    end
end

% filter each signal with [0.01-0.2] Hz if needed
srate = 12.5;
order = 4;
lobound = 0.01;
hibound = 0.2;
type = 'bandpass';
for blocki = 1:3
    rt_filt(:,:,blocki) = lz_butterworth(rt(:,:,blocki),srate,order,lobound,hibound,type);
end
clear srate order lobound hibound type blocki eo ec rt irt mix;
    
% downsampling the time-series
% rt  = downsample(rt_filt,3);
rt = rt_filt;
clear eo_filt ec_filt rt_filt irt_filt mix_filt;

% convert time-seris to VG
n = size(rt,1); % length of time-series
m = size(rt,2); % # of channel

VG_rt = zeros(n,n,m,3);

for blocki = 1:3
    VG_rt(:,:,:,blocki) = lz_VG_build_2(rt(:,:,blocki));
end
clear blocki 

% save VG
eval(['VG_',subject,'_rt = VG_rt;']);

eval(['save(''/home/lz206/Dropbox/projects/VG/VG_OneDrive/results_journal/VG_',subject,'_rt_05222017_noDownSamp.mat'',''VG_',subject,'_rt'');'])




