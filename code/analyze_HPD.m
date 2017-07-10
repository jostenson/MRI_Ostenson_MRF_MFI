%% analyze ROIs for mean and stnd dev of T1, T2

disp('Analyzing MRI system phantom results...')
tic;

% set params

nMRF = 4; % number of MRF datasets
nSl = 2; % number of slices
nROI = 14; % number of ROIs
tempMsr = 22.0; % deg C
tempStnd = 293.0 - 273.15; % deg C; ref temperature for HPD relaxation value
delTemp = tempMsr - tempStnd; % deg C

HPDlegendtxt = {'shim no MFI','shim MFI','B0 het no MFI','B0 het MFI'};

fn_1 = 'MRF_proc_HPD_top_shim.mat';
fn_2 = 'MRF_proc_HPD_mid_shim.mat';
fn_3 = 'MRF_proc_HPD_top_shim_MFI.mat';
fn_4 = 'MRF_proc_HPD_mid_shim_MFI.mat';
fn_5 = 'MRF_proc_HPD_top_B0het.mat';
fn_6 = 'MRF_proc_HPD_mid_B0het.mat';
fn_7 = 'MRF_proc_HPD_top_B0het_MFI.mat';
fn_8 = 'MRF_proc_HPD_mid_B0het_MFI.mat';

fn_rois = 'ROImaps_HPD_2017.mat';

roistats = struct();
roistats.HPDlegendtxt = HPDlegendtxt;
% HPD system phantom specs at 3 T
% T1 contrast spheres
roistats.HPD.sl1.T1.mean = [1989; 1454; 984.1; 706; 496.7; 351.5; 247.13; ...
    175.3; 125.9; 89.0; 62.7; 44.53; 30.84; 21.719];
roistats.HPD.sl1.T1.stdev = [1.0; 2.5; 0.33; 1.5; 0.41; 0.91; 0.086; 0.11; ...
    0.33; 0.17; 0.13; 0.090; 0.016; 0.0054];
roistats.HPD.sl1.NiCl2_mM = [0.299; 0.623; 1.072; 1.720; 2.617; 3.912; 5.731; ...
    8.297; 11.936; 17.070; 24.326; 34.590; 49.122; 69.680];
roistats.HPD.sl1.coeffT1_NiCl2 = 0.004; % slope; 1/(mMs)/degC

roistats.HPD.sl1.T2.mean = [1465; 1076; 717.9; 510.1; 359.6; 255.5; 180.8; ...
    127.3; 90.3; 64.3; 45.7; 31.86; 22.38; 15.83];
roistats.HPD.sl1.T2.stdev = [1.0; 1.8; 1.12; 1.36; 0.22; 0.07; 0.04; 0.14; 0.14; ...
    0.05; 0.12; 0.02; 0.02; 0.03];
roistats.HPD.sl1.coeffT2_NiCl2 = 0.0; % slope

% T2 contrast spheres
roistats.HPD.sl2.T1.mean = [2480; 2173; 1907; 1604; 1332; 1044; 801.7; 608.6; ...
    458.4; 336.5; 244.2; 176.6; 126.9; 90.9];
roistats.HPD.sl2.T1.stdev = [10.8; 14.7; 10.3; 7.2; 0.8; 3.2; 1.70; 1.03; 0.33; ...
    0.18; 0.09; 0.09; 0.03; 0.05];
roistats.HPD.sl2.MnCl2_mM = [0.013; 0.021; 0.031; 0.047; 0.069; 0.101; 0.145; ...
    0.207; 0.296; 0.421; 0.599; 0.849; 1.104; 1.704];
roistats.HPD.sl2.coeffT1_MnCl2 = [0.0029 -.28]; % quad term, linear term, 1/(mMs)/degC^2 ...

roistats.HPD.sl2.T2.mean = [581.3; 403.5; 278.1; 190.94; 133.27; 96.89; 64.07; ...
    46.42; 31.97; 22.56; 15.813; 11.237; 7.911; 5.592];
roistats.HPD.sl2.T2.stdev = [0.39; 0.55; 0.28; 0.011; 0.073; 0.049; 0.034; 0.014; ...
    0.083; 0.012; 0.0061; 0.0057; 0.0037; 0.0055];
roistats.HPD.sl2.coeffT2_MnCl2 = [0.01 -1.25]; % quad term, linear term, 1/(mMs)/degC^2 ...


%% load data

load([datadir fn_1]);
T1_d1_sl1_m = output_MRF_match.T1_map;
T2_d1_sl1_m = output_MRF_match.T2_map;
magimg_d1_sl1_m = abs(sum(input_MRF_dp.img_AR_comb,3));
load([datadir fn_2]);
T1_d1_sl2_m = output_MRF_match.T1_map;
T2_d1_sl2_m = output_MRF_match.T2_map;
magimg_d1_sl2_m = abs(sum(input_MRF_dp.img_AR_comb,3));
load([datadir fn_3]);
T1_d2_sl1_m = output_MRF_match.T1_map;
T2_d2_sl1_m = output_MRF_match.T2_map;
magimg_d2_sl1_m = abs(sum(input_MRF_dp.img_AR_comb,3));
load([datadir fn_4]);
T1_d2_sl2_m = output_MRF_match.T1_map;
T2_d2_sl2_m = output_MRF_match.T2_map;
magimg_d2_sl2_m = abs(sum(input_MRF_dp.img_AR_comb,3));
load([datadir fn_5]);
T1_d3_sl1_m = output_MRF_match.T1_map;
T2_d3_sl1_m = output_MRF_match.T2_map;
magimg_d3_sl1_m = abs(sum(input_MRF_dp.img_AR_comb,3));
load([datadir fn_6]);
T1_d3_sl2_m = output_MRF_match.T1_map;
T2_d3_sl2_m = output_MRF_match.T2_map;
magimg_d3_sl2_m = abs(sum(input_MRF_dp.img_AR_comb,3));
load([datadir fn_7]);
T1_d4_sl1_m = output_MRF_match.T1_map;
T2_d4_sl1_m = output_MRF_match.T2_map;
magimg_d4_sl1_m = abs(sum(input_MRF_dp.img_AR_comb,3));
load([datadir fn_8]);
T1_d4_sl2_m = output_MRF_match.T1_map;
T2_d4_sl2_m = output_MRF_match.T2_map;
magimg_d4_sl2_m = abs(sum(input_MRF_dp.img_AR_comb,3));



%% MRI system temp correction

% temperature correction, assumes departure around each reported relaxation
% value above matches temperature behavior in spec, which is scaled by
% concentration
r1_0_v = 1./(roistats.HPD.sl1.T1.mean/1000); % 1/s
delr1_v = delTemp * roistats.HPD.sl1.coeffT1_NiCl2 * roistats.HPD.sl1.NiCl2_mM;
r1_new_v = r1_0_v + delr1_v;
roistats.HPD.sl1.T1.meanTcor = (1./r1_new_v)*1000; % ms

r2_0_v = 1./(roistats.HPD.sl1.T2.mean/1000); % 1/s
delr2_v = delTemp * roistats.HPD.sl1.coeffT2_NiCl2 * roistats.HPD.sl1.NiCl2_mM;
r2_new_v = r2_0_v + delr2_v;
roistats.HPD.sl1.T2.meanTcor = (1./r2_new_v)*1000; % ms

r1_0_v = 1./(roistats.HPD.sl2.T1.mean/1000); % 1/s
delr1_v = (delTemp^2 * roistats.HPD.sl2.coeffT1_MnCl2(1) + ...
    delTemp * roistats.HPD.sl2.coeffT1_MnCl2(2))* roistats.HPD.sl2.MnCl2_mM;
r1_new_v = r1_0_v + delr1_v;
roistats.HPD.sl2.T1.meanTcor = (1./r1_new_v)*1000; % ms

r2_0_v = 1./(roistats.HPD.sl2.T2.mean/1000); % 1/s
delr2_v = (delTemp^2 * roistats.HPD.sl2.coeffT2_MnCl2(1) + ...
    delTemp * roistats.HPD.sl2.coeffT2_MnCl2(2))* roistats.HPD.sl2.MnCl2_mM;
r2_new_v = r2_0_v + delr2_v;
roistats.HPD.sl2.T2.meanTcor = (1./r2_new_v)*1000; % ms

%% load ROIs

load([datadir_B0 fn_rois])

%% grab ROI data and statistics

n = 1;
for ii = 1:nMRF
    for jj = 1:nSl
        for kk = 1:nROI
            
            tmp_roi = squeeze(rois_hpd(:,:,kk,n));
            txt = sprintf('roistats.d%d.sl%d.roi.nvox(kk) = sum(tmp_roi(:));',ii,jj);
            eval(txt);
            txt = sprintf('roistats.d%d.sl%d.roi.cntrd(kk) = regionprops(logical(tmp_roi),''Centroid'');',ii,jj);
            eval(txt);
            
            txt = sprintf('tmpdat_v = T1_d%d_sl%d_m(logical(tmp_roi));',ii,jj);
            eval(txt);
            txt = sprintf('roistats.d%d.sl%d.T1.mean(kk) = mean(tmpdat_v);',ii,jj);
            eval(txt);
            txt = sprintf('roistats.d%d.sl%d.T1.stdev(kk) = std(tmpdat_v);',ii,jj);
            eval(txt);
            txt = sprintf('roistats.d%d.sl%d.T1.msr%d = (tmpdat_v);',ii,jj,kk);
            eval(txt);
            
            txt = sprintf('tmpdat_v = T2_d%d_sl%d_m(logical(tmp_roi));',ii,jj);
            eval(txt);
            txt = sprintf('roistats.d%d.sl%d.T2.mean(kk) = mean(tmpdat_v);',ii,jj);
            eval(txt);
            txt = sprintf('roistats.d%d.sl%d.T2.stdev(kk) = std(tmpdat_v);',ii,jj);
            eval(txt);
            txt = sprintf('roistats.d%d.sl%d.T2.msr%d = (tmpdat_v);',ii,jj,kk);
            eval(txt);
            
            txt = sprintf('tmpdat_v = magimg_d%d_sl%d_m(logical(tmp_roi));',ii,jj);
            eval(txt);
            txt = sprintf('roistats.d%d.sl%d.img.mean(kk) = mean(tmpdat_v);',ii,jj);
            eval(txt);
            txt = sprintf('roistats.d%d.sl%d.img.std(kk) = std(tmpdat_v);',ii,jj);
            eval(txt);
            txt = sprintf('roistats.d%d.sl%d.img.snr(kk) = mean(tmpdat_v)/std(tmpdat_v);',ii,jj);
            eval(txt);
            txt = sprintf('roistats.d%d.sl%d.img.msr%d = (tmpdat_v);',ii,jj,kk);
            eval(txt);
            
        end
        n = n+1;
    end
end

%%
% comparison with HPD display results
MRFdisplay_v = [1:4];

for ii = 1:nSl
    figure(200+ii); clf;
    hold on
    for jj = MRFdisplay_v
        
        txt = sprintf('x_v = roistats.HPD.sl%d.T%d.meanTcor;',ii,ii);
        eval(txt);
        txt = sprintf('y_v = roistats.d%d.sl%d.T%d.mean;',jj,ii,ii);
        eval(txt);
        txt = sprintf('err_v = roistats.d%d.sl%d.T%d.stdev;',jj,ii,ii);
        eval(txt);
        errorbar(x_v,y_v,err_v,'x'); xlabel('spec (msec)'); ylabel('MRF (msec)')
        ylim([0 1.2*max(x_v)])
        xlim([0 1.2*max(x_v)])
        txt = sprintf('T%d comparison with HPD spec',ii);
        title(txt)
        grid on;
        legend(HPDlegendtxt(MRFdisplay_v),'Location','northwest')
%         pause()
        
    end
    plot(x_v,x_v)
    HPDlegendtxt2 = [HPDlegendtxt(MRFdisplay_v) {'HPD'}];
    legend(HPDlegendtxt2,'Location','northwest')
    hold off
    grid on
end

save([datadir fn1_stats],'roistats')

t = toc;
disp(['Analyzing MRI system phantom results complete. Elapsed time is ' num2str(t) ' s'])