%% MRF-MFI batch download code contributors, input data, process, and generate figures

clear, close all, clc

%% download contributors

txt = 'python download_contrib_all.py';
[status, cmdout] = unix(txt)
mex ./contrib/grid3_dct_11aug/grid3_MAT.c -outdir ./contrib/grid3_dct_11aug/
mex ./contrib/sdc3_nrz_11aug/sdc3_MAT.c -outdir ./contrib/sdc3_nrz_11aug/

%% download data input

%% PIQT
clear, close all

% data directory
datadir_in = '../data_in/';
datadir_out = '../data_out/';

fn_ksp_coords = 'MRF_ksp_traj_PIQT.mat';
fn_MRF_data = 'MRF_raw_PIQT.mat';
fn_B0map = 'MRF_B0_PIQT.mat'; %
fn_MRF_dict = 'MRF_dict_PIQT.mat';
fn_proc_noMFI = 'MRF_proc_PIQT.mat';
fn_proc_MFI = 'MRF_proc_PIQT_MFI.mat';
coil_img_name = 'tmpRecon';
basismod = 1.2;
dfreq = 0.1; % Hz
effMtx = 240;
permute_order = [3 2 1];
params.M = 21;

% dictionary params
input_dict.TI = 7.0; % msec
input_dict.T1_v = [20:10:3000 3200:200:5000];
input_dict.T2_v = [10:5:300 350:50:2000];
% input_dict.T1_v = [10:10:90, 100:20:1000, 1040:40:2000, 2050:100:3000];
% input_dict.T2_v = [2:2:8, 10:5:100, 110:10:300, 350:50:800];
input_dict.B1_v = [1.0];

doMagdict= 1;
doTmaps = 1;

load_recon_MRF_withMFI_final

%% Off-resonance shim exp
clear, close all, clc

% data directory
datadir_in = '../data_in/';
datadir_out = '../data_out/';

% reference "well shimmed"

fn_ksp_coords = 'MRF_ksp_traj_shim_exp.mat';
fn_MRF_data = 'MRF_raw_shim_exp_ref.mat';
fn_B0map = 'MRF_B0_shim_exp_ref.mat'; %
fn_MRF_dict = 'MRF_dict_shim_exp.mat';
fn_proc_noMFI = 'MRF_proc_shim_exp_ref.mat';
fn_proc_MFI = 'MRF_proc_shim_exp_ref_MFI.mat';
coil_img_name = 'tmpRecon';
basismod = 1.2;
dfreq = 0.1; % Hz
effMtx = 240;
permute_order = [3 2 1];
params.M = 21;

% dictionary params
input_dict.TI = 7.0; % msec
input_dict.T1_v = [20:10:3000 3200:200:5000];
input_dict.T2_v = [10:5:300 350:50:2000];
% input_dict.T1_v = [10:10:90, 100:20:1000, 1040:40:2000, 2050:100:3000];
% input_dict.T2_v = [2:2:8, 10:5:100, 110:10:300, 350:50:800];
input_dict.B1_v = [1.0];

doMagdict= 1;
doTmaps = 1;

load_recon_MRF_withMFI_final

% y shim

fn_MRF_data = 'MRF_raw_shim_exp_y.mat';
fn_B0map = 'MRF_B0_shim_exp_y.mat'; %
fn_proc_noMFI = 'MRF_proc_shim_exp_y.mat';
fn_proc_MFI = 'MRF_proc_shim_exp_y_MFI.mat';

doMagdict= 0;
doTmaps = 1;

load_recon_MRF_withMFI_final

% z shim

fn_MRF_data = 'MRF_raw_shim_exp_z.mat';
fn_B0map = 'MRF_B0_shim_exp_z.mat'; %
fn_proc_noMFI = 'MRF_proc_shim_exp_z.mat';
fn_proc_MFI = 'MRF_proc_shim_exp_z_MFI.mat';

doMagdict= 0;
doTmaps = 1;

load_recon_MRF_withMFI_final

% yz shim

fn_MRF_data = 'MRF_raw_shim_exp_yz1.mat';
fn_B0map = 'MRF_B0_shim_exp_yz1.mat'; %
fn_proc_noMFI = 'MRF_proc_shim_exp_yz1.mat';
fn_proc_MFI = 'MRF_proc_shim_exp_yz1_MFI.mat';

doMagdict= 0;
doTmaps = 1;

load_recon_MRF_withMFI_final


%% MRI system phantom
clear, close all, clc

% data directory
datadir_in = '../data_in/';
datadir_out = '../data_out/';

% T1 slice "well shimmed"

fn_ksp_coords = 'MRF_ksp_traj_HPD.mat';
fn_MRF_data = 'MRF_raw_HPD_top_shim.mat';
fn_B0map = 'MRF_B0_HPD_top_shim.mat'; %
fn_MRF_dict = 'MRF_dict_HPD.mat';
fn_proc_noMFI = 'MRF_proc_HPD_top_shim.mat';
fn_proc_MFI = 'MRF_proc_HPD_top_shim_MFI.mat';
coil_img_name = 'tmpRecon';
basismod = 1.2;
dfreq = 0.1; % Hz
effMtx = 240;
permute_order = [3 2 1];
params.M = 21;

% dictionary params
input_dict.TI = 7.0; % msec
% input_dict.T1_v = [20:10:3000 3200:200:5000];
% input_dict.T2_v = [10:5:300 350:50:2000];
input_dict.T1_v = [10:10:90, 100:20:1000, 1040:40:2000, 2050:100:3000];
input_dict.T2_v = [2:2:8, 10:5:100, 110:10:300, 350:50:800];
input_dict.B1_v = [1.0];

doMagdict= 1;
doTmaps = 1;

load_recon_MRF_withMFI_final

% T1 slice "poorly shimmed"

fn_MRF_data = 'MRF_raw_HPD_top_B0het.mat';
fn_B0map = 'MRF_B0_HPD_top_B0het.mat'; %
fn_proc_noMFI = 'MRF_proc_HPD_top_B0het.mat';
fn_proc_MFI = 'MRF_proc_HPD_top_B0het_MFI.mat';

doMagdict= 0;
doTmaps = 1;

load_recon_MRF_withMFI_final

% T2 slice "well shimmed"

fn_MRF_data = 'MRF_raw_HPD_mid_shim.mat';
fn_B0map = 'MRF_B0_HPD_mid_shim.mat'; %
fn_proc_noMFI = 'MRF_proc_HPD_mid_shim.mat';
fn_proc_MFI = 'MRF_proc_HPD_mid_shim_MFI.mat';

doMagdict= 0;
doTmaps = 1;

load_recon_MRF_withMFI_final

% T2 slice "poorly shimmed"

fn_MRF_data = 'MRF_raw_HPD_mid_B0het.mat';
fn_B0map = 'MRF_B0_HPD_mid_B0het.mat'; %
fn_proc_noMFI = 'MRF_proc_HPD_mid_B0het.mat';
fn_proc_MFI = 'MRF_proc_HPD_mid_B0het_MFI.mat';

doMagdict= 0;
doTmaps = 1;

load_recon_MRF_withMFI_final

%% in vivo brain

clear, close all, clc

% data directory
datadir_in = '../data_in/';
datadir_out = '../data_out/';

fn_ksp_coords = 'MRF_ksp_traj_invivo_long.mat';
fn_MRF_data = 'MRF_raw_invivo_sl2_long.mat';
fn_B0map = 'MRF_B0_invivo_sl2.mat'; %
fn_MRF_dict = 'MRF_dict_invivo_long.mat';
fn_proc_noMFI = 'MRF_proc_invivo_sl2_long.mat';
fn_proc_MFI = 'MRF_proc_invivo_sl2_long_MFI.mat';
coil_img_name = 'tmpRecon';
basismod = 1.2;
dfreq = 0.1; % Hz
effMtx = 240;
permute_order = [1 2 3];
flip1 = 1;
params.M = 31;

% dictionary params
input_dict.TI = 7.0; % msec
%     input_dict.T1_v = [20:100:3000 3250:250:5000];
%     input_dict.T2_v = [10:20:300 400:100:2000];
input_dict.T1_v = [20:10:3000 3200:200:5000];
input_dict.T2_v = [10:5:300 350:50:2000];
% input_dict.T1_v = [10:10:90, 100:20:1000, 1040:40:2000, 2050:100:3000];
% input_dict.T2_v = [2:2:8, 10:5:100, 110:10:300, 350:50:800];
input_dict.B1_v = [1.0];

doMagdict= 1;
doTmaps = 1;

load_recon_MRF_withMFI_final

%% generate figures

figs_MRI_2017
