%% load and recon MRF without and with MFI

% note SDC and gridding performed with functions designed for 3D, so usage
% here for 2D may not be optimal
addpath(genpath('./contrib/'))

    
%% load MRF raw data trajectory coords
load([datadir_in fn_MRF_data]);
load([datadir_in fn_ksp_coords]);

%% MRF data


data = output_MRF_raw.data;

if numel(size(data)) == 3 % multicoil
    [nC,n_ph,n_rd] = size(data); %num channels, num acqs, num readout pnts
elseif numel(size(data)) == 2 % single coil
    [n_ph,n_rd] = size(data);
    nC = 1;
else
    disp('Data dimensions unexpected')
    return
end
disp(['num channels ', num2str(nC), ...
    ' :: n phase ' num2str(n_ph) ' :: n read ' num2str(n_rd)])

if numel(size(data)) == 3
    data_comb = permute(data,[3 2 1]);
elseif numel(size(data)) == 2
    data_comb = permute(data,[2 1]);
end

ksp_data = data_comb;

% visualize
for ii = 1:nC
    figure(2); clf;
    imagesc(log(abs(ksp_data(:,:,ii))+1));
    colormap(gray);
    title(['log ksp slice, coil ' num2str(ii)]);
    drawnow
end

%% grid and recon

% k-space coords
crds = zeros(3, n_rd, n_ph);

for ii = 1:2
    txt = sprintf('k%d_coords_normalized = output_ksp_traj.ksp_norm_%d;',ii,permute_order(ii));
    eval(txt);
end

crds(1,:,:) = k1_coords_normalized;
crds(2,:,:) = k2_coords_normalized;

% roll off kernel; gridding using ISMRM Unbound, see also Zwart et al., MRM
% 2012
delta = [1.0, 0.0];
k_not = [0.0, 0.0, 0.0];
numThread = 1;
ro_kern = grid3_MAT(delta',k_not',[1.0],effMtx,numThread);
% change to complex, shift, then fft
ro_kern = squeeze(ro_kern(1,:,:,:) + 1i*ro_kern(2,:,:,:));
ro_kern = sum(ro_kern,3);
ro_kern = flipud(ifftshift(ifftn(ro_kern)));
ro_kern = abs(ro_kern);

% display roll off kernel
figure(3); clf;
imagesc(ro_kern); axis image; colormap(gray);
title('roll off kernel');

% calculate density compensation; SDC using ISMRM Unbound, see also Zwart
% et al., MRM, 2012
numIter = 25;
osf     = 2.1;
verbose = 1;
DCF = sdc3_MAT(crds,numIter,effMtx,verbose,osf);

img_tot = zeros(effMtx,effMtx,nC); % for storing summed coil images
for kk = 1:nC % loop over coils
    
    % initialize for grid recon
    numThread = 1; % or 1 for no pthread lib exec
    data_real = real(ksp_data(:,:,kk));
    data_imag = imag(ksp_data(:,:,kk));
    data = zeros(2, n_rd, n_ph);
    data(1,:) = data_real(:);
    data(2,:) = data_imag(:);
    
    img_Comb = zeros(effMtx,effMtx,n_ph);
    
    for jj = 1:n_ph % loop over MRF acqs
        
        disp(['Processing coil ' num2str(kk) ' :: acq ' num2str(jj)])    
        
        gdata = grid3_MAT(data(:,:,jj),crds(:,:,jj),DCF(:,jj),effMtx,numThread);
       
        
        % change to complex, fft, then shift
        gdata_complex = squeeze(gdata(1,:,:,:) + 1j*gdata(2,:,:,:));   
        gdata_complex = sum(gdata_complex,3);
        gdata_complex = fftshift(gdata_complex);
                
        % reconstruct image
        imgdata_complex = ifftn(gdata_complex);
        imgdata_complex = ifftshift(imgdata_complex);
                
        % apply roll-off kernel
        imgdata_complex(ro_kern~=0) = imgdata_complex(ro_kern~=0)./ro_kern(ro_kern~=0);
        
        % save image measurement
        img_Comb(:,:,jj) = imgdata_complex;
        
    end
    
    % optional flipping of image
    if exist('flip1','var')
        if flip1 == 1
            img_Comb = flip(img_Comb,1);
        end
    end
    if exist('flip2','var')
        if flip2 == 1
            img_Comb = flip(img_Comb,2);
        end
    end
    
    txtsave = sprintf(['save ' [datadir_out coil_img_name] '%d.mat img_Comb'],kk);
    eval(txtsave)
    
    figure(12)
    img_tot(:,:,kk) = sum(img_Comb,3);
    imagesc((abs(img_tot(:,:,kk))));
    axis image
    colormap(gray)
    title(['mag image of the sum of all image acquisitions, coil ' num2str(kk)])
    drawnow
    
end

%% determine noise covariance matrix

coil_noise = output_MRF_raw.noise;
nz_cov = cov(coil_noise');

%% combine coil data

% determine noise correlation matrix
[~, nNz] = size(coil_noise);
nz_cor = (coil_noise*coil_noise')./nNz;

input_AR.nz_cor = nz_cor; %eye(nC);
input_AR.imgtot = img_tot;
input_AR.Nsig = n_ph;
input_AR.var_name = 'img_Comb';
input_AR.coil_img_name = 'tmpRecon';
input_AR.datadir = datadir_out;
input_AR.coil_img_ext = '';
input_AR.plot_yes = 1;

output_AR = MRF_MC_AR(input_AR);

%% construct or load dictionary

TE = output_MRF_raw.params.TE_s;% s
nomFlip = output_MRF_raw.params.nomFlip_deg; % nominal flip angle in degrees
TRbase = output_MRF_raw.params.TRbase_s; % s

if doMagdict == 1
    
    flipCor = 1.0; % flip angle correction mulitiplier
    
    fn_MRF_seq_params = 'MRF.csv';
    % load acquisition vectors
    data_struct = importdata([datadir_in fn_MRF_seq_params]);
    
    % set dictionary construction parameters
    input_dict.FA_v = data_struct.data(:,1)*nomFlip; % deg
    input_dict.phi_v = data_struct.data(:,2); % deg
    input_dict.TR_v = TRbase*1000 + data_struct.data(:,3); % ms
    input_dict.TE_v = TE*1000 + data_struct.data(:,4); % ms
    input_dict.nreps = n_ph;
    
    input_dict.delk = 1; % step between states equal to a full dephasing imparted by crusher gradient
    input_dict.szomega = 101; % number of factors of k to include in phase history

    
    % plot sequence parameters
    figure(1); clf;    
    subplot(411)
    plot(input_dict.FA_v); ylabel('degrees'); title(['FA, TI is ' num2str(input_dict.TI) ' ms'])
    subplot(412)
    plot(input_dict.phi_v); ylabel('degrees'); title('phase')
    subplot(413)
    plot(input_dict.TR_v); ylabel('msec'); title('TR')
    subplot(414)
    plot(input_dict.TE_v); ylabel('msec'); title('TE')
    drawnow
    
    % do dictionary construction
    input_dict.reduce = 1;
    input_dict.sfrac = 0.9999;
    output_dict = MRF_dict(input_dict);
    
    % save result
    output_dict.fn = fn_MRF_dict;
    save([datadir_out fn_MRF_dict],'output_dict')
    
else
    load([datadir_out fn_MRF_dict]);
end


%% do dictionary construction/match


if doTmaps == 1
    
    input_MRF_dp.img_AR_comb = output_AR.img_AR_comb;
    input_MRF_dp.reduce = 1;
    
    output_MRF_match = MRF_dict_match(input_MRF_dp, output_dict);
    
    figure(1000)
    imagesc(output_MRF_match.T1_map)
    axis image
    colorbar()
    title('T1 map')
    
    figure(1001)
    imagesc(output_MRF_match.T2_map)
    axis image
    colorbar()
    title('T2 map')
    
    figure(1002)
    imagesc(abs(output_MRF_match.R_map))
    axis image
    colorbar()
    title('Magnitude of inner product')
    
end

% save dictionary matching results

save([datadir_out fn_proc_noMFI], 'input_MRF_dp','output_MRF_match');

%% MFI section begins

%% construct acquisition time map

tacq = output_MRF_raw.params.tacq_s; % s

input_acq_map.ksp_norm_1 = k1_coords_normalized;
input_acq_map.ksp_norm_2 = k2_coords_normalized;
input_acq_map.ksp_norm_3 = [];
input_acq_map.tacq = tacq;
input_acq_map.effMtx = effMtx;

[ output_acq_map ] = make_acq_time_map( input_acq_map );

acq_time_map = output_acq_map.acq_time_map;

figure(10); clf;
imagesc(fftshift(acq_time_map)); axis image; colorbar();
title('Acquisition time map')

%% Make interpolated B0 map

% load masked B0 map
load([datadir_in fn_B0map]);

B0_map_masked = output_B0.B0_map_masked;

% generate interpolated/extrapolated B0 map
maskidx_v = find(B0_map_masked);
[x_m,y_m] = meshgrid(1:effMtx);
F = scatteredInterpolant(x_m(maskidx_v), y_m(maskidx_v), double(B0_map_masked(maskidx_v)), 'linear', 'nearest');
B0_interp_m = F(x_m, y_m);

% display interpolated map
figure(60); clf;
imagesc(B0_interp_m); axis image; colorbar(); title('interpolated B0 map')


%% establish MFI coeffs
% following Man et al, MRM, 1997
disp('Calculating MFI coefficients...')
tic;

% get median frequency and max abs(freq - median freq) of B0 map
B0_map_MFI = round(B0_interp_m,1);
medB0 = round(median(B0_map_MFI(:)));
B0_map_MFI = B0_map_MFI - medB0;

params.method = 'DFT';
% params.M = 21;
params.t0 = TE;
params.t_end = params.t0 + tacq;
params.t1 = TE;
params.t2 = basismod*(params.t_end) - params.t_end;
params.N = size(k1_coords_normalized,1);
params.f_max = max(abs(B0_map_MFI(:)));
params.f_step = dfreq;

results_MFI_coeff = MFI_coeffs(params);

t = toc;
disp(['Calculating MFI coefficients complete. Elapsed time is ' num2str(t) ' s.'])

%% form coefficient matrix for later use during reconstruction

disp('Forming MFI coefficient matrix for use in recon...')
tic;

% form coefficient map
coeff_map = zeros(effMtx, effMtx, results_MFI_coeff.M);
for ii = 1:effMtx
    for jj = 1:effMtx
        idx = find((abs(results_MFI_coeff.delta_f_fine_v(:) - B0_map_MFI(ii,jj))) < 1e4*eps);
        
        coeff_map(ii,jj,:) = results_MFI_coeff.c_DFT(:,idx);
        
    end
end

t = toc;
disp(['Forming MFI coefficient matrix for use in recon complete. Elapsed time is ' num2str(t) ' s.'])

%% MFI recon using previously saved coil image data

disp('Doing MFI reconstruction...')
tic

img_tot = zeros(effMtx,effMtx,nC);
for ll = 1:nC
    % load given coil image set
    txt = sprintf(['load ' datadir_out coil_img_name '%d.mat'],ll); % load imgComb_m
    eval(txt);
    
    img_MFI = zeros(size(img_Comb));
    
    for mm = 1:n_ph
    
        disp(['MFI B0 adjustment coil ' num2str(ll) ', acq ' num2str(mm)])
        my_ksp = fftn(fftshift(flipud(img_Comb(:,:,mm))));
        my_ksp = repmat(my_ksp,[1 1 results_MFI_coeff.M]);
        t_m = TE + acq_time_map; % time map inclding TE
        my_img = zeros(size(my_ksp));
        
        % determine basis images
        for ii = 1:results_MFI_coeff.M
            my_ksp(:,:,ii) = my_ksp(:,:,ii).*exp(1i*2*pi*t_m*(results_MFI_coeff.delta_f_v(ii)+medB0));
            my_img(:,:,ii) = flipud(ifftshift(ifftn(my_ksp(:,:,ii))));
        end
        
        % determine MFI image
        my_img = my_img.*coeff_map;
        
        img_MFI(:,:,mm) = sum(my_img,3);
        
    end
    
    img_tot(:,:,ll) = sum(img_MFI,3);
    
    txt = sprintf(['save ' [datadir_out coil_img_name] '%d_MFI.mat img_MFI'],ll);
    eval(txt)
    
    figure(101); clf;
    subplot(121)
    imagesc(abs(sum(img_MFI,3))); axis image; colormap(gray);
    title(['Mag image of sum MFI B0 corrected images, coil ' num2str(ll)])
    subplot(122)
    imagesc(angle(sum(img_MFI,3))); axis image; colormap(gray);
    title(['Phase image of sum MFI B0 corrected images, coil ' num2str(ll)])
    drawnow
    
end

t = toc;
disp(['MFI reconstruction complete. Elapsed time is ' num2str(t) ' s.'])

%% combine MFI corrected coil data

input_AR.nz_cor = nz_cor;
input_AR.imgtot = img_tot;
input_AR.Nsig = n_ph;
input_AR.var_name = 'img_MFI';
input_AR.coil_img_name = 'tmpRecon';
input_AR.coil_img_ext = '_MFI';
input_AR.plot_yes = 1;

output_AR = MRF_MC_AR(input_AR);

%% do dictionary match

if doTmaps == 1
    
    input_MRF_dp.img_AR_comb = output_AR.img_AR_comb;
    input_MRF_dp.reduce = 1;
    
    output_MRF_match = MRF_dict_match(input_MRF_dp, output_dict);
    
    figure(1000)
    imagesc(output_MRF_match.T1_map)
    axis image
    colorbar()
    title('T1 map')
    
    figure(1001)
    imagesc(output_MRF_match.T2_map)
    axis image
    colorbar()
    title('T2 map')
    
    figure(1002)
    imagesc(abs(output_MRF_match.R_map))
    axis image
    colorbar()
    title('Magnitude of inner product')
    
end

% save processed results

save([datadir_out fn_proc_MFI], 'input_MRF_dp','output_MRF_match');
save([datadir_out fn_proc_MFI(1:end-4) '_misc.mat'], 'results_MFI_coeff','output_acq_map','B0_map_MFI');

