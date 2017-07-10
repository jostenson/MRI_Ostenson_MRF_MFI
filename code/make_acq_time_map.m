function [ output ] = make_acq_time_map( input )
% make acquisition time map for k-space coords supplied in input

%   INPUT: input.ksp_norm_1 = normalized k-space coords for dimension 1
%              ".ksp_norm_2 = normalized k-space coords for dimension 2 
%              ".tacq = acquisition time of single trajectory interleaf
%              ".effMtx = pos integer is matrix size of k-space

%   OUTPUT: output.acq_time_map = acquisition time map

addpath('/Users/JasonO/Documents/research/MRI_recon/grid3_dct_11aug')
addpath('/Users/JasonO/Documents/research/MRI_recon/sdc3_nrz_11aug')

disp('Creating acquisition time shift map...')
tic;

k1_coords_normalized = input.ksp_norm_1;
k2_coords_normalized = input.ksp_norm_2;
% k3_coords_normalized = input.ksp_norm_3;

% create matrix of acquisition times
tacq = input.tacq;
effMtx = input.effMtx;
[n_rd, n_ph] = size(k1_coords_normalized);
acqtimes_v = (linspace(0,tacq,n_rd))'; % s
acq_times = repmat(acqtimes_v,[1 n_ph]);

crds = zeros(3, n_rd,n_ph);
crds(1,:) = k1_coords_normalized(:);
crds(2,:) = k2_coords_normalized(:);

% roll off kernel
delta = [1.0, 0.0];
k_not = [0.0, 0.0, 0.0];
numThread = 1; % only have 1 data point
ro_kern = grid3_MAT(delta',k_not',[1.0],effMtx,numThread);
% change to complex, shift, then fft
ro_kern = squeeze(ro_kern(1,:,:,:) + 1j*ro_kern(2,:,:,:));
ro_kern = sum(ro_kern,3);
ro_kern = fftshift(ro_kern);

% calculate density compensation
numIter = 25;
osf     = 2.1;
verbose = 1;
DCF = sdc3_MAT(crds,numIter,effMtx,verbose,osf);

numThread = 1;
data = zeros(2, n_rd, n_ph);
data(1,:) = acq_times(:);

grid_times = grid3_MAT(data(:,:),crds(:,:),DCF(:),effMtx,numThread);
grid_times = squeeze(grid_times(1,:,:,:)) + 1i*squeeze(grid_times(1,:,:,:));
grid_times = sum(grid_times,3);
grid_times = fftshift(grid_times);
grid_times = abs(grid_times);

sortgtimes_v = sort(grid_times(:));
maxidx = round(0.999*numel(sortgtimes_v)); % use the 99.9th percentile for max
maxgtimes = sortgtimes_v(maxidx);
acq_time_map = tacq.*grid_times./maxgtimes;
acq_time_map(acq_time_map > tacq) = tacq;

output.acq_time_map = acq_time_map;

t = toc;
disp(['Creating acquisition time shift map complete. Elapsed time is ' num2str(t) ' s.'])

end

