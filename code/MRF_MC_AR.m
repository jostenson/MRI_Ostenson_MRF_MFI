function output = MRF_MC_AR(input)

%% optimally combine multi-channel data
% based on Walsh et al., MRM, 2000
% for use in conjunction with MRF processing

%   INPUT: input.nz_cor = noise correlation matrix
%              ".imgtot = Nr x Nc x nC matrix is complex sum or mean of all
%                   MRF frames for each coil
%              ".Nsig = number of MRF frames/TRs
%                   temporary .mat recon'd MRF files must be saved for each coil
%                   input.var_name = string is variable name of MRF dataset in .mat
%                   file
%              ".coil_img_name = string is base name of temporarily saved
%                   MRF uncombined coil sets
%              ".datadir = string is name of directory where temporary data
%                   is stored
%              ".coil_img_ext = string is name of extension of base
%              ".plot_yes = 1 then plot magnitude of complex
%                   sum of coil-combined  MRF set
%   OUTPUT: output.img_AR_comb = Nr x Nc x Nsig coild-combined MRF set

disp('Optimally combining coil data...');
tic;

[Nr, Nc, nC] = size(input.imgtot);
Nsig = input.Nsig;

if nC > 1
    
    Rn_m = input.nz_cor;
    % determine matched filter by location
    coeffComb = zeros(Nr,Nc,nC);
    for ii = 1:Nr
        for jj = 1:Nc
            msr_v = squeeze(input.imgtot(ii,jj,:)); % across coils at given pixel
            Rs_m = msr_v*msr_v'; % the signal correlation matrix
            P_m = Rn_m\Rs_m; % the matrix product
            [~,D_m] = eig(P_m);
            d_v = diag(D_m);
            [~,maxEigInd] = max(d_v(:));
            m_v = P_m(:,maxEigInd);
            m_v = m_v./sqrt(m_v'*(Rn_m\m_v));
            coeffComb(ii,jj,:) = m_v;
        end
    end
    
    img_AR_comb = zeros(Nr,Nc,Nsig);
    % combine MRF acquisitions
    for ii = 1:nC
        txtload = sprintf(['load ' input.datadir input.coil_img_name '%d' input.coil_img_ext '.mat'],ii);
        eval(txtload);
%         if strcmp(input.coil_img_ext,'_MFI')
%             imgComb_m = imgMFI_m;
%         end
        txt = sprintf(['img_Comb = ' input.var_name ';']);
        eval(txt);
        c_v = coeffComb(:,:,ii);
        c_v = c_v(:);
        for jj = 1:Nsig
            imgplane_v = img_Comb(:,:,jj); % get the data for the acq
            imgplane_v = imgplane_v(:);
            imgplane_v = conj(c_v).*imgplane_v; % apply the conj(coeff)
            imgplane = reshape(imgplane_v,[Nr Nc]);
            img_AR_comb(:,:,jj) = img_AR_comb(:,:,jj) + imgplane; % add coil data to current img
        end
    end
    
else
    
    txt = sprintf(['img_Comb = ' input.var_name ';']);
    eval(txt);
    img_AR_comb = imgComb;
    
end

output.img_AR_comb = img_AR_comb;

t = toc;
disp(['Coil combination complete, elapsed time is ' num2str(t) ' s.'])

if input.plot_yes == 1

    imgCombTot = abs(sum(img_AR_comb,3));
    
    figure(1); clf;
    imagesc(imgCombTot); axis image; colormap(gray);
    title('coil-combined magnitude image, sum of all acquisitions')

end

end