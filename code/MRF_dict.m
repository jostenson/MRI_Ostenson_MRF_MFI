function output = MRF_dict(input)
% Extended Phase Graph signal generation of MRF signal from Jiang et al., MRM, 2015

%   INPUT: input.nreps = number of frames (or TRs) in MRF sequence
%              ".TI = inversion time
%              ".delk = pos integer is step between states equal to a full
%                   dephasing imparted by crusher gradient
%              ".szomega = pos integer is number of factors of k to 
%                   include in phase history
%              ".T1_v = vector of T1s in dictionary
%              ".T2_v = vector of T2s in dictionary
%              ".B1_v = vector of B1 multipliers in dictionary
%              ".TR_v = vector of repetition times
%              ".TE_v = vector of echo times
%              ".FA_v = vector of flip angles in degrees
%              ".phi_v = vector of RF phases
%              ".reduce = 1 then do SVD compression of dictionary
%              ".sfrac = fraction of dictionary energy to retain if
%                   compressing

%   OUTPUT: output.dict_list = matrix describing T1, T2, B1, etc. for all 
%                       realisitic combinations
%                ".dict_norm = normalized dictionary w/columns
%                       parameterized by dict_list
%                ".dict = dictionary w/columns parameterized by dict_list
%                ".V_red = reduced right singular vectors of dictionary
%                ".dict_red = compressed dictionary
%                ". = all input params copied to output params

disp('Constructing MRF dictionary...')
tic
%% declare parameters
nreps = input.nreps;
TI = input.TI;
delk = 1; % step between states equal to a full dephasing imparted by crusher gradient
szomega = 101; % number of factors of k to include in phase history
T1_v = input.T1_v;
T2_v = input.T2_v;
B1_v = input.B1_v;
TR_v = input.TR_v;
TE_v = input.TE_v;
FA_v = input.FA_v;
phi_v = input.phi_v;

nT2 = numel(T2_v);
nT1 = numel(T1_v);
nB1 = numel(B1_v);


%% dictionary generation

dict = zeros(nreps,nT2*nT1*nB1);
dict_norm = dict;
dict_list = zeros(nT2*nT1*nB1,3);

nn = 1;
for kk = 1:nB1
    B1 = B1_v(kk);
    myFA_v = FA_v.*B1;
    for jj = 1:nT2
        T2 = T2_v(jj);
        for ii = 1:nT1
            
            T1 = T1_v(ii);
            if T1 >= T2
                disp(['MRF dictionary construction: T1 ',num2str(T1),' :: T2 ',...
                    num2str(T2), ' :: B1 ', num2str(B1)]);
                % run sequence using EPG
                tmp_v = EPG_MRF_SSFP( T1,T2,TE_v,TR_v,myFA_v,delk,nreps,szomega,phi_v,TI );
                
                dict(:,nn) = tmp_v;
                dict_norm(:,nn) = tmp_v./sqrt(tmp_v*tmp_v');
                dict_list(nn,:) = [T1, T2, B1];
                nn = nn+1;
            end
            
        end
    end
end

% remove zero rows from listT1T2_m, normFt_m, Ft_m since init matrices are
% too big
dict_list(~any(dict_list,2),:) = [];
dict_norm(:,~any(dict_norm,1)) = [];
dict(:,~any(dict,1)) = [];


% create output structure
output = input;

output.dict_list = dict_list;
output.dict_norm = dict_norm;
output.dict = dict;

%% determine reduced dictionary space

% reduce dimensionality of dictionary via SVD method by McGivney et al,
% IEEE MI, 2014

if input.reduce == 1
    disp('Compressing MRF dictionary by SVD...')
    [~,S,V] = svd(conj(dict_norm'),'econ');
    s_v = diag(S);
    fNRG_v = cumsum(s_v.^2)./sum(s_v.^2);
    nDictSpace = sum(fNRG_v <= input.sfrac);
    
    V_red = V(:,1:nDictSpace);
    dict_red = conj(dict_norm')*V_red;
    
    output.V_red = V_red;
    output.dict_red = dict_red;
    
    disp('Compression of MRF dictionary by SVD complete.')
end


%%

t = toc;
disp(['Construction of MRF dictionary complete. Elapsed time is ' num2str(t) ' s.'])

end