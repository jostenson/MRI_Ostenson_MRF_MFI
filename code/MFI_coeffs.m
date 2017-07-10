
function results = MFI_coeffs(params)

% returns MFI coefficients and associated parameters of interest

%   INPUT: params.method = 'DFT';
%               ".M = pos integer is number of basis frequencie
%               ".params.t0 = start of acq time, TE for center out
%               ".t_end = end of acquisiiton time, TE + tacq for center out
%               ".t1 = time extension on the beginning of target vectors
%               ".t2 = time extension on the end of target vectors
%               ".N = pos integer is discretizations of time
%               ".f_max = maximum frequency to consider
%               ".f_step = frequency spacing in lookup table to be
%                   generated
%   OUTPUT: results.del_f = see above
%                 ".delta_f_v = the frequency basis
%                 ".delta_f_fine_v = the lookup table freq. entries
%                 ".t_v = the total time vector
%                 ".M = see above
%                 ".L = M - 1
%                 ".n_fine = number in delta_f_fine_v
%                 ".A = system matrix
%                 ".A_cond = condition number of A
%                 ".c_DFT = the MFI coefficients
%                 ".acq_b = the target/solution vector during the
%                   acquisition window
%                 ".acq_approx_DFT = the approximation of the solution
%                   using the derived coefficients
%                 ".acq_resid_DFT = the residuals of the fits
%                 ".acq_t_v = the time vector during acquisition
%                 ".b_all = the total acquisition (currentlly commented
%                   out)

%% sign of complex exponentials
i_sign = -1;

switch params.method
    case 'DFT'
        N = params.N;
        L = params.M - 1;

        %% time endpoints and extension
        t1 = params.t1;
        t2 = params.t2;
        
        %% time vector
        T = (params.t_end + t2) - (params.t0 - t1);
        dt = T/N;
        t_v = ( 0:(N-1) )*dt;
       
        %% basis frequencies vector
        del_f = 1/T;
        if mod(params.M,2)==0,
            delta_f_v = [-params.M/2:(+params.M/2-1)] * del_f;
        else
            delta_f_v = [-(params.M-1)/2:(+(params.M-1)/2)] * del_f;
        end
        
        %% system matrix
        
        A = exp(i_sign * 1i * 2 * pi * t_v(:) * delta_f_v);
        A = A./sqrt(N);
%         pA = pinv(A); %
        
        %% calculate coefficients over fine set of frequencies
        delta_f_fine_v = [ round(-abs(params.f_max),-log10(params.f_step)):...
            params.f_step :...
            round(+abs(params.f_max),-log10(params.f_step))];
        n_fine = numel(delta_f_fine_v);
        
        %% coefficients
        c_DFT = complex(zeros(params.M,n_fine));
        
        %% approximations
        acq_idx = find( (t_v>=params.t0) & (t_v<=(params.t_end)) );
        N_acq = numel(acq_idx);
        acq_approx_DFT = complex(zeros(N_acq, n_fine));
        acq_b = complex(zeros(N_acq,n_fine));
        b_taper_all = complex(zeros(N,n_fine));
        
        
        %% loop over fine frequencies
        tic;
        for idx=1:n_fine,
            
            %% ideal phasor
            b_v = exp(i_sign * 1i * 2 * pi * delta_f_fine_v(idx) * t_v(:) )./sqrt(N);
            
            %% taper b
            b_taper = b_v;
            if isfield(params,'taper')
                if strcmp(params.taper,'off')
                    % do nothing
                else
                    bb = ( b_v(1) + b_v(end) ) / 2;
                    for idx2 = 1:acq_idx(1),
                        b_taper(idx2) = bb + (b_v(idx2) - bb) * sin( pi * ( t_v(idx2)/t1 )/2).^2;
                    end
                    for idx2 = acq_idx(end):N,
                        b_taper(idx2) = bb + (b_v(idx2) - bb) * sin( pi * ( (t_v(end)-t_v(idx2))/t2 )/2).^2;
                    end
                end
            else
                bb = ( b_v(1) + b_v(end) ) / 2;
                for idx2 = 1:acq_idx(1),
                    b_taper(idx2) = bb + (b_v(idx2) - bb) * sin( pi * ( t_v(idx2)/t1 )/2).^2;
                end
                for idx2 = acq_idx(end):N,
                    b_taper(idx2) = bb + (b_v(idx2) - bb) * sin( pi * ( (t_v(end)-t_v(idx2))/t2 )/2).^2;
                end
            end
            
            b_taper_all(:,idx) = b_taper;
            
            %% Solve linear system
%             c_DFT(:,idx) = pA*b_taper(:);
            c_DFT(:,idx) = A'*b_taper(:);
            
            %% store ideal
            acq_b(:,idx) = b_v(acq_idx);
            
            %% store approximations
            tmp_v = A * c_DFT(:,idx);
            acq_approx_DFT(:,idx) = tmp_v(acq_idx);
            
            %% status update
%             if mod(idx,100)==0,
%                 disp( sprintf('COMPLETED %04d OF %04f : %.1f SECONDS ELAPSED', idx, n_fine, toc) );
%             end
            
        end
        
        %% assign results
        results = [];
        results.del_f = del_f;
        results.delta_f_v = delta_f_v;
        results.delta_f_fine_v = delta_f_fine_v;
        results.t_v = t_v;
        results.M = params.M;
        results.L = L;
        results.n_fine = n_fine;
        results.A = A;
        results.A_cond = cond(A);
        results.c_DFT = c_DFT;
        results.acq_b = acq_b;
        results.acq_approx_DFT = acq_approx_DFT;
        results.acq_resid_DFT = acq_b - acq_approx_DFT;
        results.acq_t_v = t_v(acq_idx);
%         results.b_all = b_taper_all;
    
end
end


