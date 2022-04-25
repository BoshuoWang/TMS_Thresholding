function generateLUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Dedicated initialization for Kontsevich-Tyler, 1999
%   Look-Up table of psychometric function
%   (c) 2011, boshuo.wang@duke.edu, stefan.goetz@duke.edu

dx = 0.01;                          % Stimulation amplitude, precision is 1%  
x_vec = dx : dx : 1.5;              % 1% to 150%
size_x = length(x_vec);

da = 0.0025;                        % Threshold, precision is 0.25%
a_vec = da : da : 1;                % 0.025% to 100%
size_a = length(a_vec);
da_vec = da * ones(size_a,1);

size_b = 250;                       % Slope
tmp_vec = logspace(log10(0.001), log10(0.5), size_b*2+1);	
b_vec = tmp_vec(2:2:end-1);         % Sample at center of intervals
db_vec = diff(tmp_vec(1:2:end));    % Difference is calculated at endpoints of intervals

[x_3Dgrid, a_3Dgrid, b_3Dgrid] = ndgrid(x_vec, a_vec, b_vec);   % [x,a,b]
[a_2Dgrid, b_2Dgrid] = ndgrid(a_vec, b_vec);                    % [a,b]
[da_2Dgrid, db_2Dgrid] = ndgrid(da_vec, db_vec);                % [a,b]

psi_tensor = psychometricfun(x_3Dgrid, a_3Dgrid, b_3Dgrid); 	% [x,a,b]
% phi_tensor = p(suprathreshold | x, (a,b))

% Init. value of pt(a,b), i.e., p0_lambda
% For statistical information:
pt_ab_init = p_ab(a_2Dgrid, b_2Dgrid);                          % [a,b]
% Cleaned up by BW (03/10/2021), provides general shape. Threshold is
% normaly distributed. Slope is centered around 0.25.

% % For all equal:
% % pt_ab_init = ones(1, size_a, size_b);

pt_ab_init = pt_ab_init .* da_2Dgrid .* db_2Dgrid;  % Convert continuous pdf to discrete probabilities
pt_ab_init = pt_ab_init / sum(pt_ab_init, 'all');   % Normalization, total probably should add to 1


save(fullfile('Functions','LUT.mat'),...
		'size_x', 'size_a', 'size_b',...
        'x_vec',  'a_vec',  'b_vec',...
        'psi_tensor', 'pt_ab_init');
end


function psych_fun = psychometricfun(x, t, s)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Psychometric Function f:    r ----> [0, 1]
%                               x ----> y = P(Y > t)
%   Can be Gaussian or as defined in Kontsevich and Tyler, 1999

if (s==0)
    psych_fun = 1e-9 * ones(size(x));
else
    %   Gaussian:
    x = max(min(x, t + 6*s), t - 6*s);      % Bound x to within 6 sigma of mean, so that normcdf is never exact 1!
    psych_fun = normcdf(x, t, s);
end

end

