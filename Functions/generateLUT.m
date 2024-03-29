function generateLUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Dedicated initialization for Kontsevich-Tyler, 1999
%   Look-Up table of psychometric function
%   (c) 2011, stefan.goetz@duke.edu
%   (c) 2021, boshuo.wang@duke.edu

% Stimulation amplitude cannot exceed 100% MSO
dx = 0.01;                              % Stimulation amplitude, precision is 1% MSO  
x_min = dx; x_max = 1.3;                % 10% MSO to 130% MSO
x_vec = x_min : dx : x_max;         
size_x = length(x_vec);

% Threshold parameter can be larger than 100% MSO
da = 0.002;                             % Threshold, precision is 0.2% MSO
a_min = da; a_max = 1.3;                % 0.2% MSO to 130% MSO
a_vec = a_min : da : a_max;        
size_a = length(a_vec);
da_vec = da * ones(size_a,1);

b_min = 0.005;     b_max = 0.5;         % Two orders-of-magnitude, one above/below median (0.05)
size_b = round(log10(b_max) - log10(b_min)) * 50;     % 50 points per decade
tmp_vec = logspace(log10(b_min), log10(b_max), size_b*2+1);	
b_vec = tmp_vec(2:2:end-1);             % Sample at center of intervals
db_vec = diff(tmp_vec(1:2:end));        % Difference is calculated between endpoints of intervals

[x_3Dgrid, a_3Dgrid, b_3Dgrid] = ndgrid(x_vec, a_vec, b_vec);   % [x,a,b]
[a_2Dgrid, b_2Dgrid] = ndgrid(a_vec, b_vec);                    % [a,b]
[da_2Dgrid, db_2Dgrid] = ndgrid(da_vec, db_vec);                % [a,b]

psi_tensor = psychometricfun(x_3Dgrid, a_3Dgrid, b_3Dgrid); 	% [x,a,b]
% phi_tensor = p(suprathreshold | x, (a,b))

% Init. value of pt(a,b), i.e., p0_lambda
% For statistical information:
pt_ab_init = p_ab(a_2Dgrid, b_2Dgrid);                          % [a,b]
pt_ab_init = pt_ab_init .* da_2Dgrid .* db_2Dgrid;              % Convert continuous pdf to discrete probabilities
pt_ab_init = pt_ab_init / sum(pt_ab_init, 'all');               % Normalization, total probably should add to 1

% For uniform initial distribution:
pt_ab_uniform = ones(size_a, size_b);
pt_ab_uniform = pt_ab_uniform .* da_2Dgrid .* db_2Dgrid;        % Convert continuous pdf to discrete probabilities
pt_ab_uniform = pt_ab_uniform / sum(pt_ab_init, 'all');         % Normalization, total probably should add to 1


filepath = fileparts([pwd,filesep]);
if ~contains(filepath,'Functions')
    filename = fullfile('Functions','LUT.mat');
else
    filename = 'LUT.mat';
end
save(filename,...
    'size_x', 'size_a', 'size_b',...
    'x_vec',  'a_vec',  'b_vec',...
    'psi_tensor', 'pt_ab_init', 'pt_ab_uniform');
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