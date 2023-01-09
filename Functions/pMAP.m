function p_MAP = pMAP(amplitude_list, success_list, t, s, p_ab)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Function for all MAP methods, last parameter is function for a priori information

success_list = logical(success_list);
len_a = length(amplitude_list);

if isscalar(t) && isscalar(s)
    p_vec = normcdf(amplitude_list, t, s);                      % [a, 1]

    LH_vec = [p_vec(success_list), 1 - p_vec(~success_list)];
    p_MAP = prod(LH_vec) * p_ab(t, s);	                        % [1, 1]
elseif isvector(t) && isvector(s) && length(t) == length(s)     % t and s are vectors of same length
    [amp_grid, t_grid] = ndgrid(amplitude_list, t);             % [a, t]
    s_grid = repmat(s(:)', [len_a,1]);                          % [a, t]
    p_vec = normcdf(amp_grid, t_grid, s_grid);                  % [a, t]

    LH_vec = cat(1, p_vec(success_list, :), 1 - p_vec(~success_list, :));
    p_MAP = squeeze(prod(LH_vec, 1)).* p_ab(t, s);	            % [1, t]
    % First dimension: amplitude_list
    % Second dimension: (t,s) pair s
else
    if ismatrix(t) && ismatrix(s)                               % t and s are grid matrices of same size [t, s]
        amp_grid = repmat(amplitude_list(:), [1, size(t,1), size(t,2)]); % [a, t, s]
        t_grid = permute(repmat(t, [1, 1, len_a]), [3, 1, 2]);
        s_grid = permute(repmat(s, [1, 1, len_a]), [3, 1, 2]);
    else                                                        % t and s are vectors of different length, including one of them being a scalar
        [amp_grid, t_grid, s_grid] = ndgrid(amplitude_list, t, s);
        [t, s] = ndgrid(t, s);
    end
    p_vec = normcdf(amp_grid, t_grid, s_grid);                  % [a, t, s]
    
    LH_vec = cat(1, p_vec(success_list, :, :), 1 - p_vec(~success_list, :, :));
    p_MAP = squeeze(prod(LH_vec, 1)) .* reshape(p_ab(t, s), size(t));	% [t, s]
    % First dimension: amplitude_list
    % Second and third dimension: t and s
end

end