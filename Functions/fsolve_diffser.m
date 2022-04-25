function [x, fval, exitflag, output] = fsolve_diffser(fun, x0, fzeroval, options)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Finds zero of a monotonously increasing function (e.g., AP output)
%
%   Parameters
%   [x, fval, exitflag, output] = fsolve_diffser(fun, x0, xzeroval, options)
%
%   fun:        function handle
%   x0:         estimation of starting value (best maybe near zero)
%   fzeroval:   function value at x = 0      (for AP: -threshold)
%   options:    options set by optimset      (not everything implemented)
%
%   (c) 2021, stefan.goetz@duke.edu, boshuo.wang@duke.edu

%%%%%%%%%%%%%%%%%%%%%%%%%
%   Options
TolX = optimget(options, 'TolX', 1E-6);
TolFun = optimget(options, 'TolFun', 1E-3);
MaxIter = optimget(options, 'MaxIter', 100);
display = not(strcmpi(optimget(options, 'Display', 'off'), 'off'));

%%%%%%%%%%%%%%%%%%%%%%%%%
%   Initialization
x_iter = NaN(1, MaxIter);
f_iter = NaN(1, MaxIter);

i_cnt = 1;                          % First point is (0, f(0)) with f(0)<0
x_iter(i_cnt) = 0;
f_iter(i_cnt) = fzeroval;
x_bounds = [0, NaN];                % Lower and upper bounds of x
f_bounds = [fzeroval, NaN];         % Lower and upper bounds of f

if (display)
    fprintf('I \tx \t\t\tf \t\t\tdx \t\t\tdf\n');
    fprintf('%d \t%1.5f \t%1.5e\n', i_cnt, x_iter(i_cnt), f_iter(i_cnt));
end

i_cnt = 2;                          % Second point is (x0, f(x0))             
x_iter(i_cnt) = x0;
f_iter(i_cnt) = fun(x0);

while (f_iter(i_cnt) < 0) && (i_cnt < MaxIter)  % Upper bound not found
    if (display)
        fprintf('%d \t%1.5f \t%1.5e\n', i_cnt, x_iter(i_cnt), f_iter(i_cnt));
    end
    x_bounds(1) = x_iter(i_cnt);                % Update lower bound
    f_bounds(1) = f_iter(i_cnt);                % Update lower bound
    %   Estimation for first point right of threshold by taking linear
    %   extrapolation (tends to overestimate threshold) 
    %   offset = 0;
    %   Empiriclly: Transition between linear and nonlinear part of threshold
    %   approximately at cutting point between extrapolation of linear onset
    %   and (threshold - 20 mV)
    offset = 20; 
    nextpoint = findnextpoint([0, x_bounds(1)], [fzeroval, f_bounds(1)] + offset);
    nextpoint = max(nextpoint, 1.5 * x_bounds(1));    % To avoid block when threshold higher than at -20 mV
    
    i_cnt = i_cnt + 1;
    x_iter(i_cnt) = nextpoint;
    f_iter(i_cnt) = fun(nextpoint);
    
end

% upper bound found
x_bounds(2) = x_iter(i_cnt);     % Upper bound
f_bounds(2) = f_iter(i_cnt);     % Upper bound

%%%%%%%%%%%%%%%%%%%%%%%%%
%   Iteration
while (i_cnt < MaxIter) && ...              % Maximum iteration not reached
      (diff(x_bounds) > TolX) && ...        % x not bounded within TolX
      (abs(f_iter(i_cnt)) > TolFun)         % f not bounded within TolFun around 0
    
    if (display)
        fprintf('%d \t%1.5f \t%1.5e \t%1.5e \t%1.5e\n',  i_cnt, x_iter(i_cnt), f_iter(i_cnt), diff(x_bounds), diff(f_bounds));
    end
    
    nextpoint = findnextpoint(x_bounds, f_bounds);
    i_cnt = i_cnt + 1;
    x_iter(i_cnt) = nextpoint;
    f_iter(i_cnt) = fun(nextpoint);
    
    if (f_iter(i_cnt) > 0)
        x_bounds(2) = x_iter(i_cnt);     % Update upper bound
        f_bounds(2) = f_iter(i_cnt);     % Update upper bound
    else
        x_bounds(1) = x_iter(i_cnt);     % Update lower bound
        f_bounds(1) = f_iter(i_cnt);     % Update lower bound
    end
end

output.iterations = i_cnt;

if (i_cnt <= MaxIter)
    exitflag = 1;
    output.message = 'Success.';
else
    exitflag = 0;
    output.message = 'Maximum iteration reached.';
    if (diff(x_bounds) > TolX)
        output.message = [output.message, ' Difference in X not within TolX. '];
    end
    if (abs(f_iter(i_cnt)) > TolFun) 
        output.message = [output.message, ' FVal not within TolFun from zero. '];
    end
end

x = x_iter(i_cnt);
fval = f_iter(i_cnt);

if (display)
    fprintf('%d \t%1.5f \t%1.5e \t%1.5e \t%1.5e\n',  i_cnt, x_iter(i_cnt), f_iter(i_cnt), diff(x_bounds), diff(f_bounds));
    figure
    plot(x_iter, f_iter, 'x-')
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%
%   Evaluates zero of a linear function through the two given points
function x_next = findnextpoint(xb, fb)

x_next = xb(1) - (xb(2) - xb(1)) ./ (fb(2) - fb(1)) .* fb(1);

end

