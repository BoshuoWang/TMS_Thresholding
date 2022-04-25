function pdf_ab = p_ab(t, s)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Additional information, based on statistical data
%   (c) 2021, boshuo.wang@duke.edu

theta_log_i = [t(:), log10(s(:))];

MU =    [0.656266386641216, -1.30274896944234];
SIGMA = [ 0.00959309299392127,  0                   ; ...
          0,                    0.00157310700025391 ];

pdf_ab = mvnpdf(theta_log_i, MU, SIGMA);
pdf_ab = reshape(pdf_ab, size(t));

end


