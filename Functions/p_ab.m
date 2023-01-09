function pdf_ab = p_ab(t, s)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Probablity distribution of MLE threshold and slope parameters, based on
%   statistical data of virtual population (Subjects_25000.mat)
%   (c) 2021, boshuo.wang@duke.edu, stefan.goetz@duke.edu, 

theta_log_i = [t(:), log10(s(:))];

MU =    [0.656266386641216,   -1.30274896944234   ];
SIGMA = [0.00959309299392127,  0                  ;...
         0,                    0.00157310700025391];
     
pdf_ab = reshape(mvnpdf(theta_log_i, MU, SIGMA), size(t));

end