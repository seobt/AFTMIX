function pdf = normal_pdf ( x,mu,sigma )

% Univariate normal density with mean mu and standard srror sigma
y=x-mu;
%if (sigma<0.0000001) inv_sigma=1/sigma;
%    pdf = sqrt(inv_sigma)*exp ( -(0.5 * y.^2)*inv_sigma ) / sqrt ( 2.0 * pi);
%    pdf=0.0000001;
%end
%if (sigma>=0.0000001)
    pdf = exp (-y.^2/(2*sigma))  / sqrt ( 2.0 * pi*sigma);
%end
 