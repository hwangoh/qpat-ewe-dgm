function y=GenPrmtrsGaussian(x,Mu,Gam)

Mu=Mu(:); n=length(Mu);
if ~(size(x,1)==n)||~(numel(Gam)==n^2)
    error('gaussian:size',['Dimensions dont agree. Ensure that number of'...
        ' rows of x is the same as the dimension of Mu. Or make sure and Gam is square.'])
end
if ~(rank(Gam)==n), error('gaussian:inv','Gam is not invertible.'), end

% 
xc=bsxfun(@plus,-Mu,x);
y=(2*pi)^(-n/2)*(det(Gam))^(-1/2)*exp(-0.5*sum(xc.*(Gam\xc),1));