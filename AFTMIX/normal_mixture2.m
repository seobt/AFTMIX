function answer=normal_mixture2(x,d,loc,variance,weight)

%% univariate noraml mixture
n=size(x,1);
answer=zeros(n,1);
for i=1:length(variance)
    answer=answer+weight(i)*(d.*normal_pdf(x,loc,variance(i))+(1-d).*(1-normcdf(x,loc,sqrt(variance(i)))));
end

