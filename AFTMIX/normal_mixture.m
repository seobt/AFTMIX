function answer=normal_mixture(x,d,loc,variance,weight)

%% univariate noraml mixture
answer=0;
for i=1:length(loc)
    answer=answer+weight(i)*(d*normal_pdf(x,loc(i),variance(i))+(1-d)*(1-normcdf(x,loc(i),sqrt(variance(i)))));
end

