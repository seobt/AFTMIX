function answer=normal_mixture_cen_s(x,d,loc,variance,weight)

%% univariate noraml mixture
answer=0;
for i=1:length(loc)
    answer=answer+weight(i)*(d*normal_pdf(x,loc(i),variance(i))/variance(i)+...  
        (1-d)*normpdf(x,loc(i),variance(i))/sqrt(variance(i)));
end

