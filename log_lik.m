function answer=log_lik(data,bb,weight,beta,meanfunction,cons)

    d=size(data,2)-1;
    me=meanfunction(beta,data(:,2:d));
    
    for j=1:length(bb)
        S(:,j)=weight(j)*(data(:,d+1).*normal_pdf(data(:,1),me,bb(j)+cons)+(1-data(:,d+1)).*(1-normcdf(data(:,1),me,sqrt(bb(j)+cons))));
    end    
    answer=sum(log(sum(S,2)));

        