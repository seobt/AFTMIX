function [answer1,answer2]=aftlognormal(data,beta,sig)

mmf=@test_linear;
d=size(data,2)-1;

options = optimoptions('fmincon','Display','Off','Algorithm','interior-point');

function answer=lik(param)
    betaa=param(1:d);
    sigm=param(d+1);
    me=mmf(betaa,data(:,2:d));  
    S=data(:,d+1).*normal_pdf(data(:,1),me,sigm)+(1-data(:,d+1)).*(1-normcdf(data(:,1),me,sqrt(sigm)));        
    answer=-sum(log(max(S,1.0e-20)));
end

    [answer,fval,exi]=fmincon(@lik,[beta;sig],[],[],[],[],[-20*ones(1,d),0.0000001],[20*ones(1,d),1000],[],options);
    answer1=answer(1:d);answer2=answer(d+1);
end



