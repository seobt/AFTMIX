function [f,fval,exi]=myfun_univ2(data,bb,weight,beta,meanfunction,cons)
        
nn=size(data,1);
d=size(data,2)-1;
function g=myfun2(theta)            
    me=meanfunction(theta,data(:,2:d));
    gg=zeros(nn,1);
    for j=1:length(weight)                        
        gg=gg+weight(j)*(data(:,d+1).*normal_pdf(data(:,1),me,bb(j)+cons)+(1-data(:,d+1)).*(1-normcdf(data(:,1),me,sqrt(bb(j)+cons))));
    end        
    g=-sum(log(gg));                 
end
bound=10*max(abs(beta))*beta'./beta';
%options=optimset('Algorithm','Active-Set','MaxIter',100,'Display','off');
options = optimoptions(@fmincon,'Display','none','Algorithm','interior-point');
[f,fval,exi]=fmincon(@myfun2,beta,[],[],[],[],-bound,bound,[],options);
end
