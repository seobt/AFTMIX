function [bb,weight,W]=update_weight_univ2(data,bb,weight,beta,meanfunction,cons)

n=size(data,1);
d=size(data,2)-1;
m=length(bb);
z=zeros(n,m);
old_weight=weight;
S=zeros(n,m);
W=[];

for itr=1:40     
    itr;
    me=meanfunction(beta,data(:,2:d));
    for j=1:m
        S(:,j)=(data(:,d+1).*normal_pdf(data(:,1),me,bb(j)+cons)+(1-data(:,d+1)).*(1-normcdf(data(:,1),me,sqrt(bb(j)+cons))))./normal_mixture2(data(:,1),data(:,d+1),me,bb+cons,weight);
    end            
    %options = optimoptions('lsqlin','Display','Off','Algorithm','interior-point');
    %eps=1.0e-14;       
    %weight=lsqlin(S,2*ones(n,1),-diag(ones(m,1)),zeros(1,m),ones(1,m),1,eps*ones(1,m),(1-eps)*ones(1,m),[],options);   

    gamma=sqrt(n*10^(-6));
    Sc=[gamma*S;ones(1,m)];
    onec=2*gamma*ones(n,1);
    onec=[onec;1];
    options = optimset('TolX',eps*norm(Sc,1)*length(Sc));
    %options = optimset('TolX',1.0e-10);    
    weight=lsqnonneg(Sc,onec,options);
    weight=weight/sum(weight);
    weight=weight'; 
    %%%%%%%% Armijo rule %%%%%%%%%%%%%%%%% Not so different...
    eta=weight-old_weight;
    sig=1/2;alpha=0.3;
    old_lik=log_lik(data,bb,old_weight,beta,meanfunction,cons);
    for k=0:30        
        test_weight=abs(old_weight+sig^k*eta);
        test_weight=test_weight/sum(test_weight);
        new_lik=log_lik(data,bb,test_weight,beta,meanfunction,cons);
        if (new_lik>=old_lik+ones(1,n)*S*eta'*sig^k*alpha) 
            break
        end
    end   
    weight=test_weight;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    if (norm(new_lik-old_lik)<1.0e-12) 
        break;
    end
    old_weight=weight;
end

bb(weight<=1.1e-14)=[];
weight(weight<=1.1e-14)=[];