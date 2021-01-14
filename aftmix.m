function [beta,bb,weight,cons]=aftmix(aft_data,beta,bb,weight,cons)

mmf=@test_linear;

[n,m]=size(aft_data);
XX=[aft_data(aft_data(:,m)==1,m),aft_data(aft_data(:,m)==1,2:(m-1))];
[b,bint,r,rint,stats]=regress(aft_data(aft_data(:,m)==1,1),XX); 

[b,gab]=aftlognormal(aft_data,zeros(m-1,1),iqr(r));
if (nargin==1)
    weight=1; bb=gab/5; beta=b;cons=gab/20;
end
if (nargin==4)
    cons=gab/20; 
end
data=aft_data;

for k=1:30
    [new_sigma,ttt]=find_sigma_univ2(data,beta,weight,bb,mmf,cons);   
    if (max(ttt)<0.01 && k>1) break
    end
    bb=unique([bb,new_sigma]);
    weight=[weight,0*new_sigma];
    weight=ones(1,length(bb))/length(bb);
    [bb,weight,W]=update_weight_univ2(data,bb,weight,beta,mmf,cons);   
    for jj=1:3
        beta=myfun_univ2(data,bb,weight,beta,mmf,cons) ;  
        [bb,weight,W]=update_weight_univ2(data,bb,weight,beta,mmf,cons);                        
    end        
end


