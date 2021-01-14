function [newx,ttt]=find_sigma_univ2(data,beta,weight,variance,meanfunction,cons)

d=size(data,2)-1;
%cons=5/max(size(data));
%cons=0.05;
%cons=0;
%cons=0.001;

%options=optimset('Algorithm','Active-set','MaxIter',1000,'Display','off');
options = optimoptions('fmincon','Display','Off','Algorithm','interior-point');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%penalized_gradient%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function answer=grad(pointer)
    answer=0;
    n=size(data,1);
    tmp=zeros(n,1);
    me=meanfunction(beta,data(:,2:d));
    for j=1:length(weight)
        tmp=tmp+weight(j)*(data(:,d+1).*normal_pdf(data(:,1),me,variance(j)+cons)+(1-data(:,d+1)).*(1-normcdf(data(:,1),me,sqrt(variance(j)+cons))));        
    end    
    answer=sum((data(:,d+1).*normal_pdf(data(:,1),me,pointer+cons)+(1-data(:,d+1)).*(1-normcdf(data(:,1),me,sqrt(pointer+cons))))./max(tmp,1.0e-20*ones(n,1)));  
    answer=n-answer;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%% find new support %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tmp1=data(:,1)-meanfunction(beta,data(:,2:d));
max_range=2*(max(tmp1)-min(tmp1))^2;
min_range=max([0,(min(abs(tmp1))/2)^2]);
x=(min_range:max_range/100:max_range);
x1=zeros(1,15);
for i=1:20
    x1(i)=0.001*2^(-6+i);
end
x=unique(sort([x1,x,variance/2,50,[0.1:0.1:1]]));
newx=[];
ran=[];

for ii=1:length(x)        
    br(ii)=(gradient_sigma_univ2(data,x(ii)+1.0e-8+cons,beta,weight,variance+cons,meanfunction)-gradient_sigma_univ2(data,x(ii)-1.0e-8+cons,beta,weight,variance+cons,meanfunction))/2.0e-8;    
    if (ii>1 && br(ii-1)>0 && br(ii)<0 ) 
        newx=[newx,(x(ii-1)+x(ii))/2];
        ran=[ran;x(ii-1),x(ii)];
    end
end


if (br(length(x))>0) newx=[newx,max_range+1];
    ran=[ran;x(size(x,1),2),max(x)];
end
%%%%%%%%%%%%%%%%%%%%%%%%% Find support within a given interval%%%%%%%
final_can=[];
for t=1:length(newx)        
    [answer,fval,exi]=fmincon(@grad,newx(t),[],[],[],[],ran(t,1),ran(t,2),[],options);
    %if (exi>0 && fval<-0.0000001) final_can=[final_can,answer];
    if (fval<-0.0000001) final_can=[final_can,answer];
    end
end

if (~isempty(final_can)) newx=final_can;
else ttt=0;
    return;
end

if (isempty(newx)) ttt=0;
    return;
end

newx=unique(newx*10000000000000)/10000000000000;

newxx=[];
ttt=[];

for jj=1:length(newx)
    ttt=[ttt,gradient_sigma_univ2(data,newx(jj)+cons,beta,weight,variance+cons,meanfunction)];
    if (gradient_sigma_univ2(data,newx(jj)+cons,beta,weight,variance+cons,meanfunction)>0.0000001) newxx=[newxx,newx(jj)];
    end
end
max(ttt);
if (gradient_sigma_univ2(data,cons,beta,weight,variance+cons,meanfunction)>0) 
    ttt=[ttt,gradient_sigma_univ2(data,cons,beta,weight,variance+cons,meanfunction)];
    newxx=[newxx,0];
end

newx=unique(newxx);   
end