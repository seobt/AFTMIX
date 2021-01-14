function answer=gradient_sigma_univ2(data,phi,beta,weight,variance,meanfunction)

answer=0;
d=size(data,2)-1;
n=size(data,1);
tmp=zeros(1,n);
tmp=zeros(n,1);
me=meanfunction(beta,data(:,2:d));
for j=1:length(weight)
    tmp=tmp+weight(j)*(data(:,d+1).*normal_pdf(data(:,1),me,variance(j))+(1-data(:,d+1)).*(1-normcdf(data(:,1),me,sqrt(variance(j)))));
end    
tmp=max(tmp,1.0e-10);
    answer=sum((data(:,d+1).*normal_pdf(data(:,1),me,phi)+(1-data(:,d+1)).*(1-normcdf(data(:,1),me,sqrt(phi))))./tmp);       
answer=answer-n;

