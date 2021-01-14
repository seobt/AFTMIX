%%%%%%%%% PBC data fitting
load pbc

pbc(:,1)=log(pbc(:,1)/365.25);
pbc(:,5:7)=log(pbc(:,5:7));
pbc(:,2)=[pbc(:,2)==2];
[n,m]=size(pbc);
pbc_data=[pbc(:,1),pbc(:,3:7),pbc(:,2)];
aft_data=pbc_data;
[b,gab]=aftlognormal(aft_data,zeros(m-1,1),0.1);
[beta,bb,weight,cons]=aftmix(pbc_data);



boot_beta=[];
rng(20210101);
for boot=1:500
    boot
    boot_data=aft_data(randsample(size(aft_data,1),size(aft_data,1),true),:);
    [boot_b,a1,a2,a3]=aftmix(boot_data,beta,bb,weight);
    boot_beta=[boot_beta;boot_b'];                       
end

sqrt(var(boot_beta))

% Draw density plots
figure
hold on
mixture_graph(0,gab,1,'-k',2)
mixture_graph(0*bb,bb+cons,weight,'-r',2)
legend('PAFT','MAFT')
title('Estimated error distribution for PBC data')
hold off

sqrt(var(boot_beta))

%%%%%%%%% NWTCO data fitting
load nwtco_age

aft_data=bbb;
[beta,bb,weight,cons]=aftmix(aft_data);
[n,m]=size(aft_data);
[b,gab]=aftlognormal(aft_data,beta,0.1);

boot_beta=[];
rng(20210101);
for boot=1:500
    boot
    boot_data=aft_data(randsample(size(aft_data,1),size(aft_data,1),true),:);
    %[boot_b,a1,a2,a3]=aftmix(boot_data,beta,bb,weight,cons);
    [boot_b,a1,a2,a3]=aftmix(boot_data,beta,bb,weight,cons);
    boot_beta=[boot_beta;boot_b'];                       
end
sqrt(var(boot_beta))

% Draw density plots
figure
hold on
mixture_graph(0*bb,bb+cons,weight,'-r',2)
mixture_graph(0,gab,1,'-k',2)
legend('PAFT','MAFT')
title('Estimated error distribution for NWTCO data')
hold off



