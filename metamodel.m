density=0.01:0.0001:1;
density=density';
m=size(density,1);
complfilename=['complHRr3.txt'];
complfileID=fopen(complfilename);
compliance=textscan(complfileID,'%24.10f');
compliance=reshape(compliance{1,1},3,9901)';
cP1=compliance(:,2);

% complfilenamems=['complHRr3MS4.txt'];
% complfileIDms=fopen(complfilenamems);
% compliancems=textscan(complfileIDms,'%24.10f');
% compliancems=reshape(compliancems{1,1},5,9901)';
% cP1ms=compliancems(:,2);
% 
% complfilenamems11=['complHRr3MS11.txt'];
% complfileIDms11=fopen(complfilenamems11);
% compliancems11=textscan(complfileIDms11,'%24.10f');
% compliancems11=reshape(compliancems11{1,1},5,9901)';
% cP1ms11=compliancems11(:,2);

cP1dec=cP1;
for i=2:m
    cP1dec(i)= min(cP1(i),cP1dec(i-1));
end
cP1decms=cP1ms;
for i=2:m
    cP1decms(i)= min(cP1ms(i),cP1decms(i-1));
end
cP1decms11=cP1ms11;
for i=2:m
    cP1decms11(i)= min(cP1ms11(i),cP1decms11(i-1));
end

dif=max(cP1-cP1dec);
difms=max(cP1ms-cP1decms);
difms11=max(cP1ms11-cP1decms11);

fclose('all')

figure(1)
plot(density,cP1-cP1dec)
figure(2)
plot(density,cP1)
figure(3)
plot(density,cP1dec)
figure(4)
plot(density,cP1)
hold on
plot(density,cP1dec)
hold on
plot(density,cP1ms)
hold on
plot(density,cP1decms)
hold on
plot(density,cP1ms11)
hold on
plot(density,cP1decms11)
legend('penal=1','penal=1 dec','penal=1 multistart4','penal=1 dec multistart4','penal=1 multistart11','penal=1 dec multistart11')
hold off

%1st try : same value in d=0.5 and 1

nV=[0.2,0.5,1,2,5];
aaa5=[];
bbb5=[];
figure(5)
plot(density,cP1decms)
hold on
for i=1:size(nV,2)
    n=nV(i);
    meta1=1./density.^n;
    aaa=(cP1decms(4950)-cP1decms(end))/(meta1(4950)-meta1(end));
    bbb=cP1decms(end)-aaa*meta1(end);
    aaa5=[aaa5, aaa];
    bbb5=[bbb5, bbb];
    metaAdjusted=aaa*meta1+bbb;
    plot(density,metaAdjusted)
    hold on
end
legend('raw','0.2','0.5','1','2','5')
hold off

aaa6=[];
bbb6=[];
figure(6)
plot(density,cP1decms)
hold on
for i=1:size(nV,2)
    n=nV(i);
    meta1=1./density.*exp((density.^n)./n);
    aaa=(cP1decms(4950)-cP1decms(end))/(meta1(4950)-meta1(end));
    bbb=cP1decms(end)-aaa*meta1(end);
    aaa6=[aaa6, aaa];
    bbb6=[bbb6, bbb];
    metaAdjusted=aaa*meta1+bbb;
    plot(density,metaAdjusted)
    hold on
end
legend('raw','0.2','0.5','1','2','5')
hold off

% figure(7)
% plot(density,cP1decms)
% hold on
% for i=1:size(nV,2)
%     n=nV(i);
%     meta1=1./density.*exp((density.^n));
%     aaa=(cP1decms(4950)-cP1decms(end))/(meta1(4950)-meta1(end));
%     bbb=cP1decms(end)-aaa*meta1(end);
%     metaAdjusted=aaa*meta1+bbb;
%     plot(density,metaAdjusted)
%     hold on
% end
% legend('raw','0.2','0.5','1','2','5')
% hold off

aaa5
bbb5
aaa6
bbb6


%%relative difference

figure(8)
for i=1:size(nV,2)
    n=nV(i);
    meta1=1./density.^n;
    aaa=(cP1decms(4950)-cP1decms(end))/(meta1(4950)-meta1(end));
    bbb=cP1decms(end)-aaa*meta1(end);
    metaAdjusted=aaa*meta1+bbb;
    plot(density,(metaAdjusted-cP1decms)./cP1decms)
    hold on
end
legend('0.2','0.5','1','2','5')
hold off

figure(9)
for i=1:size(nV,2)
    n=nV(i);
    meta1=1./density.*exp((density.^n)./n);
    aaa=(cP1decms(4950)-cP1decms(end))/(meta1(4950)-meta1(end));
    bbb=cP1decms(end)-aaa*meta1(end);
    metaAdjusted=aaa*meta1+bbb;
    plot(density,(metaAdjusted-cP1decms)./cP1decms)
    hold on
end
legend('0.2','0.5','1','2','5')
hold off

%% no longer affine but linear

figure(10)
for i=1:size(nV,2)
    n=nV(i);
    meta1=1./density.^n;
    aaa=cP1decms(end)/meta1(end);
    metaAdjusted=aaa*meta1;
    plot(density,(metaAdjusted-cP1decms)./cP1decms)
    hold on
end
legend('0.2','0.5','1','2','5')
hold off

figure(11)
for i=1:size(nV,2)
    n=nV(i);
    meta1=1./density.*exp((density.^n)./n);
    aaa=cP1decms(end)/meta1(end);
    metaAdjusted=aaa*meta1;
    plot(density,(metaAdjusted-cP1decms)./cP1decms)
    hold on
end
legend('0.2','0.5','1','2','5')
hold off


%%no longer every n but find through 1 point (or more?)

nlow=0.01;
nup=100;
di=4950 % density coordinate of intersection between data and model
datadi=cP1decms(di);
while nup-nlow>0.01
    n=(nup+nlow)/2;
    model=1./density.*exp((density.^n)./n);
    aaa=cP1decms(end)/model(end);
    modelAdjusted=aaa*model;
    modeldi=modelAdjusted(di);
    if modeldi > datadi
        nup=n;
    else
        nlow=n;
    end
end
n
figure(12)
plot(density,(modelAdjusted-cP1decms)./cP1decms)


%%see influence of point position
maxDif=[];
allN=[];
for di=1:10:9900
    nlow=0.01;
    nup=100;
    %di=4950 % density coordinate of intersection between data and model
    datadi=cP1decms(di);
    while nup-nlow>0.01
        n=(nup+nlow)/2;
        model=1./density.*exp((density.^n)./n);
        aaa=cP1decms(end)/model(end);
        modelAdjusted=aaa*model;
        modeldi=modelAdjusted(di);
        if modeldi > datadi
            nup=n;
        else
            nlow=n;
        end
    end
    maxdif=max(abs((modelAdjusted-cP1decms)./cP1decms));
    maxDif=[maxDif, maxdif];
    allN=[allN, n];
end
figure(13)
plot(1:10:9900,maxDif)
figure(14)
plot(1:10:9900,allN)
