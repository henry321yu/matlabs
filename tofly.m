clear all
close all

% mo=load('C:\Users\henry chen\Desktop\magtest\proj\10 12 中油 2set compare\txt data\2m垂直\2m垂直-2.txt');
 mo=load('all\5m平行-6.txt');

[No,~]=size(mo); 
t0=[1 No];
smk=25; %移動平均係數
smk=1; %移動平均係數
testn='1'; %檔名
t=mo(:,1);

mo_cor=mo;


acc1 = mo_cor(:,5:7);
% acc2 = mo_cor(:,5:7);
m= mo_cor(:,2:4);

ax=acc1(:,1);ay=acc1(:,2);az=acc1(:,3);
% ax2=acc2(:,1);ay2=acc2(:,2);az2=acc2(:,3);
aGo=sqrt(mo_cor(:,2).^2+mo_cor(:,3).^2+mo_cor(:,4).^2);
% aG2o=sqrt(mo_cor(:,5).^2+mo_cor(:,6).^2+mo_cor(:,7).^2);
aG=sqrt(ax.^2+ay.^2+az.^2);
% aG2=sqrt(ax2.^2+ay2.^2+az2.^2);

GG=sqrt(mo_cor(:,8).^2+mo_cor(:,9).^2+mo_cor(:,10).^2)*0.01;

mm4=sqrt(m(:,1).^2+m(:,2).^2+m(:,3).^2);

x=smooth(m(:,1),smk);
y=smooth(m(:,2),smk);
z=smooth(m(:,3),smk);
smedm=[x y z];
xyz=sqrt(x.^2+y.^2+z.^2);
plot(x)
%%
clc
close all;
t0=t-t(1);
number=1; %
% number=8428; %175 2m垂直
% number=6123; %105
% number=7397; %69
% number=8085; %49
% number=7909; %100
% number=8542; %71

% number=1170; %63 4m垂直
% number=1345; %11 & 214
% number=2446; %220
% number=1942; %170
% number=3473; %127
% number=1666; %105

% number=10265; %30 4m平行
% number=8083; %50
% number=7591; %6 &206
% number=8198; %184
% number=10344; %97
% number=7058; %1 &135

% number=7675; %168 5m平行
% number=1; %1 &199
% number=7967; %161
% number=6314; %70
% number=6484; %47
number=6610; %146

tt=[];tt0=[];fixall=[];
cycletime=0.555;
temp1=floor(t0(number)/cycletime);
[~,fix]=min(abs(t0-(t0(number)-temp1*cycletime)));

waven=floor(t0(end)/cycletime);
for i=1:waven
    [~,fixall(i)]=min(abs(t0-(t0(fix)+(cycletime*i-cycletime))));
end
tt=fixall;
tt0=abs(tt-15);

t0(tt(1));
tt(1)
tt(temp1+1);

hold on
plot(x)

[~,k]=size(tt);
for i=1:k
    plot(tt0(i),x(tt0(i)),'b.')
    plot(tt(i),x(tt(i)),'r.')
end

ticks=[tt];
% ticks=[tt(temp1)];
% ticks=[tt tt0];
ticks=reshape (ticks, 1, numel(ticks));
ticks=sort(ticks,'ascend');
% set(gca, 'xtick', ticks);
% grid on
% xlabel('Time')
% ylabel('Gauss')
hold off

