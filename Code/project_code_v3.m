clc
clearvars
close all
format long
%% Part A

t=zeros(4,4);

t(1,1)=0.94;
t(1,2)=0.06;
t(2,2)=0.91;
t(2,3)=0.05;
t(2,4)=0.04;
t(3,3)=0.82;
t(3,4)=0.18;
t(4,4)=1.0;

t_d=zeros(5,5);
t_d(1:4,1:4)=t;
t_d(1:4,5)=[3.488583670663872e-04;6.989357626623678e-04;0.0028015033117678;0.00690308614890242];
t_d(5,5)=1.0;

k=1-t_d(:,end);

for i=1:4
    t_d(i,1:end-1)=t_d(i,1:end-1)*k(i);
end

%% Part B

horizon=50;

cs1=zeros(50,5);
cs2=zeros(50,5);
for i=1:horizon
    cs1(i,:)=[1 0 0 0 0]*t_d^(i);
    cs2(i,:)=[0 1 0 0 0]*t_d^(i);
end

p_f=cs1(:,end).*cs2(:,end);

yrs=(1:horizon)';

fig1=figure;
plot(yrs,p_f,'LineWidth',1.5)
xlabel('Year','FontWeight','bold','FontSize',14)
ylabel('System Failure Probability','FontWeight','bold','FontSize',14)
set(gca,'FontSize',14)
box('on');
grid('on');
print(fig1,'system failure_det_vs_yr','-dsvg')
hold off

%% Part C

load fragility_data

for ii=1:4

data_ind=find(fragility_data(:,1)==ii);
dataa=fragility_data(data_ind,:);

IM=dataa(:,2);
EDP=dataa(:,3);

% Provide IM bounds for different condition states
Bound1 = 1;

% Provide condition states numbering
numS=2;
failed = 2.0;
safe=1.0;

% Assign different condition states in X vector based on IM and provided data

    X = zeros(length(IM),1);
    
    for i=1:length(IM)
        if EDP(i)< Bound1
            X(i)=safe;
        else
            X(i)=failed;
        end
    end
    
    k=(0.001:0.001:5.0)';
    %--------------------------------NOMINAL-----------------------------------
    
    [B_nominal,~,stats] = mnrfit(log(IM),X,'model','nominal');
    [pihat_nominal,dlow_nominal,dhi_nominal] = mnrval(B_nominal,log(k),stats,'type','cumulative','model','nominal','confidence',0.95);  
    CIL_nominal = pihat_nominal-dlow_nominal;
    CIH_nominal = pihat_nominal+dhi_nominal;
    
    save(['frag' num2str(ii) '.mat'],'B_nominal','stats');
    
    fig=figure(ii+1);
    hold('on');
    for j = 1:numS-1
        plot((k),1-pihat_nominal(:,j),'-','linewidth',2);
        p = plot((k),1-CIL_nominal(:,j),':','linewidth',1.0);
        color = get(p, 'Color');
        plot((k),1-CIH_nominal(:,j),':','Color',color,'linewidth',1.0)
    end
    xlabel('PGA(g)','FontWeight','bold','FontSize',12);
    ylabel('Probability','FontWeight','bold','FontSize',12);
    set(gca,'FontSize',14)
    title(['Failure Fragility Function for Initial State ' num2str(ii)],'FontWeight','bold','FontSize',13);
    ylim([0 1]);
%     xlim([0 4]);
    box('on');
    grid('on');
    print(fig,['fragility_' num2str(ii)],'-dsvg')
    hold('off');
end

%% Part D

r1=12;
r2=3;

m=7.5;
vs=880;

[lny_m1,lny_sd1]=attenuation(m,r1,vs);
[lny_m2,lny_sd2]=attenuation(m,r2,vs);

samples=1000000;
ns=[samples 1];
pfff1=zeros(4,1);
pfff2=zeros(4,1);


pg1=5*rand(ns);
% pg2=5*rand(ns);

pff1=rand(ns); pff2=rand(ns);
lambda=(randsample([0 1],samples,true,[0.98 0.02]))';
    
for kk=1:4
    load(['frag' num2str(kk) '.mat'])    
    pfff1(kk)=sum(lambda.*normpdf(log(pg1),lny_m1,lny_sd1).*(pff1<=(1-mnrval(B_nominal,log(pg1),stats,'type','cumulative','model','nominal'))))/samples;
    pfff2(kk)=sum(lambda.*normpdf(log(pg1),lny_m2,lny_sd2).*(pff2<=(1-mnrval(B_nominal,log(pg1),stats,'type','cumulative','model','nominal'))))/samples;

end

%% Part E

p_e1=pfff1;
p_e2=pfff2;

pfu1=p_e1+t_d(1:4,end)-p_e1.*t_d(1:4,end);
pfu2=p_e2+t_d(1:4,end)-p_e2.*t_d(1:4,end);

t_u1=zeros(5,5);
t_u2=zeros(5,5);

t_u1(1:end-1,1:end-1)=t_d(1:end-1,1:end-1);
t_u2(1:end-1,1:end-1)=t_d(1:end-1,1:end-1);

t_u1(1:end-1,end)=pfu1; t_u1(end,end)=t_d(end,end);
t_u2(1:end-1,end)=pfu2; t_u2(end,end)=t_d(end,end);


% ko1=(1.0-t_u1(i,end))/sum(t_u1(i,1:end-1));
% ko2=(1.0-t_u2(i,end))/sum(t_u2(i,1:end-1));

for i=1:4
ko1(i)=(1.0-t_u1(i,end))/sum(t_d(i,1:end-1));
    ko2(i)=(1.0-t_u2(i,end))/sum(t_d(i,1:end-1));
    t_u1(i,1:end-1)=t_u1(i,1:end-1)*ko1(i);
    t_u2(i,1:end-1)=t_u2(i,1:end-1)*ko2(i);
end

%% Part F

csu1=zeros(50,5);
csu2=zeros(50,5);
for i=1:horizon
    csu1(i,:)=[1 0 0 0 0]*t_u1^(i);
    csu2(i,:)=[0 1 0 0 0]*t_u2^(i);
end

p_fu=csu1(:,end).*csu2(:,end);

fig2=figure(500);
plot(yrs,p_f,'LineWidth',1.5)
hold on
plot(yrs,p_fu,'r','LineWidth',1.5)

xlabel('Year','FontWeight','bold','FontSize',14)
ylabel('System Failure Probability','FontWeight','bold','FontSize',14)
legend('Deterioration Only','Deterioration/Earthquake','FontSize',14,'Location','northwest')
set(gca,'FontSize',14)
box('on');
grid('on');
print(fig2,'system failure_total_vs_yr','-dsvg')
hold off

%% Part G

% system state matrix
ss=zeros(25,2);
count=0;
for is=1:5
    for js=1:5
        count=count+1;
        ss(count,:)=[is js];
    end
end

tr1=[ones(5,1) zeros(5,4)];
tr2=[ones(5,1) zeros(5,4)];


t_dn=zeros(25,25); % do-nothing transition matrix
t_r1=zeros(25,25); % only repair 1
t_r2=zeros(25,25); % only repair 2
t_r1r2=zeros(25,25);

for it=1:25
    for jt=1:25
    t_dn(it,jt)=t_u1(ss(it,1),ss(jt,1)).*t_u2(ss(it,2),ss(jt,2));
    t_r1(it,jt)=tr1(ss(it,1),ss(jt,1)).*t_u2(ss(it,2),ss(jt,2));
    t_r2(it,jt)=t_u1(ss(it,1),ss(jt,1)).*tr2(ss(it,2),ss(jt,2));
    t_r1r2(it,jt)=tr1(ss(it,1),ss(jt,1)).*tr2(ss(it,2),ss(jt,2));
    end
end


% Dynamic Programming

gamma = 0.95; % discount factor  
% c_s= [0 1 2 5 10 20 30 60 100 150]; % state cost
c_s=zeros(1,25); c_s(25)=50000;
c_a= [0 80 100 140]; % action cost
horizon=50;  % total steps

V = zeros(25,horizon+1);
for i = horizon:-1:1
    for st = 1:25
        Q1 = c_s(st) + c_a(1) + gamma*t_dn(st,:)*V(:,i+1); % Q value of action 1
        Q2 = c_s(st) + c_a(2) + gamma*t_r1(st,:)*V(:,i+1); % Q value of action 2
        Q3 = c_s(st) + c_a(3) + gamma*t_r2(st,:)*V(:,i+1); % Q value of action 3
        Q4 = c_s(st) + c_a(4) + gamma*t_r1r2(st,:)*V(:,i+1); % Q value of action 4
        [V(st,i),policy(st,i)] = min([Q1;Q2;Q3;Q4]); % minQ 
    end
end


%Realization
s(1) = 2; %initial state

pmf=zeros(1,25);
pmf(s(1))=1.0;

% rng(1) % set random seed to reproduce simulation/plots
pf=zeros(horizon,1);

for t=2:horizon
    a(t-1) = policy(s(t-1),t-1);
   % a(t-1) = 1; %uncontrolled
    if a(t-1)==1
        s(t)=randsample(1:25,1,true,t_dn(s(t-1),:));
    elseif a(t-1)==2
        s(t)=randsample(1:25,1,true,t_r1(s(t-1),:));
    elseif a(t-1)==3
        s(t)=randsample(1:25,1,true,t_r2(s(t-1),:));
    else
        s(t)=randsample(1:25,1,true,t_r1r2(s(t-1),:));
    end
    
    a(t-1)=policy(s(t-1),t-1);
    if a(t-1)==1
        pmf=pmf*t_dn;
        pf(t)=pmf(end);
    elseif a(t-1)==2
        pmf=pmf*t_r1;
        pf(t)=pmf(end);
    elseif a(t-1)==3
        pmf=pmf*t_r2;
        pf(t)=pmf(end);
    else
        pmf=pmf*t_r1r2;
        pf(t)=pmf(end);
    end
end

figure (1000)       
plot(a,'-o')
% title('Action history - One realization')
xlabel('Time','FontWeight','bold','FontSize',14)
ylabel('Action','FontWeight','bold','FontSize',14)

figure(2000)
plot(s,'-o')
% title('State history - One realization')
xlabel('Time','FontWeight','bold','FontSize',14)
ylabel('State','FontWeight','bold','FontSize',14)

figure(3000)
heatmap(policy)
colorbar
% title('Optimum Policy')
xlabel('Time')
ylabel('States')

figure(4000)
heatmap(V)
colorbar
% title('Value Function')
xlabel('Time')
ylabel('States')


figure(5000)
plot(1:50,pf,'g','LineWidth',2);
box on
grid on
hold on
plot(yrs,p_f,'b','LineWidth',1.5)
plot(yrs,p_fu,'r','LineWidth',1.5)

xlabel('Year','FontWeight','bold','FontSize',14)
ylabel('System Failure Probability','FontWeight','bold','FontSize',14)
legend('Optimum Policy','Deterioration Only','Deterioration/Earthquake','FontSize',14,'Location','northwest')
set(gca,'FontSize',14)

save('MDP1.mat','a','s','policy','V','pf');

cost=0;
ccost=zeros(50,1);
for i=1:horizon
    cost=cost+V(s(i),i);
    ccost(i)=cost;
end

figure(6000);
plot(yrs,ccost,'b','LineWidth',1.5)
box on
grid on
xlabel('Year','FontWeight','bold','FontSize',14)
ylabel('Optimal Cumulative Cost of the System','FontWeight','bold','FontSize',14)

%% Part H

c_s_cb=zeros(1,25); c_s_cb(25)=50000;
c_a_cb= [0 140]; % action cost
horizon=50;  % total steps

V_cb = zeros(25,horizon+1);
cbm_rs=[13 14 15 18 19 20 23 24 25];
for i = horizon:-1:1
    for st = 1:25
        if ismember(st,cbm_rs)
        Q1 = c_s_cb(st) + c_a_cb(1) + gamma*t_dn(st,:)*V_cb(:,i+1); % Q value of action 1
        Q2 = c_s_cb(st) + c_a_cb(2) + gamma*t_r1r2(st,:)*V_cb(:,i+1); % Q value of action 2
        [V_cb(st,i),policy_cb(st,i)] = min([Q1;Q2]); % minQ
        else
        Q1 = c_s_cb(st) + c_a_cb(1) + gamma*t_dn(st,:)*V_cb(:,i+1); % Q value of action 1
        [V_cb(st,i),policy_cb(st,i)] = min(Q1); % minQ 
        end
    end
end


%Realization
s_cb(1) = 2; %initial state

pmf_cb=zeros(1,25);
pmf_cb(s_cb(1))=1.0;

% rng(1) % set random seed to reproduce simulation/plots
pf_cb=zeros(horizon,1);

for t=2:horizon
    a_cb(t-1) = policy_cb(s_cb(t-1),t-1);
   % a(t-1) = 1; %uncontrolled
    if a_cb(t-1)==1
        s_cb(t)=randsample(1:25,1,true,t_dn(s_cb(t-1),:));
    else
        s_cb(t)=randsample(1:25,1,true,t_r1r2(s_cb(t-1),:)); % repair both when system state is at prescribed states
    end
    a_cb(t)=policy_cb(s_cb(t),t);
    if a_cb(t)==1
        pmf_cb=pmf_cb*t_dn;
        pf_cb(t)=pmf_cb(end);
    else
        pmf_cb=pmf_cb*t_r1r2;
        pf_cb(t)=pmf_cb(end);
    end
end

figure (7000)       
plot(a_cb,'-o')
% title('Action history - One realization')
xlabel('Time','FontWeight','bold','FontSize',14)
ylabel('Action','FontWeight','bold','FontSize',14)

figure(8000)
plot(s_cb,'-o')
% title('State history - One realization')
xlabel('Time','FontWeight','bold','FontSize',14)
ylabel('State','FontWeight','bold','FontSize',14)

figure(9000)
heatmap(policy_cb)
colorbar
title('Optimum Policy')
xlabel('time')
ylabel('states')

figure(9000)
heatmap(V_cb)
colorbar
title('Value Function')
xlabel('time')
ylabel('states')

figure(10000)
plot(1:50,pf,'g','LineWidth',2);
box on
grid on
hold on
plot(yrs,pf_cb,'k','LineWidth',1.5)
plot(yrs,p_f,'b','LineWidth',1.5)
plot(yrs,p_fu,'r','LineWidth',1.5)

xlabel('Year','FontWeight','bold','FontSize',14)
ylabel('System Failure Probability','FontWeight','bold','FontSize',14)
legend('Optimum Policy','CBM Policy','Deterioration Only','Deterioration/Earthquake','FontSize',14,'Location','northwest')
set(gca,'FontSize',14)

save('CBM1.mat','a_cb','s_cb','policy_cb','V_cb','pf_cb')

cost_cb=0;
ccost_cb=zeros(50,1);
for i=1:horizon
    cost_cb=cost_cb+V_cb(s(i),i);
    ccost_cb(i)=cost_cb;
end

