% Proto moving horizon estimator for mass spring damper system
% Author: Adrian-Josue Guel-Cortez 2024
% In this case, we know m but we do not know k and b.
% The code uses
% fmincon from MATLAB

close all;
clearvars;
%% measurements
Ts=5e-2; tf=50;
t=0:Ts:tf;
rng(0);
F=zeros(1,length(t));
%F=0.1*sin(t);
w=3e-4;
k=1.9; b=0.5; m=1;
y=zeros(2,length(t));
y(:,1)=[1;0.1];
for j=2:length(t)
    y(:,j)=y(:,j-1)+Ts*model(y(:,j-1),F(j-1),m,b,k)+w*rand(2,1);
end

%% optimiser
lb = [-1.1,-1.1,1.5,0];
ub = [1.1,1.1,2.1,0.8];
A = [];
b = [];
Aeq = [];
beq = [];
nonlcon = [];
u0 = 0;
input=zeros(1,length(t));
state=zeros(4,length(t)); %state=[x,\dotx,k,b]
options = optimoptions('fmincon','Algorithm','interior-point');
input(1)=u0;

figure
subplot(2,2,1)
plot(t,y(1,:),'b')
ylabel('$x$','Interpreter','latex',FontSize=12)
hold on
subplot(2,2,2)
plot(t,y(2,:),'b')
ylabel('$\dot{x}$','Interpreter','latex',FontSize=12)
hold on
subplot(2,2,3)
yline(1.9,'r')
ylabel('$k$','Interpreter','latex',FontSize=12)
hold on
subplot(2,2,4)
yline(0.5,'r')
ylabel('$b$','Interpreter','latex',FontSize=12)
hold on
N=30; %window size
prev_x0=[1,0.1,0.5,0.1]';
for k=N+1:length(t)
    new_x0= fmincon(@(x) MHE(Ts,x,input(k-(N-1):k),y(:,k-(N-1):k)),prev_x0',A,b,Aeq,beq,lb,ub,nonlcon,options);
    [cost,xfinal]=MHE(Ts,new_x0,input(k-(N-1):k),y(:,k-(N-1):k));
    state(1:2,k)=xfinal(1:2)+Ts*model(xfinal(1:2),input(k),1,xfinal(4),xfinal(3));
    state(3:4,k)=xfinal(3:4);
    prev_x0=new_x0;
    subplot(2,2,1)
    plot([k*Ts,(k-1)*Ts],[state(1,k),state(1,k-1)],'k-o');
    subplot(2,2,2)
    plot([k*Ts,(k-1)*Ts],[state(2,k),state(2,k-1)],'k-o');
    subplot(2,2,3)
    plot([k*Ts,(k-1)*Ts],[state(3,k),state(3,k-1)],'k-o');
    subplot(2,2,4)
    plot([k*Ts,(k-1)*Ts],[state(4,k),state(4,k-1)],'k-o');    
    pause(0.1);
end


function aux=model(x,F,m,b,k)
aux=[x(2);
    -(k/m)*x(1)-(b/m)*x(2)+F];
end

function varargout=MHE(Ts,x,u,y)
x_est=zeros(4,length(y));
x_est(:,1)=x;
N=length(y);
    for k=2:N
        %% Euler
        x_est(1:2,k)=x_est(1:2,k-1)+Ts*model(x_est(1:2,k-1),u(k),1,x_est(4,k-1),x_est(3,k-1));
        x_est(3:4,k)=x_est(3:4,k-1);
    end
    cost=norm(x_est(1:2,:)-y)^2;
    % Handle the variable number of outputs
    switch nargout
        case 1
            varargout{1} = cost;
        case 2
            varargout{1} = cost;
            varargout{2} = x_est(:,end);
        otherwise
            error('Invalid number of output arguments');
    end    
end
