%% LKF Guel cortez 2024
% EKF for adaptive control of cart in one dimension


clearvars;
close all;

rng(0);

b=3; m=1; Ts=1e-2; L=10;
A=[[0,1];[0,-b/m]];
B=[0,1/m]';
C=[1,0];
D=0;

Ad=eye(2)+Ts*A;
Bd=Ts*B;

t=0:Ts:L;
w=1e-5;
v=1e-3;
Q=1e-5*eye(3); Q(3,3)=1e-3;
R=v^2;
P=1e-6*eye(3); P(3,3)=3e-2; 
alpha=0.03;

x_pred=zeros([3,length(t)]);
x_est=zeros([3,length(t)]);
y_pred=zeros([1,length(t)]);
Ce=[1,0,0];
x_pred(:,1)=[2,0.2,1.1]';

u_c = zeros([1,length(t)]);
k_c = [6 5];
xDeseada = [50 0];
x_c = zeros([2,length(t)]);
x_c(:,1) =[1,2]';
y_c = zeros([1,length(t)]);

p1=8;
p2=0;
for k=2:length(t)
    %% real system
    k_c(1)=(p1+p2)*m-x_est(3,k-1); k_c(2)=p1*p2*m;
    u_c(k-1) = k_c(1)*(xDeseada(1)-x_est(1,k-1)) + k_c(2)*(xDeseada(2)-x_est(2,k-1)) ;
    x_c(:,k)=Ad*x_c(:,k-1)+Bd*u_c(k-1)+sqrt(w)*randn([2,1]);
    y_c(k)=C*x_c(:,k)+sqrt(v)*randn;

    %% Kalman Filter
    x_pred(:,k)=fk(x_pred(:,k-1),u_c(k-1),m,Ts);
    y_pred(k)=Ce*x_pred(:,k);

    P_pred=Ak(x_pred(:,k),m,Ts)*P*Ak(x_pred(:,k),m,Ts)'+Q;
    P_y=Ce*P_pred*Ce'+R;
    P_xy=P_pred*Ce';

    L=P_xy/P_y;
    x_est(:,k)=x_pred(:,k)+L*(y_c(k)-y_pred(k));
    P=P_pred-L*P_y*L';
    
    %% Q learning
    Q=(1-alpha)*Q+alpha*L*(y_c(k)-y_pred(k))*(y_c(k)-y_pred(k))'*L';
end

subplot(4,1,1)
yyaxis left
plot(t,y_c,'k')
hold on
plot(t,x_est(1,:),'b.-')
ylabel('$x$','Interpreter','latex',FontSize=12)
yyaxis right 
plot(t,abs(y_c-x_est(1,:)),'r')
legend('y','x1 estimated','error')

subplot(4,1,2)
yyaxis left
plot(t,x_c(2,:),'k')
hold on
plot(t,x_est(2,:),'b.-')
ylabel('$\dot{x}$','Interpreter','latex',FontSize=12)
yyaxis right 
plot(t,abs(x_c(2,:)-x_est(2,:)),'r')
legend('x2 ','x2 estimate ','error')

subplot(4,1,3)
yyaxis left
plot(t,x_est(3,:),'b--')
hold on
yline(b);
ylabel('$b$','Interpreter','latex',FontSize=12)
yyaxis right 
plot(t,abs(b-x_est(3,:)),'r')
legend('b estimate ','b','error')

subplot(4,1,4)
plot(t,u_c,'k')
ylabel('$u$','Interpreter','latex',FontSize=12)
legend('Entrada')

function aux=Ak(x,m,T)
aux=[1 T 0;
    0 1-T*x(3)/m -T*x(2)/m;
    0 0 1];
end

function aux=fk(x,u,m,T)
aux=[x(1)+T*x(2);
    x(2)-T*x(3)*x(2)/m+T*u/m;
    x(3)];
end
