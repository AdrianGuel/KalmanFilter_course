%% LKF Guel cortez 2024
% LKF for cart model


clearvars;
close all;

b=1; m=2; Ts=1e-2; L=10;
%Cart_generative_model;

A=[[0,1];[0,-b/m]];
B=[0,1/m]';
C=[1,0];
D=0;

Ad=eye(2)+Ts*A;
Bd=Ts*B;

t=0:Ts:L;
w=1e-4;
v=1e-3;
Q=1e-5*eye(2);
R=v^2;
P=1e2*eye(2);
alpha=0.0003;

x_pred=zeros([2,length(t)]);
x_est=zeros([2,length(t)]);
y_pred=zeros([1,length(t)]);

x_pred(:,1)=[2,0]';

u_c = zeros([2,length(t)]);
k_c = [6 5];
xDeseada = [50 0];
x_c = zeros([2,length(t)]);
x_c(:,1) =[1,2]';
y_c = zeros([1,length(t)]);

for k=2:length(t)

    u_c(k-1) = k_c(1)*(xDeseada(1)-x_est(1,k-1)) + k_c(2)*(xDeseada(2)-x_est(2,k-1)) ;
    x_c(:,k)=Ad*x_c(:,k-1)+Bd*u_c(k-1)+sqrt(w)*randn([2,1]);
    y_c(k)=C*x_c(:,k)+sqrt(v)*randn;

    x_pred(:,k)=Ad*x_pred(:,k-1)+Bd*u_c(k-1);
    y_pred(k)=C*x_pred(:,k);

    P_pred=Ad*P*Ad'+Q;
    P_y=C*P_pred*C'+R;
    P_xy=P_pred*C';

    L=P_xy/P_y;
    x_est(:,k)=x_pred(:,k)+L*(y_c(k)-y_pred(k));
    P=P_pred-L*P_y*L';
    
    %% Q learning
    %Q=(1-alpha)*Q+alpha*L*(y_c(k)-y_pred(k))*(y_c(k)-y_pred(k))'*L';
end

subplot(3,1,1)
plot(t,y_c,'k')
hold on
plot(t,x_est(1,:),'r--')
legend('y','x1 estimated')

subplot(3,1,2)
plot(t,x_c(2,:),'k')
hold on
plot(t,x_est(2,:),'b--')
legend('x2 ','x2 estimate ')

subplot(3,1,3)
plot(t,u_c,'k')
legend('Entrada')
