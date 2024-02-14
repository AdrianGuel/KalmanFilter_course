%Extended kalman filter example
% Aplication in a simple mass-spring-damper system
%Author: Adrian-Josue Guel-Cortez 2024
%In this case, we know m but we do not know k and b.

%Solution to the model
Ts=1e-2;
t=0:Ts:10;
rng(0);
F=zeros(1,length(t));
%F=0.1*sin(t);
w=1e-4;
k=1.9; b=0.5; m=1;
y=zeros(2,length(t));
y(:,1)=[1;0.1];
for j=2:length(t)
    y(:,j)=y(:,j-1)+Ts*model(y(:,j-1),F(j-1),m,b,k,w);
end

%Kalman filter
n=4;
P_previous=1e3*eye(n,n);
Q=1e3*eye(n,n);
R=0.01; H=[1,0,0,0]; V=1;
W=eye(n,n);

l=1;
x_previous=[0.1,0.5,0.2,0.1];
x=zeros(l*length(t),n);
k=2;
for i=2:l*length(t)
    x_model=f(x_previous,F(i-1),[0 0 0 0 0],m,Ts);
    P_kminus=A(x_model,m,Ts)*P_previous*A(x_model,m,Ts)'+W*Q*W';
    K_k=P_kminus*H'*((H*P_kminus*H'+V*R*V')^(-1));
    x_k=x_model+K_k*(y(1,i)-x_model(1));
    P_k=(eye(n,n)-K_k*H)*P_kminus;
    x(i,:)=x_k';
    x_previous=x_k;
    P_previous=P_k;
    Q=(1-0.03)*Q+0.03*K_k*(y(1,i)-x_model(1))*(y(1,i)-x_model(1))'*K_k';
end

figure
set(gcf,'color','w');
subplot(3,2,[1 2])
plot(t,y(1,:),'b')
hold on
plot(t,x(:,1),'--r','LineWidth',2)
ylabel('$x$','Interpreter','latex',FontSize=12)
legend('truth','estimated')

subplot(3,2,[3 4])
plot(t,y(2,:),'b')
hold on
plot(t,x(:,2),'--r','LineWidth',2)
ylabel('$\dot{x}$','Interpreter','latex',FontSize=12)
legend('truth','estimated')

subplot(3,2,5)
k=1.9; b=0.5;
plot(t,x(:,3),'--r','LineWidth',2)
yline(k,'b')
ylabel('$k$','Interpreter','latex',FontSize=12)
legend('estimated','truth')

subplot(3,2,6)
plot(t,x(:,4),'--r','LineWidth',2)
yline(b,'b')
ylabel('$b$','Interpreter','latex',FontSize=12)
legend('estimated','truth')
k_estimated=x(end,3);
b_estimated=x(end,4);

function aux=A(x,m,T)
aux=[1 T 0 0;
    -T*x(3)/m 1-T*x(4)/m -T*x(1)/m -T*x(2)/m;
    0 0 1 0;
    0 0 0 1];
end

function aux=f(x,F,w,m,T)
aux=[x(1)+T*x(2)+T*w(1);
    x(2)-T*x(3)*x(1)/m-T*x(4)*x(2)/m+T*F/m;
    x(3);
    x(4)];
end

function aux=model(x,F,m,b,k,w)
aux=[x(2)+w*randn;
    -(k/m)*x(1)-(b/m)*x(2)+F+w*randn];
end
