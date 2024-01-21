%Information length from stochastic simulation of Kramers equation
%      dx = y dt + sqrt(D_11) dW,   
%      dy = (-omega²*x-gamma*y) dt + sqrt(D_22) dW,
%      x(0) = x0, y(0) = y0,
% using an initial gaussian distribution p(x,y,0)~N(0,Sigma)

close all
clear all
clc
%%
figure
set(gcf,'color','w');
ax1 = subplot(2,2,[2 4]);
%set(ax1,'YScale','log')
%set(ax1,'XScale','log')
hold(ax1,'on')
grid(ax1,'on')
xlabel(ax1,'$x\sim \mathcal{N}(\mu_x, \mu_y,\Sigma)$','Interpreter','Latex','FontSize', 14)
ylabel(ax1,'$y\sim \mathcal{N}(\mu_x, \mu_y,\Sigma)$','Interpreter','Latex','FontSize', 14)
xlim(ax1,[-1 3])

ax2 = subplot(2,2,1);
%set(ax2,'YScale','log')
%set(ax1,'XScale','log')
hold(ax2,'on')
grid(ax2,'on')
xlabel(ax2,'$t$','Interpreter','Latex','FontSize', 14)
%ylabel(ax2,'$\mathcal{E}(t)$','Interpreter','Latex','FontSize', 14)
ylabel(ax2,'$\mathcal{E}(t)$','Interpreter','Latex','FontSize', 14)

ax3 = subplot(2,2,3);
%set(ax3,'YScale','log')
%set(ax1,'XScale','log')
hold(ax3,'on')
grid(ax3,'on')
xlabel(ax3,'$t$','Interpreter','Latex','FontSize', 14)
%ylabel(ax3,'$\mathcal{L}(t)$','Interpreter','Latex','FontSize', 14)
ylabel(ax3,'$\mathcal{L}(t)$','Interpreter','Latex','FontSize', 14)
%%
timerVal = tic;
%Stochastic simulation
omega=1;
gamma=2;
n=2;
A=[0,1;-omega^2,-gamma];
sigma0=[0.01 0; 0 0.01];
D=[0.1 0; 0 0.1];
x0=1;
y0=2;
X0=[x0,y0];
N=1e3;
tmax=5;
[t,x,y]=Kramers_euler_maruyama (omega,gamma, D, x0, y0, sigma0, tmax, N);    
opts = odeset('RelTol',1e-8,'AbsTol',1e-12);
[t,yx] = ode45(@(t,y) Sigmax(t,y,A,D,n), t,  sigma0, opts);
[t,M] = ode45(@(t,y) Mux(t,y,A), t, X0, opts);

tEnd = toc(timerVal); 

plot(ax1,x,y,'k.',M(:,1),M(:,2),'b')
%save('KramersSS','t','x','y','yx','M')
Ts=diff(t);
E=zeros(1,length(t));
for k=2:length(t)
    GMModel0 = fitgmdist([x(:,k-1),y(:,k-1)],1);
    GMModel1 = fitgmdist([x(:,k),y(:,k)],1);
    gmPDFt = @(x,y) arrayfun(@(x0,y0) (sqrt(pdf(GMModel1,[x0 y0]))-sqrt(pdf(GMModel0,[x0 y0]))).^2,x,y);    
    %fcontour(ax1,gmPDFt)
    E(k)=4*integral2(gmPDFt,-inf,inf,-inf,inf)/(Ts(1)^2);
end

IL=cumtrapz(t,E);
plot(ax2,t,E,'k')
plot(ax3,t,IL,'k')

%% 
function dydt = Sigmax(t,y,A,D,n)
   At=A';
   aux2=reshape(y,2,2);
   dydt=kron(eye(n,n),A)*aux2(:)+kron(eye(n,n),aux2)*At(:)+2*D(:);
end
%  
function dydt = Mux(t,y,A)
   dydt=A*y;
end
