%% simple plots of gaussian random variables
% Ref Gregory plett
figure
ybar = [1; 2]; covar = [2, 0.75; 0.75, 1];
A = chol(covar,'lower');
for k = 1:5000,
x = randn([2, 1]);
y = ybar + A*x;
plot(y(1),y(2),'k.'); hold on
end

figure
ybar = [1; 2]; covar = [2, 0.75; 0.75, 1];
[L,D] = ldl(covar);
for k = 1:5000,
x = randn([2, 1]);
y = ybar + (L*sqrt(D))*x;
plot(y(1),y(2),'b.'); hold on
end
