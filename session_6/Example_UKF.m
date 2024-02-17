%% UKF example from 
% Gregory plett Plett, G. L. (2015). Battery management systems, Volume II: Equivalent-circuit methods. Artech House.

% Define size of variables in model with w = 1 and v = 2. Nx = 1; % state = 1x1 scalar
Nxa = 3; % augmented state has also w(k) and v(k) contributions
Nz = 1; % output = 1x1 scalar
% Some constants for the SPKF algorithm. Use standard values for
% cases with Gaussian noises. (These are the weighting matrices
% comprising the values of alpha(c) and alpha(m) organized in a way to
% make later computation efficient).
h = sqrt(3);
Wmx(1) = (h*h-Nxa)/(h*h); Wmx(2) = 1/(2*h*h); Wcx=Wmx;
Wmxz = [Wmx(1) repmat(Wmx(2),[1 2*Nxa])]';
% Initialize simulation variables
SigmaW = 1; % Process noise covariance
SigmaV = 2; % Sensor noise covariance
maxIter = 40;
xtrue =2+ randn(1); % Initialize true system initial state
xhat = 2; % Initialize Kalman filter initial estimate
SigmaX = 1; % Initialize Kalman filter covariance
% Reserve storage for variables we might want to plot/evaluate
xstore = zeros(maxIter+1,length(xtrue)); xstore(1,:) = xtrue;
xhatstore = zeros(maxIter,length(xhat));
SigmaXstore = zeros(maxIter,length(xhat)^2);

for k = 1:maxIter
    % SPKF Step 1a: State estimate time update
    % 1a-i: Calculate augmented state estimate, including ...
    xhata = [xhat; 0; 0]; % process and sensor noise mean
    % 1a-ii: Get desired Cholesky factor
    Pxa = blkdiag(SigmaX,SigmaW,SigmaV);
    sPxa = chol(Pxa,'lower');
    % 1a-iii: Calculate sigma points (strange indexing of xhat to avoid
    % "repmat" call, which is very inefficient in Matlab)
    X = xhata(:,ones([1 2*Nxa+1])) + h*[zeros([Nxa 1]), sPxa, -sPxa];
    % 1a-iv: Calculate state equation for every element
    % Hard-code equation here for efficiency
    Xx = sqrt(5+X(1,:)) + X(2,:);
    xhat = Xx*Wmxz;
    % SPKF Step 1b: Covariance of prediction
    Xs = (Xx(:,2:end) - xhat(:,ones([1 2*Nxa])))*sqrt(Wcx(2));
    Xs1 = Xx(:,1) - xhat;
    SigmaX = Xs*Xs' + Wcx(1)*Xs1*Xs1';
    % [Implied operation of system in background, with
    % input signal u, and output signal z]
    w = chol(SigmaW)'*randn(1);
    v = chol(SigmaV)'*randn(1);
    ztrue = xtrue^3 + v; % z is based on present x and u
    xtrue = sqrt(5+xtrue) + w; % future x is based on present u
    % SPKF Step 1c: Create output estimate
    % Hard-code equation here for efficiency
    Z = Xx.^3 + X(3,:);
    zhat = Z*Wmxz;
    % SPKF Step 2a: Estimator gain matrix
    Zs = (Z(:,2:end) - zhat*ones([1 2*Nxa])) * sqrt(Wcx(2));
    Zs1 = Z(:,1) - zhat;
    SigmaXZ = Xs*Zs' + Wcx(1)*Xs1*Zs1';
    SigmaZ = Zs*Zs' + Wcx(1)*Zs1*Zs1';
    Lx= SigmaXZ/SigmaZ;
    % SPKF Step 2b: Measurement state update
    xhat = xhat + Lx*(ztrue-zhat); % update prediction to estimate
    % SPKF Step 2c: Measurement covariance update
    SigmaX = SigmaX - Lx*SigmaZ*Lx';
    % [Store information for evaluation/plotting purposes]
    xstore(k+1,:) = xtrue;
    xhatstore(k,:) = xhat;
    SigmaXstore(k,:) = SigmaX(:);
end

figure(1); clf;
plot(0:maxIter-1,xstore(1:maxIter),'k-',0:maxIter-1,xhatstore,'b--', ...
0:maxIter-1,xhatstore+3*sqrt(SigmaXstore),'m-.',...
0:maxIter-1,xhatstore-3*sqrt(SigmaXstore),'m-.'); grid;
legend('true','estimate','bounds'); xlabel('Iteration'); ylabel('State');
title('Sigma-point Kalman filter in action');
figure(2); clf;
plot(0:maxIter-1,xstore(1:maxIter)-xhatstore,'b-',0:maxIter-1, ...
3*sqrt(SigmaXstore),'m--',0:maxIter-1,-3*sqrt(SigmaXstore),'m--');
grid; legend('Error','bounds');
title('SPKF Error with bounds');
xlabel('Iteration'); ylabel('Estimation Error');
