%% Kalman Filter
% Filtering the measured parameter
% using MATLAB(R).
%% Basic data
% $$z_k = x_k + h_k$$
%
% *z[2]* - observation data
%
% *w[2]* - measuring error
%
% *x[2]* - our hidden variable, which we want to filter;
% *x(1)* - coordinate;
% *x(2)* - velocity
%
% $$x_{k+1} = Fx_k+w_k$$
%% Intitialization
x = zeros(2,100);
z = zeros(2,100);
%x = zeros(1, 100);
%z = zeros(1, 100);
sigma = 1;
sigma0 = 1;
T = 1;

k = zeros(1,100);
for i=2:100
    k(i) = k(i-1)+T;
end

F = [1 T; 0 1]
P = [sigma*sigma sigma*sigma/T;sigma*sigma/T 2*sigma*sigma/(T*T)]
%Q_k = sample_gaussian(0, Q, 99);
%H = [1 0; 0 1]
H = [1 0]
%R = randn(2);
R = sigma;
Q = zeros(2);
q = 0;
I = [1 0; 0 1];
%K = zeros(1,100);
%R_k = sample_gaussian(0, R, 99);

%% Test
z(:,1) = [1 8000];
z(:,2) = [2 1];
%z(1) = 1;
%z(2) = 2;
for i=1:100
    if i==1
        x(:,1) = z(:,1)
        %x(1)=z(1);
    elseif i==2
        x(1,2)=z(1,2);
        %x(2)=z(2);
        x(2,2)=(x(1,2)-x(1,1))/T;
        i
        %Kaa = 1;
        %K = 2;
        %P = [sigma0 Kaa;Kaa K]
        P
    else
        %prediction
        z(:,i) = F*z(:,i-1);
        x(:,i) = F*x(:,i-1);
        %z(i) = z(i-1)+ (z(i-1)-z(i-2))*T;
        %x(i) = x(i-1)+ (x(i-1)-x(i-2))*T;
        
        P = F*P*F'+Q;
        %sigma0 = (sigma0 + 2*Kaa*T + K*T*T);
        %Kaa = Kaa + K*T;        
        
        %correction
        K = P*H'/(H*P*H'+R);
        %W = sigma0/(sigma0+sigma);
        %Wa = Kaa/(sigma0+sigma);
        
        x(:,i) = x(:,i) + K*(z(1,i) - H*x(:,i));
        %x(i) = x(i) + W*(z(i)-x(i));
        
        i
        P = (I-K*H)*P
        %sigma0 = W*sigma;
        %Kaa = Wa*sigma;
        %K = W*K;
        %P = [sigma0 Kaa; Kaa K]
    end    
end

z_k = z(1,:);
x_k = x(1,:);
plot(k, z_k, 'r')
axis auto;
title('Kalman filter');
xlabel('t, c');
ylabel('coordinate, m');
hold on;
plot(k, x_k, 'y')
legend('meassurement', 'filter');

%x(:,:)
%z(:,:)