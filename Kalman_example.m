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
%% Base
x = zeros(2, 100);
z = zeros(2, 100);
k = zeros(1, 100);
sigma = 0.1;
T = 1;

F = [1 T; 0 1]
Q = [sigma*sigma sigma*sigma/T;sigma*sigma/T 2*sigma*sigma/(T*T)]
%Q_k = sample_gaussian(0, Q, 99);
H = [1 0; 0 1]
h = [1 0]
%R = randn(2);
R = zeros(2);
%R_k = sample_gaussian(0, R, 99);

%% Intitialization
z(:,1) = [1, 8000];
z(1,2) = 2;
for i=1:99
    if i==1
        x(:,1)=z(:,1);
    elseif i==2
        x(1,2)=z(1,2);
        x(2,2)=(x(1,2)-x(1,1))/T;
    else
        z(:,i+1) = F*z(:,i);
        x(1,i+1) = h*z(:,i+1) + Q(1,1);
    end    
end

%% Prediction
% *coordinate prediction*
%
% $$x_{k+1} = Fx_k+Bu_k$$
%
% *error prediction:*
%
% $$P_k=FP_{k-1}F'+Q_k$$

%%
% P = zeros(2);
% for i=2:100
%     x(:,i)=F*x(:,i-1)+B*u(:,i);
%     P=F*P*F'+Q;
% end

%% Correction
% y_k = z_k + Hx_k
%%
y = zeros(2,100);
% *Optimal Kalman gain*
%
% $$K_k = P_kH'/(HP_kH'+R)$$
%
% *R* - measuring error covariance matrix
%
% *coordinate correction:*
%
% $$x_k = x_k + K_k(z_k-Hx_k)$$
%
% *error correction:*
%
% $$P_k=(I-K_kH)P_k$$
%
% *I* - Identity matrix
%
I = eye(2)
P = randn(2);

%% Test Example
for i=1:99
    [x(:,i+1), P] = Kalman_filter(z(:,i), F, P , H, Q, x(:,i), R);
end
y_k = y(1,:);
z_k = z(1,:);
x_k = x(1,:);
plot(k, y_k);
axis auto;
title('Фильтрация координаты');
xlabel('t, c');
ylabel('координата, м');
hold on;
plot(k, z_k, 'r')
hold on;
plot(k, x_k, 'y')
legend('отфильтрованная','измеренная', 'реальная');

x(1,:)
z(1,:)
%y(1,:)


%y_k = y(2,:);
z_k = z(2,:);
x_k = x(2,:);

% plot(k, y_k);
% axis auto;
% title('Фильтрация скорости');
% xlabel('t, c');
% ylabel('скорость, м/c');
% hold on;
% plot(k, z_k, 'r');
% hold on;
% plot(k, x_k, 'y');
% legend('отфильтрованная','измеренная', 'реальная');
