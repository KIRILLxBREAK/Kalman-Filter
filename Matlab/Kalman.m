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
x = zeros(3,100);
z = zeros(3,100);
%x = zeros(1, 100);
%z = zeros(1, 100);
sigma = 1;
sigma0 = 1;
T = 1;

k = zeros(1,100);
for i=2:100
    k(i) = k(i-1)+T;
end

F = [1 T T*T/2; 0 1 T; 0 0 1];
P = [sigma*sigma 3*sigma*sigma/(2*T) sigma*sigma/(T*T); 3*sigma*sigma/(2*T) 13*sigma*sigma/(2*T*T) 6*sigma*sigma/(T*T*T); sigma*sigma/(T*T) 6*sigma*sigma/(T*T*T) 6*sigma*sigma/(T*T*T*T)];
H = [1 0 0];
R = sigma;
Q = zeros(3);
q = 0;
I = [1 0 0; 0 1 0; 0 0 1];
K = zeros(3,1);

%% Test
z(:,1) = [0 0 2];
z(:,2) = [1 2 2];
z(:,3) = [4 4 2];

x_k = zeros(3, 1);
x_k_1 = zeros(3, 1);

for i=1:100
    if i==1
        x_k = z(:,1);
        x(:,1) = x_k;
    elseif i==2
        x_k = z(:,2);
        x(:,2) = x_k;
    elseif i==3
        x_k_1 = z(:,3);
        x(:,3) = x_k_1;
    else
        %prediction
        %z(:,i) = F*z(:,i-1);
        [t,Y] = ode45(@Kalman_moving, [0 i], [0 0 2]);
        z(:,i) = Y(end,:);
        x_k = F*x_k_1;
        x(:,i) = F*x_k;
        
        P = F*P*F'+Q;
        
        %correction
        K = P*H'/(H*P*H'+R);
        
        x_k = x_k + K*(z(1,i) - H*x_k);
        P = (I-K*H)*P;       
        
        i
        n=i;
        P
        %P = vpa(P, 10)
        
        k11 = 3*(3*n*n - 3*n + 2)/(n*(n+1)*(n+2));
        k12 = - 18*(2*n-1)/(n*(n+1)*(n+2));
        k13 = 60/(n*(n+1)*(n+2));
        k22 = (12*(2*n-1)*(8*n-11))/(n*(n*n-4)*(n*n-1));
        k23 = 360/(n*(n*n-4)*(n+1));
        k33 = 720/(n*(n*n-4)*(n*n-1));
        
        Pnew = [k11 k12 k13; k12 k22 k23; k13 k23 k33]
        %Pnew = vpa(Pnew, 10)
        x_k_1 = x_k;
    end    
end

subplot(1,3,1);
plot(k, x(1,:), '-', k, z(1,:), 'r-.');
title('Coordinate Filter');
xlabel('t, c');
ylabel('coordinate, m');
legend('meassurement', 'filter');

subplot(1,3,2);
plot(k, x(2,:), '-', k, z(2,:), 'r-.');
title('Velocity Filter');
xlabel('t, c');
ylabel('velocity, m/c');
legend('meassurement', 'filter');

subplot(1,3,3);
plot(k, x(3,:), '-', k, z(3,:), 'r-.');
title('Acceleration Filter');
xlabel('t, c');
ylabel('acceleration, m/c*c');
legend('meassurement', 'filter');