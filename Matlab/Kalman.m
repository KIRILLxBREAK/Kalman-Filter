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
N = 100;
x = zeros(6,N);
xF = zeros(6,N);
z = zeros(6,N);

x_k = zeros(6,1);
x_kF = zeros(6,1);
x_k_1 = zeros(6,1);
x_k_1F = zeros(6,1);

sigma = 0.1;
sigmaF = 0.1;
T = 1;
%% moving
tspan = 1:T:N*T; 
[t, Y] = ode45(@Kalman_moving, tspan, [0 0 2 0 0 2]);
plot(t, Y(:,1));
[fi, Rad] = cart2pol(Y(:,1), Y(:,4));

F = [1 T T*T/2 0 0    0;
     0 1   T   0 0    0;
     0 0   1   0 0    0;
     0 0   0   1 T  T*T/2;
     0 0   0   0 1    T;
     0 0   0   0 0    1;];
sigmaX  = (sigma*cos(fi(3)))^2     + (Rad(3)*sigmaF*sin(fi(3)))^2;
sigmaY  = (sigma*sin(fi(3)))^2     + (Rad(3)*sigmaF*cos(fi(3)))^2;
sigmaXY = (sin(2*fi(3))*sigma^2)/2 - (sin(2*fi(3))*(Rad(3)*sigmaF)^2)/2;

PF = [sigmaX         3*sigmaX/(2*T)     sigmaX/(T*T)       sigmaXY         3*sigmaXY/(2*T)    sigmaXY/(T*T);
     3*sigmaX/(2*T)  13*sigmaX/(2*T*T)  6*sigmaX/(T*T*T)   3*sigmaXY/(2*T) 13*sigmaXY/(2*T*T) 6*sigmaXY/(T*T*T);
     sigmaX/(T*T)    6*sigmaX/(T*T*T)   6*sigmaX/(T*T*T*T) sigmaX/(T*T)    6*sigmaX/(T*T*T)   6*sigmaX/(T*T*T*T);
     sigmaXY         3*sigmaXY/(2*T)    sigmaXY/(T*T)      sigmaY          3*sigmaY/(2*T )    sigmaY/(T*T);
     3*sigmaXY/(2*T) 13*sigmaXY/(2*T*T) 6*sigmaXY/(T*T*T)  3*sigmaY/(2*T)  13*sigmaY/(2*T*T)  6*sigmaY/(T*T*T);
     sigmaX/(T*T)    6*sigmaX/(T*T*T)   6*sigmaX/(T*T*T*T) sigmaY/(T*T)    6*sigmaY/(T*T*T)   6*sigmaY/(T*T*T*T)];

P = [sigma^2        3*sigma^2/(2*T)    sigma^2/(T*T)       0               0                  0;
    3*sigma^2/(2*T) 13*sigma^2/(2*T*T) 6*sigma^2/(T*T*T)   0               0                  0;
    sigma^2/(T*T)   6*sigma^2/(T*T*T)  6*sigma^2/(T*T*T*T) 0               0                  0;
    0               0                  0                   sigma^2         3*sigma^2/(2*T )   sigma^2/(T*T);
    0               0                  0                   3*sigma^2/(2*T) 13*sigma^2/(2*T*T) 6*sigma^2/(T*T*T);
    0               0                  0                   sigma^2/(T*T)   6*sigma^2/(T*T*T)  6*sigma^2/(T*T*T*T)];

H = [1 0 0 0 0 0;
     0 0 0 1 0 0];
R = [sigma^2    0;
     0          sigma^2];
RF = [sigmaX    sigmaXY;
      sigmaXY   sigmaY];
Q = zeros(6);
q = [0 0];
I = [1 0 0 0 0 0;
     0 1 0 0 0 0;
     0 0 1 0 0 0;
     0 0 0 1 0 0;
     0 0 0 0 1 0;
     0 0 0 0 0 1];
K = zeros(6,2);
KF = zeros(6,2);

%% filtering
format long g;

for i=1:N
    if i==1
        %measurement
        [z(1,i), z(4,i)] = pol2cart(fi(i),Rad(i));
        
        x_k = z(:,i);
        x(:,i) = z(:,i);
        x_kF = z(:,i);
        xF(:,i) = z(:,i);
    elseif i==2
        %measurement
        [z(1,i), z(4,i)] = pol2cart(fi(i),Rad(i));
        
        x_k_1 = z(:,i);
        x(:,i) = z(:,i);
        x_k_1F = z(:,i);
        xF(:,i) = z(:,i);
    elseif i==3
        [z(1,i), z(4,i)] = pol2cart(fi(i),Rad(i));
        z(3,i) = (z(1,i) - 2*z(1,i-1) + z(1,i-2))/(T*T);
        z(2,i) = (z(1,i) - z(1,i-2))/(2*T) + z(3,i)*T;
        z(6,i) = (z(4,i) - 2*z(4,i-1) + z(4,i-2))/(T*T);
        z(5,i) = (z(4,i) - z(4,i-2))/(2*T) + z(6,i)*T;
        
        x_k_1 = z(:,i);
        x(:,i) = z(:,i);
        x_k_1F = z(:,i);
        xF(:,i) = z(:,i);
    else    
        %measurement      
        [z(1,i), z(4,i)] = pol2cart(fi(i),Rad(i));
        z(3,i) = (z(1,i) - 2*z(1,i-1) + z(1,i-2))/(T*T);
        z(2,i) = (z(1,i) - z(1,i-2))/(2*T) + z(3,i)*T;
        z(6,i) = (z(4,i) - 2*z(4,i-1) + z(4,i-2))/(T*T);
        z(5,i) = (z(4,i) - z(4,i-2))/(2*T) + z(6,i)*T;
        
        i
        %prediction 
        x_k = F*x_k_1;
        x(:,i) = x_k;
        x_kF = F*x_k_1F;
        xF(:,i) = x_kF;
        
        P = F*P*F'+Q;
        PF = F*PF*F'+Q;
        
        %correction
        K = P*H'/(H*P*H'+R);
        KF = PF*H'/(H*PF*H'+RF);
        
        x_kF = x_kF + KF*([z(1,i); z(4,i)] - H*x_kF)
        x_k = x_k + K*([z(1,i); z(4,i)] - H*x_k)
        
        P = (I-K*H)*P;
        PF = (I-KF*H)*PF;
        
%         disp('KF=');
%         disp(num2str(KF,'% f'));
%         disp('K=');
%         disp(num2str(K,'% f'));
%         PF = vpa(PF, 10)
%         P = vpa(P, 10)
        disp('PF=');
        disp(num2str(PF,'% f'));
        disp('P=');
        disp(num2str(P,'% f'));
    end    
end

% subplot(2,3,1);
% plot(k, x(1,:), '-', k, z(1,:), 'r-.');
% title('x Coordinate Filter');
% xlabel('t, c');
% ylabel('coordinate, m');
% legend('meassurement', 'filter');
% 
% subplot(2,3,2);
% plot(k, x(2,:), '-', k, xF(2,:), 'r-.');
% title('x Velocity Filter');
% xlabel('t, c');
% ylabel('velocity, m/c');
% legend('meassurement', 'filter');
% 
% subplot(2,3,3);
% plot(k, x(3,:), '-', k, xF(3,:), 'r-.');
% title('x Acceleration Filter');
% xlabel('t, c');
% ylabel('acceleration, m/c*c');
% legend('meassurement', 'filter');
% 
% subplot(2,3,4);
% plot(k, x(4,:), '-', k, xF(4,:), 'r-.');
% title('y Coordinate Filter');
% xlabel('t, c');
% ylabel('coordinate, m');
% legend('meassurement', 'filter');
% 
% subplot(2,3,5);
% plot(k, x(5,:), '-', k, xF(5,:), 'r-.');
% title('y Velocity Filter');
% xlabel('t, c');
% ylabel('velocity, m/c');
% legend('meassurement', 'filter');
% 
% subplot(2,3,6);
% plot(k, x(6,:), '-', k, xF(6,:), 'r-.');
% title('y Acceleration Filter');
% xlabel('t, c');
% ylabel('acceleration, m/c*c');
% legend('meassurement', 'filter');