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
sigma = 1;
sigmaR = 1;
T = 1;

F = [1 T T*T/2; 0 1 T; 0 0 1];
P = [sigmaR        3*sigmaR/(2*T)       sigmaR/(T*T);
    3*sigmaR/(2*T) 13*sigmaR/(2*T*T)    6*sigmaR/(T*T*T);
    sigmaR/(T*T)   6*sigmaR/(T*T*T)     6*sigmaR/(T*T*T*T)];
H = [1 0 0];
R = sigma;
Q = zeros(3);
I = [1 0 0; 0 1 0; 0 0 1];

x = zeros(3,N);
xF = zeros(3,N);
z = zeros(3,N);
Krv = 0;
Kvv = 0;
Kra = 0;
Kva = 0;
Kaa = 0;
pfk = zeros(3);

K = zeros(3, 1);
Gz = zeros(3, N);

k = zeros(1,N);
for i=2:N
    k(i) = k(i-1)+T;
end
%% Test
z(:,1) = [0 0 2];
z(:,2) = [1 2 2];
z(:,3) = [4 4 2];

x_k = zeros(3, 1);
x_k_1 = zeros(3, 1);

for i=1:N
    if i==1
        x(1,i) = z(1,i);
        xF(1,i) = z(1,i);
    elseif i==2
        x(1,i) = z(1,i);
        xF(1,i) = z(1,i);
    elseif i==3
        x(1,i) = z(1,i);
        xF(1,i) = z(1,i);
        %завязка
        x(3,i) = (x(1,i) - 2*x(1,i-1) + x(1,i-2))/(T*T);
        x(2,i) = (x(1,i) - x(1,i-2))/(2*T) + x(3,i)*T;
        %x(6,i) = (x(4,i) - 2*x(4,i-1) + x(4,i-2))/(T*T);
        %x(5,i) = (x(4,i) - x(4,i-2))/(2*T) + x(6,i)*T;
        
        xF(3,i) = (xF(1,i) - 2*xF(1,i-1) + xF(1,i-2))/(T*T);
        xF(2,i) = (xF(1,i) - xF(1,i-2))/(2*T) + xF(3,i)*T;
        
        Krv = 3*sigmaR/(2*T);
        Kra = sigmaR/(T*T);
        Kvv = 13*sigmaR/(2*T*T);
        Kva = 6*sigmaR/(T*T*T);
        Kaa = 6*sigmaR/(T*T*T*T);
    else
        %prediction
        P = F*P*F'+Q;
        sigmaR = sigmaR + 2*T*Krv + T*T*(Kra+Kvv) + T*T*T*Kva + T*T*T*T*Kaa/4;
        Krv = Krv + T*(Kra+Kvv) + 3*T*T*Kva/2 + T*T*T*Kaa/2;
        Kra = Kra + T*Kva + T*T*Kaa/2;
        Kvv = Kvv + 2*T*Kva + T*T*Kaa;
        Kva = Kva + T*Kaa;
        pfk = [sigmaR Krv Kra; Krv Kvv Kva; Kra Kva Kaa] + Q;
        
        z(:,i) = F*z(:,i-1);
        %[t,Y] = ode45(@Kalman_moving, [0 i], [0 0 2]);
        %z(:,i) = Y(end,:);
        %x_k = F*x_k_1;
        x(:,i) = F*x(:,i-1);
        xF(:,i) = F*xF(:,i-1); 
        
        i
        %correction
        if R == zeros(2)
            %G->H'
            P = 0
            pfk = 0
            
            %x->H'*z
            x(1,i) = z(1,i);
            xF(1,i) = z(1,i);
        elseif P == zeros(3)
            %G->0
            P
            pfk
        else
            %Kalman gain(blending factor)
            G = P*H'/(H*P*H'+R);
            W = 1/(sigmaR+R);
            Wv = Krv/(sigmaR+R);
            K = [W*sigmaR; W*Krv; W*Kra];
            
                       
            P = (I-K*H)*P
            sigmaR = W*sigmaR*R;
            Kaa = Kaa - W*Kra*Kra;
            Kvv = Kvv  - Wv*Krv;
            Kva = Kva - W*Krv*Kra;
            Krv = W*R*Krv;
            Kra = W*R*Kra;
            pfk = [sigmaR Krv Kra; Krv Kvv Kva; Kra Kva Kaa]
               
            n=i;  
            k11 = 3*(3*n*n - 3*n + 2)/(n*(n+1)*(n+2));
            k12 = - 18*(2*n-1)/(n*(n+1)*(n+2));
            k13 = 60/(n*(n+1)*(n+2));
            k22 = (12*(2*n-1)*(8*n-11))/(n*(n*n-4)*(n*n-1));
            k23 = 360/(n*(n*n-4)*(n+1));
            k33 = 720/(n*(n*n-4)*(n*n-1));
        
            Pnew = [k11 k12 k13; k12 k22 k23; k13 k23 k33]
%         %Pnew = vpa(Pnew, 10)
%         x_k_1 = x_k;
%x_k = x_k + K*(z(1,i) - H*x_k);
            x(:,i) = x(:,i) + G*(z(1,i) - H*x(:,i));
			xF(:,i) = xF(:,i) + K*(z(1,i)-H*xF(:,i));
        end
    end    
end

% subplot(1,3,1);
% plot(k, x(1,:), '-', k, z(1,:), 'r-.');
% title('Coordinate Filter');
% xlabel('t, c');
% ylabel('coordinate, m');
% legend('meassurement', 'filter');
% 
% subplot(1,3,2);
% plot(k, x(2,:), '-', k, z(2,:), 'r-.');
% title('Velocity Filter');
% xlabel('t, c');
% ylabel('velocity, m/c');
% legend('meassurement', 'filter');
% 
% subplot(1,3,3);
% plot(k, x(3,:), '-', k, z(3,:), 'r-.');
% title('Acceleration Filter');
% xlabel('t, c');
% ylabel('acceleration, m/c*c');
% legend('meassurement', 'filter');
