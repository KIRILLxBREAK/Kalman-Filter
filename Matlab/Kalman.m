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
sigma = 0;
sigmaR = 1;
T = 1;

F = [1 T; 0 1];
P = [sigmaR sigmaR/T; sigmaR/T 2*sigmaR/(T*T)];
H = [1 0];
R = sigma;
Q = zeros(2);
I = [1 0; 0 1];

x = zeros(2,N);
xF = zeros(2,N);
z = zeros(2,N);
Krv = 0;
Kvv = 0;
pfk = zeros(2);

K = zeros(2, 1);


k = zeros(1,N);
for i=2:N
    k(i) = k(i-1)+T;
end
%% Test
z(:,1) = [1 8000];
z(:,2) = [2 1];
for i=1:N
    if i==1
        x(:,1) = z(:,1);
        xF(:,1) = z(:,1);
    elseif i==2
        x(1,2)=z(1,2);
        x(2,2)=(x(1,2)-x(1,1))/T;
        xF(1,2)=z(1,2);
        xF(2,2)=(xF(1,2)-xF(1,1))/T;
        i
        Krv = sigmaR/T;
        Kvv = 2*sigmaR/T*T;
        pfk = [sigmaR Krv; Krv Kvv]
    else
        %prediction
        P = F*P*F' + Q; 
        sigmaR = (sigmaR + 2*Krv*T + Kvv*T*T);
        Krv = Krv + Kvv*T;  
        pfk = [sigmaR Krv; Krv Kvv] + Q;
        
        z(:,i) = F*z(:,i-1);
        x(:,i) = F*x(:,i-1);
        xF(:,i) = F*xF(:,i-1);

        %correction
        i
        if R == zeros(2)
            %G->H'
            P = 0
            pfk = 0
            
            %x->H'*z
            x(1,i) = z(1,i);
            xF(1,i) = z(1,i);
        elseif P == zeros(2)
            %G->0
            P
            pfk
        else
            %Kalman gain(blending factor)
            G = P*H'/(H*P*H'+R)
            W = 1/(sigmaR+R);
            Wv = Krv/(sigmaR+R);
        
            P = (I-G*H)*P
            sigmaR = W*sigmaR*R;
            Kvv = Kvv - Wv*Krv;
            Krv = W*R*Krv;        
            pfk = [sigmaR Krv; Krv Kvv]
%         K = [W*sigmaR; Wv]
%         pfk = (I-K*H)*pfk         
%         sigmaR = pfk(1,1);
%         Krv = pfk(1,2);
%         Kvv = pfk(2,2);

        x(:,i) = x(:,i) + G*(z(1,i) - H*x(:,i));
        xF(:,i) = xF(:,i) + K*(z(1,i)-H*xF(:,i));
        end
    end
end

plot(k, z(1,:), 'r')
title('Kalman filter');
xlabel('t, c');
ylabel('coordinate, m');
hold on;
plot(k, x(1,:), 'y')
legend('meassurement', 'filter');