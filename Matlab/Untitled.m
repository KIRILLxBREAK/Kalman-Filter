x = zeros(6, 100);
z = zeros(6, 100);

sigma = 1;
T = 1;

R = [sigma 0; 0 sigma];
P = zeros(6);
H = [1 0 0 0 0 0; 0 1 0 0 0 0];
F = [ 1 0 T 0 T*T/2 0; 0 1 0 T 0 T*T/2; 0 0 1 0 T 0; 0 0 0 1 0 T; 0 0 0 0 1 0; 0 0 0 0 0 1];
I = eye(6);

z(:,1) = [1 1 8000 8000 0 0];
z(:,2) = [2 2 1 1 1 1];

for i=1:100
    if i==1
        x(:,1) = z(:,1);
    elseif i==2
        x(1,2)=z(1,2);
        x(2,2)=(x(1,2)-x(1,1))/T;
    else
        %prediction
        z(:,i) = F*z(:,i-1);
        x(:,i) = F*x(:,i-1);
        
        P = F*P*F';       
        
        %correction
        K = P*H'/(H*P*H'+R);       
        
        x(:,i) = x(:,i) + K*(z(1,i) - H*x(:,i));
        
        P = (I-K*H)*P
    end    
end