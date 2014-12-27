function[x, P] = Kalman_filter(z, F, P , H, Q , x, R)
    % prediction
    x = F*x;
    P = F*P*F' + Q;
    % correction
    K= P*H'/(H*P*H' + R);
    x = x + K*(z- H*x);
    I = eye(2);
    P = (I - K*H)*P;
end

%%
% for k=1:99
%     % prediction
%     x(k+1) = x(k)*F;
%     P(k) = F*P(k)*F' + Q
%     % correction
%     K(k)= P(k+1)*H'/(H*P(k+1)*H' + R)
%     x(k+1)=x(k+1) + K(k)*(z(k)- H*x(k+1))
%     P(k)= (I - K(k)*H)*P(k)
% end