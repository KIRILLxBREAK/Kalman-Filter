function dy = Kalman_moving(t,y)
    dy = zeros(3,1);
    dy(3) = 0;
    dy(2) = y(3);
    dy(1) = y(2);
end