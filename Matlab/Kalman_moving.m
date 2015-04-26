function dy = Kalman_moving(t,y)
    dy = zeros(6,1);
    dy(3) = 0;
    dy(2) = y(3);
    dy(1) = y(2);
    dy(6) = 0;
    dy(5) = y(6);
    dy(4) = y(5);
end

