close all
clear all
clc

x= -250:0.01:250;
k = [0.003,0.005,0.01,0.015,0.02];

figure
y1 = (1/k(1))*log(exp(k(1)*x)+1);
plot(x,y1)
grid on
hold on
y2 = (1/k(2))*log(exp(k(2)*x)+1);
plot(x,y2)
y3 = (1/k(3))*log(exp(k(3)*x)+1);
plot(x,y3)
y4 = (1/k(4))*log(exp(k(4)*x)+1);
plot(x,y4)
y5 = (1/k(5))*log(exp(k(5)*x)+1);
plot(x,y5)
plot(0:0.01:250,0:0.01:250,'k')
legend('0.001','0.005','0.01','0.015','0.02','f(x) = x')
