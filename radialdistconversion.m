x = [0:1:10];
y1 = [10:-1:0];
y2 = [0:1:10];
y12 = [0,0,0,0,0,0,0,0,0,0,0];

for i=1:(x+1)
    y12(i) = y1(i)*y2(i);
end

figure(1);
hold on
plot(x,y1);
plot(x,y2);
plot(x,y12);
hold off