clear all;
close all;

x = [-3:.1:10];
norm = normpdf(x,5,1);

figure(1);
plot(x,norm);