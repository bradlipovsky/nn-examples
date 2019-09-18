clear;

% Neural network estimation of glacier thickness in three lines of code:
load glaciers
net = feedforwardnet(10);
[net,tr] = train(net,X,T);

% mean(abs(net(X)-T))