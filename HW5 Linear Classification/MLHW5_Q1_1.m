%%
dist=@(s) 1/(1+exp(-s));
figure()
subplot(2,1,1)
fplot(dist,[0,10])
title('logistic funtion f(s)')
xlabel('s')
ylabel('f(s)')
subplot(2,1,2)
logdis=@(s) log(1/(1+exp(s)));
fplot(logdis,[0,10])
title('llog logistic funtion f(s)')
xlabel('s')
ylabel('f(s)')