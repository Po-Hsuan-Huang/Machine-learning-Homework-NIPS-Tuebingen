%% machine learning ex3.

%% Task 2.2 numerically estimate the cross-entropy of Gaussian Maximum likelihood model
%  generative model is Uniform model [-1 1]


%% initial parameter of Gaussian model in 1-dimension,

n = 500000; % sampling number

Rnd =  -1+ 2*rand(1,n,1);  


Mean = mean(Rnd); % maximum likelihood mean. equal to real mean

Var = var(Rnd); % maximu likelihood covariance matrix is equal to real covariance matrix. 


%%
p_model = @(x)sqrt(2*pi*Var)^-1  * exp((-1/2*Var)* (x-Mean)^2); 
 
figure()  
histogram(Rnd,100);
% hold on
% c = -1:0.01:1;
% %plot(c,sqrt(2*pi*Var)^-1  * exp((-1/2*Var)* (c-Mean).^2))
% hold off
p_true =  1/2;



%cross-entropy.
fun = @(x)-log(sqrt(2*pi*Var).^-1  * exp((-1/(2*Var)).* (x-Mean).^2))*p_true ;
CE = integral(fun ,-1,1);

fprintf('cross-entropy is %d',CE);
