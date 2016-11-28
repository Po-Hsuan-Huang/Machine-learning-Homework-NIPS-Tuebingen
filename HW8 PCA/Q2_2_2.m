%% machine learning ex3.

%% Task 2.2 numerically estimate the cross-entropy of Gaussian Maximum likelihood model
%  generative model is Uniform model [-1 1]


%% initial parameter of Gaussian model in d-dimension,
clear;
d = 3; 
n = 5000; % sampling number


for i =1 :d
Rnd(i,:) =  -1+ 2*rand(1,n,1);  
end

Mean = mean(Rnd,2); % maximum likelihood mean. equal to real mean

%Cov =cov(Rnd'); % maximu likelihood covariance matrix is equal to real covariance matrix. 
sigma = zeros(d);
for i = 1:n
sigma = sigma + (Rnd(:,i)-Mean)*(Rnd(:,i)-Mean)';
end
Cov = (1/n) * (Rnd(:,i)-Mean)*(Rnd(:,i)-Mean)';

%%
p_model = @(x,y,z)sqrt(2*pi*Cov)^-d  * exp((-1/2)* ([x;y;z]-Mean)'*Cov*([x;y;z]-Mean)); 
 
figure()  
histogram(Rnd,100);
% hold on
% c = -1:0.01:1;
% %plot(c,sqrt(2*pi*Var)^-1  * exp((-1/2*Var)* (c-Mean).^2))
% hold off
p_true =  (1/2)^d;



%cross-entropy.
syms x y z
A = ([x;y;z]-Mean)'*Cov*([x;y;z]-Mean);
fun = @(x,y,z)  -log(sqrt(2*pi*det(Cov)).^-1  * exp((-1/2)* A))*p_true ;
CE =integral3(fun,-1,1,-1,1,-1,1);

fprintf('cross-entropy is %d',CE);




