%% Competition 2
%% something I learned, enforce positive definite condition to your covariance,
% or there will not be possible to production density.
%% don't forget it is the 'inverse' covariance matrix on the denominator 
%% only for 2 fold cross-validation
clear

load('C_train.mat')  % load given training data called X;

figure()

neu= reshape(X(1:99,1:200),3,200*33);
x1= neu(1,:);
x2= neu(2,:);
x3 =neu(3,:);
scatter3(x1,x2,x3);
%% Subsets and Subtrain
n  = 10 ;% n fold estimation

width = size(X,2)/n;

SubX= reshape(X,100,width,n) ;


for fold = 1:n   
Seg{fold} = SubX(:,:,fold) ;  
end


for fold  = 1:n
    tag= ones(1,n); tag(fold) =0; 
     
    Subtrain(:,:,fold)= [Seg{tag==1}];   %% I am sooo smart!
    
end



%% as required by exercise, mean is zero.

for fold = 1:n
    
    C_tr(:,:,fold)= cov(Subtrain(:,:,fold));  %% True covariance.

end


%% Multivariate Gaussian kernel model.

edge = -45:0.1:45;

sigma= 0.5;

d = 50; % how many kernels


for fold = 1:n
    
    
feature= zeros(d,length(edge3));
for i = 1:d
    feature(1,:)= edge3;
    xlin = linspace(-45,45,d);
    feature(i+1,:)=  (sqrt(2*pi*sigma^2))^-1.*exp(-0.5*(edge-xlin(i)).^2/sigma^2);
    
end
zTrain = feature(2:d+1,:)';


alpha = 1;
beta = 1; 

Covpost= inv(alpha*eye(d)+beta*(zTrain'*zTrain));  % covariance matrix
Meal= beta*Covpost*zTrain'*c3';             % posterior mean


Meanpost(fold,:)= Meal'./sum(Meal);             % normalization

%% weighted predictive mean using xspace
  
    Model(fold,:)= zeros(1,length(edge3));

 for i=1:d   
    Model(fold,:)=Model(fold,:) + Meanpost(i).*(sqrt(2*pi*sigma^2))^-1.*exp(-0.5*(edge3-xlin(i)).^2/sigma^2);    
 end   
 


end







