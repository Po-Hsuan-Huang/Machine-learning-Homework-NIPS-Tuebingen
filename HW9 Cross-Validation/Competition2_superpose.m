%% Competition 2   superpose
%% variable description:
% mean cross-validation test error is Err_test.
%  err_test is the cross-validation test error for each fold.
% mean log-likelyhood for Maximum Likelyhood method is LL_ml_tt
% mean log-likelyhood for Superpoision method is LL_super_tt
% mean log-likelyhood for Rotation method(MAP) is LL_rotate_tt
% difference of log-likelyhood between LL_ml and LL_rotate is alles


%% something I learned, enforce positive definite condition to your covariance,
% or there will not be possible to production density.
%% only for 2 fold cross-validation
clear

load('C_train.mat')  % load given training data called X;






n  = 2 ;% n fold estimation

width = size(X,2)/n;

SubX= reshape(X,100,width,n) ;

%%  emprirical covariance

for fold = 1:n
    
    C_tr(:,:,fold)= cov(SubX(:,:,fold)');  %% True covariance.

for i = 1:100
    for j = 1:100
        
     C_num(i,j,fold)= mean(( SubX(i,:,fold)-mean(SubX(i,:,fold),2)).*  (SubX(j,:,fold)-mean(SubX(j,:,fold),2)));%-sum(X(i,:))*sum(X(j,:));
     
    end
end
end
  


%% MLE
for fold = 1:n
Mean_ml(:,fold) = mean(SubX(:,:,fold),2);


for i = 1:width
cov(:,:,i)= (SubX(:,i,fold)-Mean_ml(:,fold))*(SubX(:,i,fold)-Mean_ml(:,fold))';
end


C_ml(:,:,fold)= mean(cov,3);  % Maximum likelihood covariance.
end


%% diagonal composition model
for fold = 1:n
    
Co(:,:,fold) = diag(diag(C_ml(:,:,fold)));

alpha = 0.6;

C_super(:,:,fold)= (1-alpha)*Co(:,:,fold) +alpha*C_ml(:,:,fold);
end


%% Cross-validation


T = [2,1];

for fold = 1:n
    
    test = T(fold);
   sum1= mean(( C_tr(:,:,test)-C_super(:,:,fold)).^2,1);
    err_test(fold)= mean(sum1); 
    
end
Err_test= mean(err_test);


%% loglikelihood on test set


for fold = 1:n
    A = abs(C_ml(:,:,fold))^-1; % not real inverse. but the inverse of its elements.
    B = abs(C_super(:,:,fold))^-1;
    mu = Mean_ml(:,fold);

for i = 1:width
    
   test = T(fold); 
   ll_ml(i) =  -0.5*(SubX(:,i,test)-mu)'*A*(SubX(:,i,test)-mu);
   ll_super(i)= -0.5*(SubX(:,i,test)-mu)'*B*(SubX(:,i,test)-mu);
end
LL_super(fold)=( -100*log(sqrt(2*pi)) -0.5*log(det(A))  + sum(ll_super))/width;
LL_ml(fold)   =( -100*log(sqrt(2*pi)) -0.5*log(det(B))+ sum(ll_ml))/width;

end

LL_super_tt= mean(LL_super);
LL_ml_tt= mean(LL_ml);

alles = LL_super_tt-LL_ml_tt;




% 
% figure(1)
% surf(C_tr(:,:,1));
% figure(2)
% surf(C_ml(:,:,1)-C_tr(:,:,1));
% V0 = C_tr(1,4,1);
% V2 = ((width-1)/width)^-1*C_ml(1,4,1);
% V1 = C_num(1,4,1);
% 

%% My lienar regression model,




