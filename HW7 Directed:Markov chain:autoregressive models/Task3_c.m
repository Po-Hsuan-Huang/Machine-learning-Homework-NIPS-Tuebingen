%% Machine learning Lab bethege 2

%% Autoregressive model

%% Task 3(c) With historoy dependent varaince Gaussian pdf.

% Our model assuems that var = W*X_km +C ;

%% This function use the whole xTrain(1,30000) as Training set
%% Thsis function use xTest(1:100) as Test phase
%% it only predicts the next time step

% one can adjust:





clear all;

steps = 100 ;% prediction period. 


load('training_data.mat')  % load given training data called X;
load('test_data.mat')
xTrain= X(1:30);
xTest = Y (1:steps);


order = [0 1 2 4 8 16];

for j = 1 :length(order)
    
    
    
    m = order(j);

%% train weight of xTrain

 %  X_km= zeros(length(xTrain)-m,m+1); 
     clear X_km

for i = 1 : length(xTrain)-m
    X_km(i,1:m) = fliplr(  xTrain(1,i:i+m-1) );
end
    X_km(:,m+1) = ones(length(xTrain)-m,1);    
    

t= xTrain(1,m+1:end)';
% dude=(X_km'*X_km);
% sis = X_km'\(X_km'*X_km);
W =(  X_km'\(X_km'*X_km) )' * t;
Weight{j}= W/sum(W);

b(j) = W(end,1);  % offset constant




end

%% Prediction error for each xTrain point.

% to generate eTrain, (pre-true)^2 training set,  to train linear model of
% histroy dependent variance:  variance =  V'*Xkm+c.
% Here we assume variance is squar error of each datapoint.

for j = 1 : length(order)
    
    m = order(j);
    
    % design matrix
           clear Xkm

    for i = 1 : length(xTrain)-m
         Xkm(i,1:m) = fliplr(xTrain(1,i:i+m-1)); % (N-M)xM matrix
    end
         Xkm(:,m+1) = ones(length(xTrain)-m,1);% (N-M)x(M+1) matrix

        pre(j,1:m) = xTrain(1,1:m);

        pre(j,m+1:length(xTrain)) = (Xkm  * Weight{j})'; 
        
        eTrain(j,:) = (xTrain-pre(j,:)).^2; 
end
%% Train Weight of eTrain


 %  X_km= zeros(length(xTrain)-m,m+1); 
     clear X_km   % our model is  variance =  V'*Xkm+c 
for j = 1 : length(order)
    
    m = order(j);
    clear X_km
    
for i = 1 : length(xTrain)-m
    X_km(i,1:m) = fliplr(  xTrain(1,i:i+m-1) );
end
    X_km(:,m+1) = ones(length(xTrain)-m,1);  % the last term of Weight V is c.    
    

t= eTrain(j,m+1:end)';
% dude=(X_km'*X_km);
% sis = X_km'\(X_km'*X_km);
V =(  X_km'\(X_km'*X_km) )' * t;
Weight_v{j}= V/sum(V);

end


%% test

for j = 1 : length(order)
    
  m = order(j);
 %X_test= zeros(length(xTest)-m,m+1);
 predict = zeros(1,steps);
 predict(1,1:m) = xTrain(1,1:m);
    
 clear X_test

%  X_test(1,1:m) = fliplr(  xTest(1,1:m) );
%  X_test(1,m+1) = 1; 

 for i= 1: steps-m

    X_test(i,1:m) = fliplr(  predict(1,i:i+m-1) );

    X_test(i,m+1) = 1;    


pre =  X_test(i,:) * Weight{j} ;  % predicted mean of order(j)
pre_var = X_test(i,:) * Weight_v{j}; % predicted variance of order(j)
predict(1,m+i) =  pre + sqrt(pre_var)*randn(1,1);  % prediction with predictive vairance. 


 end

Mean(1,1:steps,j) = predict;

end

%% plot 

figure (1)


xline = 1: steps;
plot(xline,xTest,xline,Mean(:,:,1),xline,Mean(:,:,2),xline,Mean(:,:,4),xline,Mean(:,:,4),xline,Mean(:,:,5),xline,Mean(:,:,6));
title('trend prediction Gaussian m-th order')
legend('real', '0th','1st','2nd','4th','8th','16th');
xlabel('timesteps')
ylabel('value')
%%
figure(2)

for i = 1:length(order)  % plot m subplots for each order of Gaussian process 
    
    
subplot(6,1,i)


xline = 1: 6+order(i) ;
u(1,:) = Mean(:,:,i) ;
% plot(xline,u);
% hold on
% plot(xTest(1: steps));
% hold off
plot(xline,xTest(1,1: 6+order(i)),xline,u(1,1: 6+order(i)));
Str = sprintf('trend prediction Gaussian %d-th order', order(i));
Str2 = sprintf('blue is real, red is prediction');
title({Str,Str2});
%legend('real', 'predict');
xlabel('timesteps')
ylabel('value')
axis tight;

end 









% 
% figure (1)
% 
% 
% xline = 1: steps;
% plot(xline,xTest,xline,Mean(:,:,1),xline,Mean(:,:,2),xline,Mean(:,:,4),xline,Mean(:,:,4),xline,Mean(:,:,5),xline,Mean(:,:,6));
% title('trend prediction_Gaussian m-th order')
% legend('real', '0th','1st','2nd','4th','8th','16th');
% xlabel('timesteps')
% ylabel('value')
% 
% figure(2)
% 
% for i = 1:length(order)  % plot m subplots for each order of Gaussian process 
%     
%     
% subplot(6,1,i)
% 
% 
% xline = 1: 6+order(i);
% u = Mean(:,:,i);
% plot(xline,xTest(1: 6+order(i)),xline,u(1: 6+order(i)));
% title('trend prediction_Gaussian m-th order')
% %legend('real', 'predict');
% xlabel('timesteps')
% ylabel('value')
% 
% 
% end 

