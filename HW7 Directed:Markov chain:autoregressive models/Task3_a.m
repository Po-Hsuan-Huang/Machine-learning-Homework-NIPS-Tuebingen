%% Machine learning Lab bethege 2
%% Autoregressive model

%% This function use the whole xTrain(1,30000) as Training set
%% Thsis function use xTest(1:100) as Test phase
%% it only predicts the next time step

% one can adjust:





clear all;

steps = 100 ;% prediction period. 


load('training_data.mat')  % load given training data called X;
load('test_data.mat')
xTrain= X(1:30000);
xTest = Y (1:steps);


order = [0 1 2 4 8 16];

%% calculate the variance from the xTrain.

S = sqrt( var(xTrain) ) ;   % standart deviatin of the training set probability distribution.
 






for j = 1 :length(order)
    
    
    
    m = order(j);

%% train 

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


%% test
 %X_test= zeros(length(xTest)-m,m+1);
 predict = zeros(1,steps);
 predict(1,1:m) = xTest(1,1:m);
    
 clear X_test

%  X_test(1,1:m) = fliplr(  xTest(1,1:m) );
%  X_test(1,m+1) = 1; 

 for i= 1: steps-m

    X_test(i,1:m) = fliplr(  predict(1,i:i+m-1) );

    X_test(i,m+1) = 1;    


pre =  X_test(i,:) * Weight{j} ;  % predicted mean of order(j)

predict(1,m+i) =  pre + S.*randn(1,1);  % prediction with Gaussain Noice 


 end

Mean(1,1:steps,j) = predict;

end

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




% subplot(6,1,2)
% 
% xline = 1: 6+order(2);
% 
% plot(xline,xTest,xline,Mean(:,:,2));
% title('trend prediction_Gaussian m-th order')
% legend('real', '1st');
% xlabel('timesteps')
% ylabel('value')
% 
% subplot(6,1,3)
% 
% xline = 1: 6+order(3);
% 
% plot(xline,xTest,xline,Mean(:,:,3));
% title('trend prediction_Gaussian m-th order')
% legend('real', '2nd');
% xlabel('timesteps')
% ylabel('value')
% 
% subplot(6,1,4)
% xline = 1: 6+order(4);
% 
% plot(xline,xTest,xline,Mean(:,:,4));
% title('trend prediction_Gaussian m-th order')
% legend('real', '4th');
% xlabel('timesteps')
% ylabel('value')
% 
% subplot(6,1,5)
% xline = 1: 6+order(5);
% 
% plot(xline,xTest,xline,Mean(:,:,5));
% title('trend prediction_Gaussian m-th order')
% legend('real', '8th');
% xlabel('timesteps')
% ylabel('value')
% 
% subplot(6,1,6)
% xline = 1: 6+order(6);
% 
% plot(xline,xTest,xline,Mean(:,:,6));
% title('trend prediction_Gaussian m-th order')
% legend('real', '16th');
% xlabel('timesteps')
% ylabel('value')
% 
% 
