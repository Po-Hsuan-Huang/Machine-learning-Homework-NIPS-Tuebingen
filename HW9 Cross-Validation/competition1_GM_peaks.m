%% Machine learning, Bethege 4.  Competition1.

%% variables discription
% mean log-likelihood of Kernel method is GMM_PK.
% mean cross-validation test error is CV_tt. 

%% This program does following tasks:
% 1. get the loglikelihood difference between Histogram model and Gaussian Model.
%    which set a bench mark for comparison.

% 2. chunking raw trainingset into k folds


% 3.Use finpeaks to get GM mean and weight, adjust variance manually.


% 4. Cross-validation for each fold, with its mean 'CV_tt'.

% 5. calculate loglikelihood for each fold, with its mean 'LOGLIKELIHOOD'.

% 6. plot out the Gaussian, GMM, and banchmark Histogram on one graph.

% Basically, THE RESULT IS VERY BAD.


clear

load('X_train.mat')  % load given training data called X;


%% pdf on Training set
[pdf1, c1,edge1] =Gaussian_pdf(X_train(1,:),X_train(1,:));
[pdf2, c2,edge2]= histogram_pdf(X_train(1,:),X_train(1,:));



%%  Gaussian Mixture Model:
% use the biggest 4 modes to construct the model. 8 parameters in total.

n = 10; % n fold corss-validation
d = 5; % number of Gaussian components.
width= length(X_train)/n;
edge2 = edge2(1,1:end-1)+diff(edge2,1)/2;
[psor,lsor] = findpeaks(c2,edge2,'MinPeakProminence',0.1,'SortStr','descend'); % use peaks as initial mean guess.

     




   

%%  produce subtraining set 'Subtrain'.
% 
for i = 1:n 
   
Subset(i,1:width) =  X_train(1,width*(i-1)+1:width*i);  % chunk trianingset into n segment
end

   
    regroup= 99*ones(n, length(X_train));
    Subtrain= 99*ones(n,length(X_train)*((n-1)/n));
   for i = 1 :n
  
   for j = 1:n
   if j ~= i    
   regroup(i,1:j*width) =  cat(2, regroup(i,1:(j-1)*width), Subset(j,:)); % concacnate segments
   end
   end
      Subtrain(i,:)=regroup(i,(regroup(i,:)~=99)); % kickout unreplaced '99' values   

   end


%% Findpeaks for Gaussian Mixture Model.

for fold = 1:n
    
mu(fold,:) = lsor(1:d); % the mean in descending order.

var(fold,:) = (2*pi)^-1*psor.^-2; % use peakheight to infer variance

wei(fold,:) = psor(1:d)./sum(psor(1:d));

end

%% getting estimated  probability on testset for each fold
clear Z
for fold = 1: n
        
for j = 1:d
P(j,:)= wei(fold,j).*(1/(sqrt(2*pi*var(fold,j)))).*exp(-0.5*(Subset(fold,:)-mu(fold,j)).^2./var(fold,j));
end

Z(fold,:)= sum(P,1);

end

%% Compute loglikelilhood for each subset, and compare with histogram model.
for fold = 1:n
for i = 1 : length(Z(fold,:)) 
gmm(fold) = sum(log(Z(fold,:)))/length(Subset(fold,:)); % nats/sample
end
end
GMM_PK = mean(gmm);

%% Cross-validation. Use left-out segment as label. yield 10 CV values
for fold = 1:n % n folds
   [pdf4,c4,edge4]=histogram_pdf(Subset(fold,:),Subset(fold,:));
    edge4 = edge4(1,1:end-1)+diff(edge4,1)/2; 

for j = 1:d
Y(j,:)= wei(fold,j).*(1/(sqrt(2*pi*var(fold,j)))).*exp(-0.5*(edge4-mu(fold,j)).^2./var(fold,j));
end
Y_tt= sum(Y,1);
CV(fold) = sum( ((c4-Y_tt).^2).*Y_tt);
end

CV_tt = mean(CV);


%% Plot GMM pdf and well-chosen test-set histogram model
% choose one fold to plot
fold = 4; 
[pdf5,c5,edge5]=histogram_pdf(Subset(fold,:),Subset(fold,:));
linex =-4:0.1:4; 
clear P
for j = 1:d
P(j,:)= wei(fold,j).*(1/(sqrt(2*pi*var(fold,j)))).*exp(-0.5*(linex-mu(fold,j)).^2./var(fold,j));
end


figure(5)
[hAx,hLine1,hLine2] = plotyy(edge4,c5,linex,sum(P,1),'bar','plot');
hLine2.Color = 'red';
hAx(1).YLim= [0 1];
hAx(2).YLim= [0 1];%
hAx(1).YTick = 0: 0.1: 1;
hAx(1).YGrid = 'on';
legend('True','GMM');
str1 = sprintf('test set %d hist & %d th subtrain GMM',fold,fold);
title({str1});



% comparison between Gaussian and GMM and Histogram
figure(6)

xlin= [linex,edge1];
ylin = [sum(P,1),c1];
[hAx,hLine1,hLine2] = plotyy(edge2,c2,xlin,ylin,'bar','plot');
% set(hLine2(1), 'Color', 'r'); % for example...salt to suit...
% set(hLine2(2), 'Color', 'b'); % for example...salt to suit...

hAx(1).YLim= [0 1];
hAx(2).YLim= [0 1];
hAx(1).YTick = 0: 0.1: 1;
hAx(1).YGrid = 'on';
legend('hist','GMM', 'Gaussian');
str1 = sprintf('Xtrain hist. & %d th fold GMM & Gaussian ',fold);
title({str1});





