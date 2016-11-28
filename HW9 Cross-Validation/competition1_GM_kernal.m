%% Machine learning, Bethege 4.  Competition1.

%% variables discription
% mean log-likelihood of Kernel method is GM_Ker.
% mean cross-validation test error is CV_tt. 


clear

load('X_train.mat')  % load given training data called X;


%% pdf on Training set
[pdf1, c1,edge1] =Gaussian_pdf(X_train(1,:),X_train(1,:));
[pdf2, c2,edge2]= histogram_pdf(X_train(1,:),X_train(1,:));



%%  Linear regression with Gaussian Kernel:
% use the biggest 4 modes to construct the model. 8 parameters in total.

n = 10; % n fold corss-validation
d = 50; % number of Gaussian components. It had better be more than 50.
width= length(X_train)/n;
edge2 = edge2(1,1:end-1)+diff(edge2,1)/2;
[psor,lsor] = findpeaks(c2,edge2,'MinPeakProminence',0.1,'SortStr','descend'); % use peaks as initial mean guess.

   

%%  produce subtraining set 'Subtrain'.
% 
% for i = 1:n 
%    
% Subset(i,1:width) =  X_train(1,width*(i-1)+1:width*i);  % chunk trianingset into n segment
% end


Subset= (reshape(X_train,width,n))';
for fold = 1:n
    
Seg{fold} = Subset(fold,:) ;  

end
 
 
%    regroup= 99*ones(n, length(X_train));
%     Subtrain= 99*ones(n,length(X_train)*((n-1)/n));
%    for i = 1 :n
%   
%    for j = 1:n
%    if j ~= i    
%    regroup(i,1:j*width) =  cat(2, regroup(i,1:(j-1)*width), Subset(j,:)); % concacnate segments
%    end
%    end
%       Subtrain(i,:)=regroup(i,(regroup(i,:)~=99)); % kickout unreplaced '99' values   
% 
%    end
  

 for fold = 1:n
 tag = ones(1,n); tag(fold) = 0;   
 Subtrain(fold,:)= [Seg{tag==1}];   %% I am sooo smart!
  
 end


%% use alpha = beta= 1, calculate the posterior mean

sigma= 0.5;

for fold =1 :n
[pdf3, c3,edge3]= histogram_pdf(Subtrain(fold,:),Subtrain(fold,:));
 edge3 = edge3(1,1:end-1)+diff(edge3,1)/2;

feature= zeros(d,length(edge3));
for i = 1:d
    feature(1,:)= edge3;
    xlin = linspace(-4,4,d);
    feature(i+1,:)=  (sqrt(2*pi*sigma^2))^-1.*exp(-0.5*(edge3-xlin(i)).^2/sigma^2);
    
end
%% Calculate the Nx50 matrix ztrain for which the n-th row is  Z(Xn) and produce 

zTrain = feature(2:d+1,:)';



alpha = 1;
beta = 1; 

Covpost= inv(alpha*eye(d)+beta*(zTrain'*zTrain));  % covariance matrix
Meal= beta*Covpost*zTrain'*c3';             % posterior mean


Meanpost(fold,:)= Meal'./sum(Meal);             % normalization


%% Q2_3 weighted predictive mean using xspace
  
  Model(fold,:)= zeros(1,length(edge3));

 for i=1:d   
    Model(fold,:)=Model(fold,:) + Meanpost(fold,i).*(sqrt(2*pi*sigma^2))^-1.*exp(-0.5*(edge3-xlin(i)).^2/sigma^2);    
 end   
 


end
%% getting estimated  probability on testset for each fold
for fold = 1: n
        
% for i = 1 : length(Subset(fold,:))
% 
% Binz  = ceil((Subset(fold,:)-edge3(1))./ diff(edge3(1:2)));
% 
% check1= Binz <= length(edge3); 
% check2 = Binz > 0;
% check3= check1.*check2;
% if prod(check3)==0
%   str1= sprintf('%d th test set has datapoint outside difine domain of %d fold predictive function ',fold,fold);
%     display(str1);
% end
% 
% 
% j=1;
% for i  =  1 : length(check3)
% if check3(i)==1
% Bin(j) = Binz(i);
% j= j+1;
% end
% end
Meanpredict(fold,:)= zeros(1,length(Subset(fold,:)));

 for i=1:d   
    Meanpredict(fold,:)= Meanpredict(fold,:) + Meanpost(fold,i).*(sqrt(2*pi*sigma^2))^-1.*exp(-0.5*(Subset(fold,:)-xlin(i)).^2/sigma^2);    
 end   
 



%% Compute loglikelilhood for each subset, and compare with histogram model.

gmk(fold) = sum( log(Meanpredict(fold,:)))/length(Subset(fold,:));  

  

end

GM_Ker = mean(gmk);





%% Cross-validation. Use left-out segment as label. yield 10 CV values
for fold = 1:n % n folds
   [pdf4,c4,edge4]=histogram_pdf(Subset(fold,:),Subset(fold,:));
    edge4 = edge4(1,1:end-1)+diff(edge4,1)/2; 


CV(fold) = sum( ((c4-Model(fold,:)).^2).*Model(fold,:));
end

CV_tt = mean(CV);


%% Plot Kernel pdf and well-chosen test-set histogram model
% choose one fold to plot
fold = 4; 
[pdf5,c5,edge5]=histogram_pdf(Subset(fold,:),Subset(fold,:));

figure(5)
[hAx,hLine1,hLine2] = plotyy(edge4,c5,edge3,Model(fold,:),'bar','plot');
hLine2.Color = 'red';
hAx(1).YLim= [0 1];
hAx(2).YLim= [0 1];%
hAx(1).YTick = 0: 0.1: 1;
hAx(1).YGrid = 'on';
legend('True','GMM');
str1 = sprintf('test set %d hist & %d th subtrain Kernel',fold,fold);
title({str1});



% comparison between Gaussian and GMM and Histogram
figure(6)
y1=Model(fold,:);
xlin= [edge3,edge1];
ylin = [y1,c1];
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





