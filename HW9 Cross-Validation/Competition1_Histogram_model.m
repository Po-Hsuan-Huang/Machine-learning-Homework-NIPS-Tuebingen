
%% Machine learning, Bethege 4.  Competition1.Histogram model
%% variables discription
% mean log-likelihood of Kernel method is DELTA.
% mean cross-validation test error is CV_tt.

% this function use external function histogram_pdf,which is included in the
% folder.
%% Competition1 histogram model


clear

load('X_train.mat')  % load given training data called X;

%% Singel Gaussian Model:

n = 10; % n fold corss-validation

width= length(X_train)/n;
    
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



 %% Cross-validation. Use left-out segment 'Subset' as label.

for fold = 1:n % n folds
   [pdf4,c4,edge4]=histogram_pdf(Subset(fold,:),Subset(fold,:));
  
    edge4 = edge4(1,1:end-1)+diff(edge4,1)/2; 
    
    [pdf1, c1, edge1]= histogram_pdf(Subset(fold,:),Subtrain(fold,:));    

    CV(fold) = sum( ((c4-c1).^2).*c1); %% use subset histogram as benchmark.
end

CV_tt = mean(CV);


%% Compute loglikelilhood for each subset, and compare with histogram model.
for fold = 1:n

    
delta(fold) = sum (log(pdf1))/length(pdf1);  % nats/sample


end

DELTA = mean(delta);



