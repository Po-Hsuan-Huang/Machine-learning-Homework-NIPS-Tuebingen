%% Competition1 Single Gaussian Model

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

%% Getting parameter based on Subtrain

for fold = 1:n   %% how many subsets

    
mu(fold) = mean(Subtrain(fold,:));

Var(fold,:) = var( Subtrain(fold,:));

end

 %% Cross-validation. Use left-out segment 'Subset' as label.

for fold = 1:n % n folds
   [pdf4,c4,edge4]=histogram_pdf(Subset(fold,:),Subset(fold,:));
  
    edge4 = edge4(1,1:end-1)+diff(edge4,1)/2; 

Y(1,:)= (1/(sqrt(2*pi*Var(fold)))).*exp(-0.5*(edge4-mu(fold)).^2./Var(fold));

CV(fold) = sum( ((c4-Y(1,:)).^2).*Y(1,:)); %% use subset histogram as benchmark.
end

CV_tt = mean(CV);

%% getting estimated datapoint probability based on test set

for fold = 1: n
       
P(fold,:)= 1./(sqrt(2*pi*Var(fold))).*exp(-0.5*(Subset(fold,:)-mu(fold)).^2./Var(fold));

end

%% Compute loglikelilhood for each subset, and compare with histogram model.
for fold = 1:n

    
gamma(fold) = sum (log(P(fold,:)))/length(P(fold,:));  % nats/sample


end

GAMMA = mean(gamma);



