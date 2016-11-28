%% Machine learning, Bethege 4. .



%% competition1
clear

load('X_train.mat')  % load given training data called X;


[pdf1, c1,edge1] =Gaussian_pdf(X_train);
[pdf2, c2,edge2]= histogram1(X_train);

% figure(1)
% bar( edge2,c2);
% hold on
% plot(edge2, val);
% hold off
gamma = 0;
Delda = 0;
for i = 1 : length(pdf1) 
gamma= gamma+log(pdf1(i));
Delda =Delda+log(pdf2(i));
end
ans1= (Delda - gamma)/length(pdf1);
%% Gausssina Kernel

%% Mixture of Gaussian Model:
% use the biggest 4 modes to construct the model. 8 parameters in total.

n = 10; % n fold corss-validation
d =  10; % number of Gaussian components.
width= length(X_train)/n;
   

 % find 3 modes of Gaussian mixture
     edge2 = edge2(1,1:end-1)+diff(edge2,1)/2;

     [psor,lsor] = findpeaks(c2,edge2,'SortStr','descend');

     



%% plot leave one out subsets
for i = 1:n 
   
   Subset(i,1:width) =  X_train(1,width*(i-1)+1:width*i);
   [pdfsub,Values,Edges]= histogram1(Subset(i,:));
   Subval(i,:)= Values;
end
%% plotting each subtraining set 
% 
    edge3 = Edges;
    edge3 = edge3(1,1:end-1)+diff(edge3,1)/2;
    regroup= 99*ones(n, length(X_train));
    Subtrain= 99*ones(n,length(X_train)*((n-1)/n));
   for i = 1 :n
   Subhist(i,:) = sum(Subval)-Subval(i);
  
   for j = 1:n
   if j ~= i    
   regroup(i,1:j*width) =  cat(2, regroup(i,1:(j-1)*width), Subset(j,:));
   end
     end
      Subtrain(i,:)=regroup(i,(regroup(i,:)~=99));

%    figure (i) 
%    bar(edge3,S(i,:),1);  %% choose bar width =1 (full width)
%    str1= sprintf('leave %d subset out', i);
%    title({str1});
   end
 


%% EM for Gaussian Mixture Model.

% initial guess
for fold = 1:n   %% how many subsets

mu = lsor; % the mean in descending order.

var = 10*rand(1,length(mu)); 

wei = psor./sum(psor);

%%






for iter = 1:3  %% EM iteration

% E step
    
    
    
    for j = 1:d
     F(j,:) = wei(j).*(1/(sqrt(2*pi*var(j)))).*exp(-(Subtrain(fold,:)-mu(j)).^2./var(j));
    end
    
    for j = 1:d
     Z(1,:)= sum(F,1);
     R(j,:) = F(j,:)./ Z(1,:);   %% rik responsibility: how much percentage of data point i is k mode responsible for
     
     G(j) = sum(R(j,:),2);       %%  how much percentage of data is K mode responsible for 
    end
    
% M step


for j = 1:d
wei(j) = G(j)/(length(X_train)*((n-1)/n));

mu(j) = sum(R(j,:).* Subtrain(i,:),2)./G(j);

var(j) = sum( R(j,:).*Subtrain(fold,:).*Subtrain(fold,:)./G(j),2);
end


end



end
%% estimate pdf
for i = 1: length(X_train)
for j = 1:d
P(j,i)= wei(j).*(1/(sqrt(2*pi*var(j)))).*exp(-(X_train(i)-mu(j)).^2./var(j));
end
end
Z= sum(P,1);


% figure()
% [hAx,hLine1,hLine2] = plotyy(edge3,c2,linex,Z(1,:),'bar','plot')
% hLine2.Color = 'red'
% hAx(1).YLim= [0 1];
% hAx(2).YLim= [0 1];% 
% hAx(1).YTick = 0: 0.1: 1;
% hAx(1).YGrid = 'on'
% legend('True','GMM')

% plot(linex, sum(P,1),'color','b');
% hold on
% bar()
% 
% hold off


%%


GMM= 0;
for i = 1 : length(Z) 
GMM= GMM+log(Z(i));
end
ans2= (GMM-gamma)/length(Subtrain(1,:));


