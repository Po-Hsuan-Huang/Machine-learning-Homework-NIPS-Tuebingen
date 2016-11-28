%% machine learning ex3.

%% Task 2.5  compare model before and after PCA model A = [1 -1 ;1 1]/sqrt(2),
clear
% before PCA 
n = 1000;

u1 = -1 + 2*rand(1,n,1);
u2 = -1 + 2*rand(1,n,1);

%figure(1)
scatter(u1,u2,'.');

xlim([-2, 2]);
ylim([-2, 2]);

% after PCA

A =[1 -1 ;1 1]/sqrt(2);
U(1,:)= u1;
U(2,:)= u2;

for i = 1: n
V = A*U;
end
%figure(2)
scatter(V(1,:),V(2,:),'.');

xlim([-2, 2]);
ylim([-2, 2]);
%figure(3)
subplot(2,1,1)
Px1 = histogram(V(1,:),1000);
Vax1 = var(V(1,:));
subplot(2,1,2)
Px2 = histogram(V(2,:),1000);

Vax2 = var(V(2,:));

%%  our model: Maximum Likelihood Gaussian 

VarML = 1/3; 

meanML = 0;
% 
% W(1,:) = meanML+sqrt(VarML)*randn(1,n);
% W(2,:) = meanML+sqrt(VarML)*randn(1,n);

pd = makedist('Triangular','a',-sqrt(2),'b',0,'c',sqrt(2));

W(1,:) =trirnd(-sqrt(2),0,sqrt(2),n);
W(2,:) = trirnd(-sqrt(2),0,sqrt(2),n);


figure(4)
 scatter(V(1,:),V(2,:),'b.');
hold on
scatter(W(1,:),W(2,:),'r.');
xlim([-2, 2]);
ylim([-2, 2]);

hold off


%% estimate cross-emtropy.
figure(5)
sling= -sqrt(2):0.01:sqrt(2);
P_model= pdf(pd,sling);

plot(sling,P_model);


% x ,y >0
Pmx = (sqrt(2)-x)/2;
Pmy = (sqrt(2)-y)/2;

integral(log(P_model_model));














