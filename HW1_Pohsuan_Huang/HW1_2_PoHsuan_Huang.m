% Machine Learnign Exercise 1 
% Problem 3(d)
% Joint distribution of exponetial distribution
% Po-Hsuan Huang 2014.10.30 created.



%plot joint distribution
[x,y] = meshgrid(0:0.05:1,0:0.05:1);
q = (5*exp(-5*x)+2*exp(-2*y)); %Z = X+Y
r =  q*5*exp(-5*x);  %non-normalized probability X*Z
N = cumsum(r);       % calculate the cumulative sum of each column vectors.
M = cumsum(N(20,:)); % calculate the cumulative sum of the row vector.
S = M(1,20); % pick the last element, which is the sum of all elements of matrix r
R =r/S;  % normalized probability X*Z
figure(1)
subplot(2,1,1)
surface(x,y,R);
xlabel('x');
ylabel('y');
zlabel('X*Z');
title('joint distribution of X and Z');
axis tight;
view(45,45);
% plot marginalized distribution
subplot(2,1,2)    
D = cumsum(R'); % sum over all rows of x,equivalent to marginalization.
E = D(20,:);  % pick the last row, which is the cumulative sum of all x rows of each y.
plot(y,E,'b-',y,R(1,:),'r-.',y,R(20,:),'r--');
legend('marginalized','joint-distribution at x =0','joint-distribution at x =1');
xlabel('y');
ylabel('marginalized X*Z');




