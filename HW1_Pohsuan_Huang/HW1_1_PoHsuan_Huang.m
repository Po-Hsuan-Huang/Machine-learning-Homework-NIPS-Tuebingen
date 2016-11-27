% Machine Learnign Exercise 1 
% Problem 2(b)
% conditional probability of medical tests.
% Po-Hsuan Huang 2014.10.29 created.

function HW1_1_PoHsuan_Huang()
% plot a diagram showing the relation between P(deseased)
% and P(deseased|positive ).
figure(1);
a = @D1T1;
subplot(2,1,1);
fplot(a,[0,0.1]);
%legend('');
ylabel('P(deseased|positive )');
xlabel(' P(deseased)');
title('P(deseased|positive ) distribution over P(diseased)');


% plot a diagram showing the relation between P(deseased)
% and P(healthy|negative ).
subplot(2,1,2);
b= @D0T0;
fplot(b,[0.9,1]);
%legend('');
ylabel(' P(healthy|negative)');
xlabel(' P(deseased)');
title(' P(healthy|negative) distribution over P(diseased)');


% if the prevalence of the desease is 0.1%, calculating the 
% probability of deseased if test is positive.
answer = D1T1(0.001);
fprintf(' deseased probability if test is positive is %g',answer);
% ans = 0.5000.

end

% Construct a function D0T0 for conditional probablity P(healthy|negative )
function y = D0T0 (d) % d denotes P(deseased)

% set up a matrix A indicating the probability of different condisitons
A = [0.999,0.001,0.999,0.001];
% the first element is P(negative|healthy);
% the 2nd element is P(positive|healthy)
% the 3rd element is P(positive|deseased);
% the 4th element is P(negative|deseased)
y = (A(1)*(1-d))/ (A(1)*(1-d)+ A(4)*d);

end

% Construct a function D1T1 for conditional probablity P(deseased|positive )
function [y] = D1T1(d)

A = [0.999,0.001,0.999,0.001];
y = (A(3)*d)/ (A(2)*(1-d)+ A(3)*d);

end


% It may be better to draw in loglog plot, but I don't know hwo to draw it.