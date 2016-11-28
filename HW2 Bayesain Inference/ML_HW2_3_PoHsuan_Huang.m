% Machine Learnign Exercise 2
% Problem 2(c)
% Po-Hsuan Huang 2014.11.7



% Plot posterior variance as function of n


function ML_HW2_3_PoHsuan_Huang(alpha, beta,data)

 
 %%posterior mean is simply the quotient of alpha/beta in the posterior
 %probability Gamma function. It is excatly the simplicity of calculation
 %mean that we should assume exponential family as our data distribution,
 %although it is not always the case.
 
 % cumulative sum. This function stores cumulative sum at each index of
 % a 1 by n array.
 summation = cumsum(data);
 
 length = 1:size(data,2);
 
 % posterior variance is simply alpha/beta
 posterior_var = (alpha+summation)./(beta+length).^2;
 
 %%plotting and printing results
 result1=posterior_var(10);
 fprintf('posterior variance when n=10 is %g.',result1);
 
 figure(1)
 plot(posterior_var);
 hold on
 legend('posterior var');
 title('posterior variance of theda over n');
 xlabel('number of data points n');
 ylabel('variance value of theda');
 hold off
end