% Machine Learnign Exercise 2
% Problem 2(b)
% Po-Hsuan Huang 2014.11.7


% calculate btoh the Maximum Likelihood Estimation and the posterior mean
% using the first n points.


function ML_HW2_2_PoHsuan_Huang(alpha, beta,data)

 
 %%posterior mean is simply the quotient of alpha/beta in the posterior
 %probability Gamma function. It is excatly the simplicity of calculation
 %mean that we should assume exponential family as our data distribution,
 %although it is not always the case.
 
 % cumulative sum. This function stores cumulative sum at each index of
 % a 1 by n array.
 summation = cumsum(data);
 
 length = 1:size(data,2);
 
 % posterior mean is simply alpha/beta
 posterior_mean = (alpha+summation)./(beta+length);
 
 %% Maximal Likelihood Estimation is simply the mean of data, which is proved in question 1.a
 
 %  MLE is simply summation/number
 MLE = (summation)./(length);
 
 %% Plotting two functions, and print reusults when n=10.
 result1=posterior_mean(10);
 fprintf('posterior mean when n=10 is %g.',result1);
 result2=MLE(10);
 fprintf('MLE when n=10 is %g.',result2);
 
 figure(1)

 plot(posterior_mean);
 hold on
 plot(MLE);
 legend('posterior mean','MLE');
 title('posterior mean and MLE for theda over n');
 xlabel('number of data points n');
 ylabel('predictive value of theda');
 hold off

end




