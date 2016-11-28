%% Gaussian Model

%% second input is the training set,first input is the test set.
% output is the esimated value on testset based on training set.


function [pdf,y,tline]=Gaussian_pdf(X,Train)

mean_i = mean(Train);
var_i  = var(Train);
figure() 
tline = -4:0.1:4;

y= (1/sqrt(2*pi*var_i)) .* exp(-0.5*(tline-mean_i).^2./var_i);
plot(tline,y);

for i = 1:length(X)
    
pdf(i)= (1/sqrt(2*pi*var_i)) * exp(-0.5*(X(i)-mean_i)^2/var_i);

end

end
