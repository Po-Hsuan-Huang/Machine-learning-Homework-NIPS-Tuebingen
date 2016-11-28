%% Histogram model based on X_train
%% second input is the training set,first input is the test set.
% output is the esimated value on testset based on training set.

function [pdf,val,edge]=histogram_pdf(X,Train) 

numbins = 60;
his = histogram(Train,numbins,'Normalization','pdf');
val = his.Values; 
edge = his.BinEdges;
width = his.BinWidth;

for i = 1 : length(X)

Binz  = ceil((X(i)-edge(1))/ width);

if Binz <=numbins && Binz > 0
   pdf(i) = val(Binz) ;
else
 display('test data outside function histogram_model define domain !!');   
end

  

end
end