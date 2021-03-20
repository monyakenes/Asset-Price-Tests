clc
clear
clf

% 2005-01-01 ~ 2019-12-31
Datatable1 = readtable('Amazon.xlsx','ReadRowNames',true);
Datatable2 = readtable('IBM.xlsx','ReadRowNames',true);
Datatable3 = readtable('JPM.xlsx','ReadRowNames',true);

Num = 500;  %% Write the number of samples that desires to choose, e.g., 250, 500, 1000 
Data(:, 1) = Datatable1{2:Num,{'PriceOrBid_AskAverage'}}
Data(:, 2) = Datatable2{2:Num,{'PriceOrBid_AskAverage'}}
Data(:, 3) = Datatable3{2:Num,{'PriceOrBid_AskAverage'}}

%% Return
for i=1:3
    for l=1:Num-2
        Data_return(l, i) = (Data(l+1, i)/Data(l, i)) - 1;
    end
end

%% Normalization
for i=1:3
    mu(i) = mean(Data_return(:, i));
    std(i) = sqrt(var(Data_return(:, i)));
    
    temp = Data_return(:, i) - mu(i);
    Y(:, i) = temp.*(1/std(i));
end


%%% Kernel density estimation
dx = 0.05;
centers = [-5: dx: 5];
M = length(Data_return(:, 1));

for i=1:3
    X = Y(:, i);
    N = hist(X, centers)
    
    figure(2*i - 1) 
    bar(centers, N/(M * dx), 'blue')
    hold on
    plot(centers, 1/sqrt(2*pi) * exp(-0.5*centers.^2),'r','LineWidth',1.5)     
    legend('Empirical Distribution', 'Standard Normal')
    title('Kernel density estimation')

    %%% Quantile of normal distribution
    X_inc = sort(X);
    for k = 1:M
        p(k) = k/(1+M);
        Theoretical(k) = norminv(p(k), 0 , 1);
    end
    figure(2*i) 
    scatter(X_inc, Theoretical, 'cyan')
    hold on
    plot(centers, centers,'r','LineWidth',1.5)
    legend('Pairs of sorted data and quantile points', 'y=x line')
    title('QQ plot')

end