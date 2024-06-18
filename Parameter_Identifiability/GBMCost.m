% Cost Function for the GBM model
% Original Source: Marisa Eisenberg (marisae@umich.edu)
% Since we had less T cell data than the other cell types, we had to modify how the T cell data 
% was inputted into the error equations. We used the weighted ordinary least squares (OLS) error, 
% but this script constains multiple other error functions that can be used: Poisson ML, 
% root mean square error (RMSE), mean square error (MSE), and relative error.

function NLL = GBMCost(tspan,fixed,params,indexchooseparam,data,initialcondition)
params = abs(params);

[t,x] = ode45(@GBMFuncidentifiable,[0 tspan(:,:,1)'],initialcondition,[],fixed,params,indexchooseparam);

Tcelldata = data(:,:,2);
Tcelldata = Tcelldata(~isnan(Tcelldata))';
indices=zeros(length(t),1);
for i = 1:length(t)
    for j = 1:length(Tcelldata)
        if t(i) == tspan(j,:,2)
            indices(i) = 1;
        end
    end
end
xTcell = nonzeros(x(:,2).*indices)';

% ERROR FUNCTIONS (NLL = ...)

% Poisson ML
%tumorerror = (x(2:end,1))'*ones(length(x(2:end,1)),1) - data(:,:,1)'*log(x(2:end,1));
%Tcellerror = (xTcell)*ones(length(xTcell),1) - Tcelldata*log(xTcell)';
%MDSCerror = (x(2:end,3))'*ones(length(x(2:end,3)),1) - data(:,:,3)'*log(x(2:end,3));
%NLL = abs(tumorerror) + abs(Tcellerror) + abs(MDSCerror);  


% Root mean square error
%NLL = (sqrt(sum((data(:,:,1) - x(2:end,1)).^2)/length(tspan(:,:,1)))) + (sqrt(sum((Tcelldata - xTcell).^2)/length(Tcelldata))) + (sqrt(sum((data(:,:,3) - x(2:end,3)).^2)/length(tspan(:,:,3)))); 


% Mean square error (MSE) 
%NLL = (sum((data(:,:,1) - x(2:end,1)).^2))/length(tspan(:,:,1)) + (sum((Tcelldata - xTcell).^2))/length(Tcelldata) + (sum((data(:,:,3) - x(2:end,3)).^2))/length(tspan(:,:,3));


% Relative error
% NLL = ((1/length(data(2:end,:,1))*sum((abs(data(2:end,:,1) - x(2:end,1)))./data(2:end,:,1))) ...
% + (1/length(data(2:end,:,2))*sum((abs(data(2:end,:,2) - x(2:end,2)))./data(2:end,:,2)))...
%    + (1/length(data(2:end,:,3))*sum((abs(data(2:end,:,3) - x(2:end,3)))./data(2:end,:,3))))/(length(tspan)-1); 
    

% Ordinary least squares (OLS) -- weighted
NLL = ((sum((data(:,:,1) - x(2:end,1)).^2)) + (6/4)*(sum((Tcelldata - xTcell).^2)) + (sum((data(:,:,3) - x(2:end,3)).^2)))/3; % essentially like MSE, but MSE is like weight version of OLS (just with predetermined weights)

