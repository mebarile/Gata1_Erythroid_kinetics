function main_find_best_fit
% This code finds the best fit for Smim1 spliced and unspliced counts

%%
fileID = fopen('gene_name.txt','r');
gene = fscanf(fileID,'%s');
fclose(fileID);



data = dlmread(['./fit_s_u_',gene,'.txt']);


data = sortrows(data,1);



t0 = 0.5;

data(data(:,1)<t0,:) = [];


size(data)

data(1:6,:) = [];

n = size(data,1);

measured=data(:,3:-1:2);

m = max(measured);

measured = measured./repmat(m,n,1);

%
st = reshape(measured(:,2),n/10,10);
ut = reshape(measured(:,1),n/10,10);



mst = mean(st);
mut = mean(ut);



mst = mst/max(mst);
mut = mut/max(mut);



tt = reshape(data(:,1),n/10,10);
mtt = mean(tt);

t0 = mtt(1);




measured = [mut',mst'];
% error = [eu',es'];

t_end=0.9;



t_plot = t0:0.01:t_end;

%% how many parameters

number_parameters = 10;

%% optimisation


start = zeros(1,number_parameters);

stop = ones(1,number_parameters) * 1000;

stop(10) = 0.8;

stop(9) = .2;



options=optimset('Display','Iter','TolX',1e-8,'TolFun',1e-8,'MaxFunEvals',1000,'MaxIter',6000);

name_best = ['./best_',gene,'.txt'];

if exist(name_best,'file')==2
    
    data2 = dlmread(name_best);
    chibest = data2(end);
    
    guess = data2(1:number_parameters);
    
    
    
else
    
    guess = rand(1,number_parameters)*10;
    
    
    
    
    chibest=exp(20);
    
    
end



[thetaRes,chisq] = lsqnonlin(@fitFun,guess(:),start,stop,options);


if chisq<chibest
    dlmwrite(name_best,[thetaRes(:);chisq ])
    
    
end



%% cost function
    function result=fitFun(theta)
        
        alfa1 = theta(1);
        alfa2 = theta(2);
        
        beta1 = theta(3);
        beta2 = theta(4);
        
        gam1 = theta(5);
        gam2 = theta(6);
        
        
        
        alfa = alfa1;
        beta = beta1;
        gam = gam1;
        
        
        I = theta(7:8);
        
        err_u = theta(9);
        err_s = theta(9);
        
        
        ts = theta(10);
        
        
        sol1 = ode45(@ODE,[t0,ts],I);
        
        
        model1 = deval(sol1,mtt(mtt<= ts))';
        
        alfa = alfa2;
        beta = beta2;
        gam = gam2;
        
        
        
        
        I = deval(sol1,ts);
        
        sol2=ode45(@ODE,[ts,t_end],I);
        
        model2 = deval(sol2,mtt(mtt>ts))';
        
        
        
        
        function dxdt=ODE(~,x)
            
            dxdt=ones(2,1);
            
            
            dxdt(1) = alfa - beta * x(1); %unspliced
            
            dxdt(2) = beta * x(1) - gam * x(2);  %spliced
            
            
            
        end
        
        
        
        model = [model1;model2];
        
        
        error = [err_u' * ones(length(mtt),1), err_s' * ones(length(mtt),1)];
        
        
        result= sqrt(((measured-model)./error).^2 + 2  * log(error) +20);
        
        
        
        
    end



end