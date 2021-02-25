function prediction_profile_likelihood(pop,tt)
% This code computed the 95% confidence bounds on the model for Smim1
% spliced (pop = 2) or unspliced (pop = 1) for each time point before the switch(tt between 1 and 18)

%% read data
measured_best=dlmread('smim1_mod.txt');

point_tt = measured_best(tt,pop+1);
time_tt = measured_best(tt,1);

steps = 100; % change this to a larger value if the profile isn't smooth

up = 1; % upper bound
down = 1; % lower bound


data2best = dlmread('best.txt');

chimin = data2best(end);
chibest = chimin;

nome_file = ['res_', num2str(pop),'_', num2str(tt)];

if exist(nome_file,'file')==2
    
    data=dlmread(nome_file);
    %
    [point_min,id_min] = min(data(:,1));
    [point_max,id_max] = max(data(:,1));
    
    best_min = data(id_min,2:12)';
    
    best_max = data(id_max,2:12)';
    
else
    point_min = point_tt;
    point_max = point_tt;
    
    best_min = data2best(1:11);
    best_max = data2best(1:11);
    
    
end


gene ='smim1';


data = dlmread(['fit_s_u_',gene,'.txt']);


data = sortrows(data,1);


t0 = 0.5;

data(data(:,1)<t0,:) = [];
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




measured_base = [mut',mst'];

t_end=0.9;



t_end=41;

%% how many parameters


number_parameters = 11;
%% optimization

max_counts = 1;
threshold = 3.8;

small_steps = 0.01;

options=optimset('display','iter','TolX',1e-8,'TolFun',1e-8,'MaxFunEvals',15000,'MaxIter',6000);%'display','iter',


start = zeros(1,number_parameters);
stop = ones(1,number_parameters)*10000;

stop(11) = time_tt;
start(11) = 0.6;

stop(9:10) = .03;

fName1 = nome_file;



if down
    point = point_min;
    best = best_min;
    
    
    
    stepd=abs(point)/steps;
    
    conta=0;
    
    if  abs(point)<0.001
        
        
        stepd=small_steps;
    end
    
    
    
    point_change = point - stepd
    
    
    
    
    
    while conta<max_counts && point_change >= 0
        
        
        [theta,chisq]=lsqnonlin(@fitFun,best,start,stop,options);
        
        
        
        
        if chisq>chibest+threshold
            conta=conta+1;
            
            
        else
            conta=0;
            best=theta;
            
            dlmwrite(fName1,[point_change;theta;chisq]','-append');
            
        end
        
        if chisq<chimin
            chimin=chisq;
            
            
        end
        
        
        point_change=point_change-stepd;
        
    end
    
end


if up
    
    point = point_max;
    best = best_max;
    
    stepu=abs(point)/steps;
    
    
    
    if  abs(point)<0.001
        
        
        stepu=small_steps;
    end
    
    
    
    conta=0;
    
    point_change=point+stepu
    
    
    
    
    
    
    while conta<max_counts
        
        
        [theta,chisq]=lsqnonlin(@fitFun,best,start,stop,options);
        
        
        if chisq>chibest+threshold
            conta=conta+1;
            
            
        else
            dlmwrite(fName1,[point_change;theta;chisq]','-append');
            
            
            conta=0;
            best=theta;
            
        end
        
        if chisq<chimin
            chimin=chisq;
            
            
        end
        
        
        point_change=point_change + stepu;
        
        
    end
    
    
    
end






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
        
        
        ts = theta(11);
        
        
        sol1 = ode45(@ODE,[t0,ts],I);
        
        
        model1 = deval(sol1,mtt(mtt<= ts))';
        
        
        
        
        alfa = alfa2;
        beta = beta2;
        gam = gam2;
        
        
        
        
        I = deval(sol1,ts);
        
        sol2=ode45(@ODE,[ts,t_end],I);
        
        model2 = deval(sol2,mtt(mtt>ts))';
        
        
        model_tt = deval(sol2,time_tt);
        model_point = model_tt(pop);
        
        
        
        function dxdt=ODE(~,x)
            dxdt=ones(2,1);
            
            
            
            
            dxdt(1) = alfa - beta * x(1); %unspliced
            
            dxdt(2) = beta * x(1) - gam * x(2);  %spliced
            
        
        end
        
        
        
        model = [model1;model2];
        
        
        error = [err_u' * ones(length(mtt),1), err_s' * ones(length(mtt),1)];
        
        
        
        
        result= sqrt(((measured_base-model)./error).^2 + 2  * log(error) +20);
        
        result = result(:);
        
        
        result(end+1) = (model_point-point_change)/0.005;
        
        
        
        
        
    end


end
