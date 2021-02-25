function [chibest, chisq] = main_all_comb(gene, t0, keep)
% This codes fits the data with all possible models for each MURK gene

%%
data = dlmread(['./input/fit_s_u_',gene,'.txt']);


data = sortrows(data,1);



data(data(:,1)<t0,:) = [];

qq = mod(size(data,1),10);

data(1:qq,:) = [];


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




%% how many parameters

number_parameters = 9;

%% optimisation


start_base = zeros(1,number_parameters);

stop_base = ones(1,number_parameters) * 10340;

stop_base(3) = 0.8;
start_base(3) = 0.6;

if t0>=0.6
    
    start_base(3) = mtt(2);
    stop_base(3) = 0.89;
end

name_best = ['./bests/best_',gene,'.txt'];
data2 = dlmread(name_best);

err = data2(9);




temp = start_base(7:9);


start = [start_base(1:6),temp(find(keep))];


temp = stop_base(7:9);

stop = [stop_base(1:6),temp(find(keep))];

number_parameters = 6 + sum(keep);

options=optimset('TolX',1e-8,'TolFun',1e-8,'MaxFunEvals',1000,'MaxIter',6000);

name_best = ['./best_',gene,'_',num2str(keep(1)),num2str(keep(2)),num2str(keep(3)),'.txt'];


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

    dlmwrite(name_best,[ thetaRes(:);chisq ])
    
    
end


%% cost function
    function result=fitFun(theta)
        
        store = 6;
        
        alfa1 = theta(4);
        beta1 = theta(5);
        
        gam1 = theta(6);
        
        if keep(1)
        alfa2 = theta(store+1);
        store = store+1;
        else
            alfa2 = alfa1;
        end
        
        
        if keep(2)
        beta2 = theta(store+1);
        store = store+1;
        else
            beta2 = beta1;
        end
        
        
         if keep(3)
        gam2 = theta(store+1);
        else
            gam2 = gam1;
        end
        

        
 
        
        alfa = alfa1;
        beta = beta1;
        gam = gam1;
        
        
        I = theta(1:2);
        
        
        
        ts = theta(3);
        
        
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
            

%             
        end
        
        
        
        model = [model1;model2];
        



        result= (measured-model)./err;
            
        
        
    end



end
