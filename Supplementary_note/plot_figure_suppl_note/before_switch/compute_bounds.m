function compute_bounds(pop)
%% This code retrieves the 95% confidence bounds on the model. Pop = 1 or 2

%%
best = dlmread('best.txt');

chibest = best(end);

measured_best = dlmread('smim1_mod.txt');


for j = 1:18
    
    xbest = measured_best(j, pop + 1);
    
    
    fName = ['res_', num2str(pop),'_', num2str(j)];
    
    
    if exist(fName,'file') == 2
        
        
        data = dlmread(fName);
        
        
        
        data = sortrows(data, 1);
        
        x = data(:,1);
        chi = data(:,end);
        
        x(chi > chibest + 3.8) = [];
        
        
        chi(chi > chibest + 3.8) = [];
        
        
        
        [massi,imax] = max(x);
        [mini,imini] = min(x);
        
        dlmwrite(['bounds', num2str(pop), '.txt'], [xbest massi mini chi(imax) chi(imini)], 'delimiter', '\t', '-append');
        
    end
    
    
end

end
