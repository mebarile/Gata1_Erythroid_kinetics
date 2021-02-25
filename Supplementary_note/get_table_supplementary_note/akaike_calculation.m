function akaike_calculation
% This code computes the corrected Akaike index for each model for each MURK gene

%%
data = fileread('list_murk_genes.txt');
x = strsplit(data);

n = 20; % number data points


for h = 1:89
    
    gene = x{h};
    
    aic = zeroes(8,5);
    
    p = 1;
    

    for i = 0:1
        for j = 0:1
            for m = 0:1
                
                keep = [i, j, m];
                
                nome_best = ['./combinations/best_',gene,'_',num2str(keep(1)), num2str(keep(2)), num2str(keep(3)), '.txt'];
                
                data = dlmread(nome_best);
                
                k = sum(keep) + 6;
                
                if k == 6
                    
                    k = 5;
                    
                end
                
                chi = data(end);
                
                
                aic(p,:) = [keep, chi, akaike(chi, n, k)];
                
                p = p + 1;
                
                
            end
        end
    end
    
    dlmwrite(['akaike_out/', gene, '_aic.txt'], aic)
    
    
end

end
