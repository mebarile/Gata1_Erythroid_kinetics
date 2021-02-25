function launch_main_all_comb
% This codes launches the code to fit the data with the all models model for each MURK gene

%%
data = fileread('list_murk_genes.txt');
x = strsplit(data);

t0 = 0;



for h = 1:89
    
    gene = x{h};
    
    
    
    
    for i = 0:1
        for j = 0:1
            for k = 0:1
                
                main_all_comb(gene, t0, [i,j,k]);
                
                
            end
        end
    end
    
    
    
end


end
