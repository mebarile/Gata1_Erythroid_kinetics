function launch_main_best_fit
% This codes launches the code to fit the data with the basic model for each MURK gene

%%
data = fileread('list_murk_genes.txt');
x = strsplit(data);

t0 = 0;

for i = 1:89
    
    gene = x{i};
    
    
    main_find_best(gene, t0);
        
    
end

end
