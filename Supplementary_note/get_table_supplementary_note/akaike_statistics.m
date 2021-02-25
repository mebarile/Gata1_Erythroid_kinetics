function akaike_statistics
%% This code ranks the models according to theis corrected Akaike index

%%
data = fileread('list_murk_genes.txt');
x = strsplit(data);

mat = zeros(8,1);

for h = 1:89
    
    gene = x{h};
    
    
    
    
    aic = dlmread(['./akaike_out/',gene,'_aic.txt']);
    vv = 1:8;
    aic = [vv',aic];
    
    aic_sort = sortrows(aic,6);
    
    
    for i = 1:8
        mat(i) = mat(i) + find(aic_sort(:,1)==i);
    end
    
end

model_specifics = {'all parameters constant',...
    'only degradation changes',...
    'only splicing changes',...
    'splicing and degradation change',...
    'only transcription changes',...
    'transcription and degradation change',...
    'transcription and splicing change',...
    'all parameters change',...
    };

[~,index] = sort(mat,'ascend');

'Ranking of models:'
model_specifics{index}
end
