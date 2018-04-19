% b: In Brainstorm export to Matlab (cortex of ~15000 vertices)

Atlas = {b.Atlas.Name}';


Nthres = 0.7:0.1:0.8;

Resultados.Nscouts =[];
Resultados.Ntres =[];


for ik = 1:numel(Nthres)
    
    miatlas = ['MLF:',num2str(Nthres(ik)),'thres'];
    Nik = strcmp(Atlas,miatlas);
    
    Resultados(ik).Nscouts = numel(b.Atlas(Nik).Scouts);
    Resultados(ik).Nthres = Nthres(ik);
    
end

% PLOT FIGURE

figure ()

plot([Resultados.Nthres], [Resultados.Nscouts],'or-');
axis ([0 1 1 15000]);

ylabel;
xlabel,



%%%

% For PipelineReduce.m

disp('*** Process: Cluster IK')
% Process: Cluster IK ==========================================

Thresholds = 0.1:0.05:1;

for iT = 1:numel(Thresholds)
    
    THRES = Thresholds(iT);
    
    % tic toc, save time!!
    
    % Process: Cluster IK
    bst_process('CallProcess', 'process_cluster_ik', ...
        sFilesSou, [], ...
        'scouts',    scoutsRED, ...
        'scoutfunc', 1, ...  % Mean
        'threshold', THRES, ...
        'isnorm',    0);
    
end