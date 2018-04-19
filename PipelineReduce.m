% TEST 'MRA' ALGORITHM
%
% 1) Generate random brain scouts
% 2) Simulate AR signals
% 3) Simulate recordings from scouts (noise)
% 4) Compute sources
% 5) FW model scouts 
%    - ALL
%    - ONE BY ONE
% 6) Reduce Dimensionality
% 7) ROC curve
% 
% Author: Guimar Niso, 2014
% Author: Cristina Gil Avila, 2018



clc; clear; close all

SubjectNames = {'Subject02'};

% Start a new report
bst_report('Start', []);

SIM = 1:2:50; % Simulations

nscouts = 6;
samples = 1000;
srate = 1000; 
noise1 = 0.3;
noise2 = 0.2;
niter = 10;
ninic = 50;

condition = 'Micha';
scoutsRAND = {'RAND', {'01', '02', '03', '04', '05', '06'}};
scoutsRED = {'Destrieux', {'G_Ins_lg_and_S_cent_ins L', 'G_Ins_lg_and_S_cent_ins R', 'G_and_S_cingul-Ant L', 'G_and_S_cingul-Ant R', 'G_and_S_cingul-Mid-Ant L', 'G_and_S_cingul-Mid-Ant R', 'G_and_S_cingul-Mid-Post L', 'G_and_S_cingul-Mid-Post R', 'G_and_S_frontomargin L', 'G_and_S_frontomargin R', 'G_and_S_occipital_inf L', 'G_and_S_occipital_inf R', 'G_and_S_paracentral L', 'G_and_S_paracentral R', 'G_and_S_subcentral L', 'G_and_S_subcentral R', 'G_and_S_transv_frontopol L', 'G_and_S_transv_frontopol R', 'G_cingul-Post-dorsal L', 'G_cingul-Post-dorsal R', 'G_cingul-Post-ventral L', 'G_cingul-Post-ventral R', 'G_cuneus L', 'G_cuneus R', 'G_front_inf-Opercular L', 'G_front_inf-Opercular R', 'G_front_inf-Orbital L', 'G_front_inf-Orbital R', 'G_front_inf-Triangul L', 'G_front_inf-Triangul R', 'G_front_middle L', 'G_front_middle R', 'G_front_sup L', 'G_front_sup R', 'G_insular_short L', 'G_insular_short R', 'G_oc-temp_lat-fusifor L', 'G_oc-temp_lat-fusifor R', 'G_oc-temp_med-Lingual L', 'G_oc-temp_med-Lingual R', 'G_oc-temp_med-Parahip L', 'G_oc-temp_med-Parahip R', 'G_occipital_middle L', 'G_occipital_middle R', 'G_occipital_sup L', 'G_occipital_sup R', 'G_orbital L', 'G_orbital R', 'G_pariet_inf-Angular L', 'G_pariet_inf-Angular R', 'G_pariet_inf-Supramar L', 'G_pariet_inf-Supramar R', 'G_parietal_sup L', 'G_parietal_sup R', 'G_postcentral L', 'G_postcentral R', 'G_precentral L', 'G_precentral R', 'G_precuneus L', 'G_precuneus R', 'G_rectus L', 'G_rectus R', 'G_subcallosal L', 'G_subcallosal R', 'G_temp_sup-G_T_transv L', 'G_temp_sup-G_T_transv R', 'G_temp_sup-Lateral L', 'G_temp_sup-Lateral R', 'G_temp_sup-Plan_polar L', 'G_temp_sup-Plan_polar R', 'G_temp_sup-Plan_tempo L', 'G_temp_sup-Plan_tempo R', 'G_temporal_inf L', 'G_temporal_inf R', 'G_temporal_middle L', 'G_temporal_middle R', 'Lat_Fis-ant-Horizont L', 'Lat_Fis-ant-Horizont R', 'Lat_Fis-ant-Vertical L', 'Lat_Fis-ant-Vertical R', 'Lat_Fis-post L', 'Lat_Fis-post R', 'Pole_occipital L', 'Pole_occipital R', 'Pole_temporal L', 'Pole_temporal R', 'S_calcarine L', 'S_calcarine R', 'S_central L', 'S_central R', 'S_cingul-Marginalis L', 'S_cingul-Marginalis R', 'S_circular_insula_ant L', 'S_circular_insula_ant R', 'S_circular_insula_inf L', 'S_circular_insula_inf R', 'S_circular_insula_sup L', 'S_circular_insula_sup R', 'S_collat_transv_ant L', 'S_collat_transv_ant R', 'S_collat_transv_post L', 'S_collat_transv_post R', 'S_front_inf L', 'S_front_inf R', 'S_front_middle L', 'S_front_middle R', 'S_front_sup L', 'S_front_sup R', 'S_interm_prim-Jensen L', 'S_interm_prim-Jensen R', 'S_intrapariet_and_P_trans L', 'S_intrapariet_and_P_trans R', 'S_oc-temp_lat L', 'S_oc-temp_lat R', 'S_oc-temp_med_and_Lingual L', 'S_oc-temp_med_and_Lingual R', 'S_oc_middle_and_Lunatus L', 'S_oc_middle_and_Lunatus R', 'S_oc_sup_and_transversal L', 'S_oc_sup_and_transversal R', 'S_occipital_ant L', 'S_occipital_ant R', 'S_orbital-H_Shaped L', 'S_orbital-H_Shaped R', 'S_orbital_lateral L', 'S_orbital_lateral R', 'S_orbital_med-olfact L', 'S_orbital_med-olfact R', 'S_parieto_occipital L', 'S_parieto_occipital R', 'S_pericallosal L', 'S_pericallosal R', 'S_postcentral L', 'S_postcentral R', 'S_precentral-inf-part L', 'S_precentral-inf-part R', 'S_precentral-sup-part L', 'S_precentral-sup-part R', 'S_suborbital L', 'S_suborbital R', 'S_subparietal L', 'S_subparietal R', 'S_temporal_inf L', 'S_temporal_inf R', 'S_temporal_sup L', 'S_temporal_sup R', 'S_temporal_transverse L', 'S_temporal_transverse R'}};

a = 1;

for i = SIM
        
disp('*** Process: Generate random brain scouts')
% Process: Generate random brain scouts ===================================
bst_process(...
    'CallProcess', 'process_random_scouts_generator', [], [], ...
    'subjectname', SubjectNames{1}, ...
    'condition', condition, ...
    'nscouts', nscouts, ...
    'nvertices', [10, 100]);

disp('*** Process: Simulate AR signals')
% Process: Simulate AR signals ============================================
sFilesAR = bst_process(...
    'CallProcess', 'process_simulate_ar', [], [], ...
    'subjectname', SubjectNames{1}, ...
    'condition', condition, ...
    'samples', samples, ...
    'srate', srate, ...
    'A', ['Nnodes=6;' 10 'A1 = zeros (Nnodes,Nnodes);' 10 'A2 = zeros (Nnodes,Nnodes);' 10 'A3 = zeros (Nnodes,Nnodes);' 10 'A4 = zeros (Nnodes,Nnodes);' 10 '% t=1' 10 'A1(1,1) = 1.3393;' 10 'A1(4,4) = 0.25*sqrt(2);' 10 'A1(4,5) = 0.25*sqrt(2);' 10 'A1(5,4) = -0.25*sqrt(2);' 10 'A1(5,5) = 0.25*sqrt(2);' 10 '% t=2' 10 'A2(1,1) = -0.5823;' 10 'A2(2,1) = 0.5;' 10 'A2(4,1) = -0.5;' 10 '% t=3' 10 'A3(3,1) = 0.4;' 10 'A3(6,6) = -0.25*sqrt(2);' 10 '% t=4' 10 'A4(6,6) = 0.25*sqrt(2);' 10 '' 10 '% Matrix of coefficients' 10 'A = [A1 A2 A3 A4];'], ...
    'b', 'b =[0.1 0.11 0.12 0.13 0.14 0.15]; ', ...
    'C', ['C = eye(6,6);' 10 '']);


disp('*** Process: Simulate recordings from scouts (noise)')
% Process: Simulate recordings from scouts (noise) ========================
sFilesRec = bst_process('CallProcess', 'process_simulate_recordings_noise', sFilesAR, [], ...
    'scouts',      scoutsRAND, ...
    'savesources', 0, ...
    'noise1',      noise1, ...
    'noise2',      noise2);


disp('*** Process: Compute sources')
% Process: Compute sources ================================================
sFilesSou = bst_process(...
    'CallProcess', 'process_inverse', sFilesRec, [], ...
    'comment', '', ...
    'method', 1, ...  % Minimum norm estimates (wMNE)
    'wmne', struct(...
         'NoiseCov', [], ...
         'InverseMethod', 'wmne', ...
         'ChannelTypes', {{}}, ...
         'SNR', 3, ...
         'diagnoise', 0, ...
         'SourceOrient', {{'fixed'}}, ...
         'loose', 0.2, ...
         'depth', 1, ...
         'weightexp', 0.5, ...
         'weightlimit', 10, ...
         'regnoise', 1, ...
         'magreg', 0.1, ...
         'gradreg', 0.1, ...
         'eegreg', 0.1, ...
         'ecogreg', 0.1, ...
         'seegreg', 0.1, ...
         'fMRI', [], ...
         'fMRIthresh', [], ...
         'fMRIoff', 0.1, ...
         'pca', 1), ...
    'sensortypes', 'MEG', ...
    'output', 1);  % Kernel only: shared


% Process: FW model scouts ALL=============================================
bst_process(...
    'CallProcess', 'process_fw_scouts', sFilesSou, [], ...
    'atlas', 'RAND', ...
    'all', 1);

% Process: FW model scouts ONE BY ONE =====================================
bst_process(...
    'CallProcess', 'process_fw_scouts', sFilesSou, [], ...
    'atlas', 'RAND', ...
    'all', 0);


% for i = SIM
% niter = i;

disp('*** Process: Reduce Dimensionality')
% Process: Reduce Dimensionality ==========================================
sFiles = bst_process('CallProcess', 'process_reduce_dimension', sFilesSou, [], ...
    'scouts',    scoutsRED, ...
    'scoutfunc', 1, ...  % Mean
    'ninic',     ninic, ...
    'niter',     niter, ...
    'isnorm',    0);


disp('*** Process: ROC curve')
% Process: ROC curve ======================================================
bst_process(...
    'CallProcess', 'process_roc_curve', sFilesSou, [], ...
    'subjectname', SubjectNames{1}, ...
    'condition', condition, ...
    'GT', 'RAND', ...
    'RES1', ['MRA:',num2str(niter),'iter']);

% Save ROC results
roc_result = importdata('/home/bic/guiomar/.brainstorm/process/roc_result_mra.mat');
ROC_orig{a} = roc_result;

% Process: ROC curve ======================================================
bst_process(...
    'CallProcess', 'process_roc_curve', sFilesSou, [], ...
    'subjectname', SubjectNames{1}, ...
    'condition', condition, ...
    'GT', 'FWM', ...
    'RES1', ['MRA:',num2str(niter),'iter']);

% Save ROC results
roc_result = importdata('/home/bic/guiomar/.brainstorm/process/roc_result_mra.mat');
ROC_fw1{a} = roc_result;

% Process: ROC curve ======================================================
bst_process(...
    'CallProcess', 'process_roc_curve', sFilesSou, [], ...
    'subjectname', SubjectNames{1}, ...
    'condition', condition, ...
    'GT', 'FWM_02', ...
    'RES1', ['MRA:',num2str(niter),'iter']);

% Save ROC results
roc_result = importdata('/home/bic/guiomar/.brainstorm/process/roc_result_mra.mat');
ROC_fwall{a} = roc_result;

% %%%% Remove this remaining thing
% sStudy = bst_get('Study', sFilesRec.iStudy);



% DELETE ALL FILES AND START AGAIN =====================================

% ==== Delete Atlas 

% Get subject
sSubject = bst_get('Subject', SubjectNames{1});
% Get surface file
CortexFile = sSubject.Surface(sSubject.iCortex).FileName;

% Load surface
sSurface = in_tess_bst(CortexFile);

% Find the index of the altas
iAtlas4 = find(strcmpi({sSurface.Atlas.Name}, 'RAND'));
iAtlas3 = find(strcmpi({sSurface.Atlas.Name}, 'FWM'));
iAtlas2 = find(strcmpi({sSurface.Atlas.Name}, 'FWM_02'));
iAtlas1 = find(strcmpi({sSurface.Atlas.Name}, ['MRA:',num2str(niter),'iter']));

% Remove the atlas  %% becareful order!!
sSurface.Atlas(iAtlas1) = [];
sSurface.Atlas(iAtlas2) = [];
sSurface.Atlas(iAtlas3) = [];
sSurface.Atlas(iAtlas4) = [];

% Save the modifications back to the surface file
bst_save(file_fullpath(CortexFile), sSurface, 'v6');

close all;
a = a+1;

% ==== Delete other files

disp('*** Process: Delete data files')
% Process: Delete data files
bst_process(...
    'CallProcess', 'process_delete', ...
    [sFilesAR, sFilesRec, sFilesSou], [], ...
    'target', 1);  % Delete data fles


end

% % sSurface.Atlas(iAtlas2) = [];
% % sSurface.Atlas(iAtlas3) = [];
% % sSurface.Atlas(iAtlas4) = [];




% =========================================================================

% % % Save and display report
% % ReportFile = bst_report('Save', sFiles); %%% sFiles??
% % bst_report('Open', ReportFile);



%% PLOT ROC CURVE ==========================================================

COLOR = jet(length(SIM));
markerSize = 10;
lineWidth = 1;

figure(1)
hold on;
for i=1:length(SIM)
    TPR1 = ROC_orig{i}(1,:);
    FPR1 = ROC_orig{i}(2,:);
    
    plot( [ 0 FPR1 ], [ 0 TPR1 ], '.-','color',COLOR(i,:),'MarkerSize',markerSize,'LineWidth', lineWidth);% 'MarkerFaceColor',
end
plot( [ 0 1 ], [ 0 1 ], 'k:');
xlabel('FPR');
ylabel('TPR');
axis square;
title('MRA vs RAND')
legend([strsplit(num2str(SIM),' '),'rand'])


figure(2)
hold on;
for i=1:length(SIM)
    TPR1 = ROC_fwall{i}(1,:);
    FPR1 = ROC_fwall{i}(2,:);
    
    plot( [ 0 FPR1 ], [ 0 TPR1 ], '.-','color',COLOR(i,:),'MarkerSize',markerSize,'LineWidth', lineWidth);% 'MarkerFaceColor',
end
plot( [ 0 1 ], [ 0 1 ], 'k:');
xlabel('FPR');
ylabel('TPR');
axis square;
title('MRA vs FWall')
legend([strsplit(num2str(SIM),' '),'rand'])


figure(3)
hold on;
for i=1:length(SIM)
    TPR1 = ROC_fw1{i}(1,:);
    FPR1 = ROC_fw1{i}(2,:);
    
    plot( [ 0 FPR1 ], [ 0 TPR1 ], '.-','color',COLOR(i,:),'MarkerSize',markerSize,'LineWidth', lineWidth);% 'MarkerFaceColor',
end
plot( [ 0 1 ], [ 0 1 ], 'k:');
xlabel('FPR');
ylabel('TPR');
axis square;
title('MRA vs FW1')
legend([strsplit(num2str(SIM),' '),'rand'])

