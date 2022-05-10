function [casoX,tot_gen,tot_branch,tot_con,contingency] = loadSCOPF_data(args)
% Esta función carga la información del sistema a correr
% se usa para organizar el código
% DirName especifica la ruta donde se encuentra el archivo .mat con los datos
if nargin == 0
    % caso 14 nodos
    casoX = loadcase('case14');
    casoX.isshunt = false; %no tener en cuenta compensación
    casoX.gen(:,21) = [12 19 24 21 24]'; % agrega factores participacion
    % modificar limites para rutina
    casoX.branch(:,6) = 30;
    % parte de la estructura de contingencias
    nodosgen = [1;2;8]; % nodos a los cuales están conectados los generadores que participan
    tot_gen = size(casoX.gen,1);
    % tot_gen = 1;
    tot_branch = size(casoX.branch,1);
    % tot_branch = 0;
    tot_con = tot_gen+tot_branch;
    % cell_out =mat2cell([1],ones(tot_con,1),1);
    cell_out = mat2cell([(1:size(casoX.gen,1))';(1:size(casoX.branch,1))'],ones(tot_con,1),1);
    cell_Gk = mat2cell(repmat(nodosgen,tot_con,1),length(nodosgen)*ones(tot_con,1),1);
    contingency = struct('out',cell_out,'G_k',cell_Gk);
else
    % cargar los de la competencia de arpa
    DirName = 'C:\Users\Camilo Acosta\Documents\GitHub\SCOPF_benders\Data_sets\';
%     DirName = 'C:\Users\pirit\Dropbox\FINAL ARPA\FINAL_FINAL-master\';
%     DirName = 'C:\Users\Camilo Acosta\Dropbox\FINAL ARPA\FINAL_FINAL-master\';
    % DirName = '/home/camilo/Dropbox/FINAL ARPA/FINAL_FINAL-master/'; % linea para linux
%     DirName = 'C:\Users\Usuario UTP\Documents\Data\';
    
    FileName = 'network02O_23-scenario1.mat'; % red de 500 nodos
%     FileName = 'network_75O-040_scenario17.mat'; % red de 2712 nodos
%     FileName = 'network_12O-050_scenario39.mat'; % red de 9591 nodos
    MatData = [DirName FileName];
    CaseData = fullfile(DirName,FileName);
    casoX = loadcase(CaseData);
    load(MatData); % carga sin necesidad
    %contingency = contingency([1:21 155:233]); % sistema 500 nodos
    % contingency = contingency([1:21]);
    % contingency = contingency([156]);
    % contingency = contingency([155:164]);
    %tot_gen = 21; % sistema 500 nodos
    %tot_branch = 79; % sistema 500 nodos
    %tot_con = tot_gen + tot_branch;
    % Sistema 2742 nodos
    % contingency = contingency([1:16 181:364]); % sistema 2742
    % contingency = contingency([1:16]);
    % tot_gen = 16; % 16 para sistema 2742 nodos
    % tot_branch = 184; % 184 nodos para sistema 2742 nodos
    % tot_con = tot_gen + tot_branch;
    % sistema de 9591 nodos
%     contingency = contingency([1:250 2000:2250]);
%     tot_gen = 250;
%     tot_branch = 250;
%     tot_con = tot_branch + tot_gen;

end

%% ****************************************************************************
%                               SOFTLIMS SETTINGS
casoX.softlims.VMIN.hl_mod = 'none'; % Not vmin soft limits
casoX.softlims.VMAX.hl_mod = 'none'; % Not vmax soft limits
casoX.softlims.PMIN.hl_mod = 'none'; % Not pmin soft limits
casoX.softlims.PMAX.hl_mod = 'none'; % Not pmax soft limits
casoX.softlims.QMIN.hl_mod = 'none'; % Not qmin soft limits
casoX.softlims.QMAX.hl_mod = 'none'; % Not qmax soft limits
casoX.softlims.ANGMIN.hl_mod = 'none'; % Not angmin soft limits
casoX.softlims.ANGMAX.hl_mod = 'none'; % Not angmax soft limits
% Only branch limits are important
casoX.softlims.RATE_A.hl_mod = 'remove'; % change limit to infinity
casoX.softlims.RATE_A.idx = find(casoX.branch(:,6)>0); % index of lines to be replaced limit.
casoX.softlims.RATE_A.cost = 500e3*ones(size(casoX.softlims.RATE_A.idx,1),1); % violation cost (this cost function is changed internally)

