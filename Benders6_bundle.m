clear all;clc; close all
% Archivo principal rutina de benders para el SCOPF
tic
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
 MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
 QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
%% Load Case
[casoX,tot_gen,tot_branch,tot_con,contingency] = loadSCOPF_data(1);

%% **********************************************************************************************************
% 											MPOPT OPTIONS
% **********************************************************************************************************
mpopt = mpoption('opf.ac.solver','IPOPT','verbose',0,'out.all',0);
% modelo AC
mpopt.model = 'AC';
optionsipopt = ipopt_options; % cargar configuracion default
optionsipopt.print_level = 3;
optionsipopt.linear_solver = 'ma57';
optionsipopt.fixed_variable_treatment = 'make_constraint';
% change ipopt parameters only for networks with 1000 buses or larger
if size(casoX.bus,1)>=100
	optionsipopt.mu_strategy = 'adaptive';
	% optionsipopt.alpha_for_y = 'min-dual-infeas';
	optionsipopt.corrector_type = 'affine';
	optionsipopt.linear_system_scaling = 'slack-based';
	optionsipopt.bound_relax_factor = 1e-8;
	optionsipopt.print_level = 3;
    optionsipopt.dual_inf_tol = 1;
    optionsipopt.compl_inf_tol = 1e-4;
    optionsipopt.ma57_automatic_scaling = 'yes';
    optionsipopt.mehrotra_algorithm = 'no';
    % optionsipopt.mu_strategy = 'adaptive';
    optionsipopt.mu_oracle = 'quality-function';
    optionsipopt.fixed_mu_oracle = 'quality-function';
end
optionsipopt.max_iter = 3000;
optionsipopt.max_cpu_time = 3600;
optionsipopt.fixed_variable_treatment = 'relax_bounds';
% optionsipopt.acceptable_tol = 1;
optionsipopt.print_level = 3;
optionsipopt.acceptable_constr_viol_tol = 1e-4;
mpopt = mpoption(mpopt,'ipopt.opts',optionsipopt,...
    'opf.start',0);
% ---------------------------------------------------------
%% ********************************************************
%                       CARGAR TOGGLES
% *********************************************************
% shunt devices
casoX = toggle_shuntgen(casoX);
% nodal mismatch power equations
casoX = toggle_softlims_ext1(casoX,'on'); % bus balance equations with slack constraints
% Fixed feasible variables and obtain initial cost
now1 = tic();
results0 = runopf(casoX,mpopt); % base case OPF
timeOPFbase = toc(now1);
fprintf('el OPF se demora: %d s',timeOPFbase);
% *** bundle ****
variables = results0.var;
% *** fin bundle ***
basecost = results0.f; %
basecos_inicial = basecost;
fprintf('El primer OPF es %i\n',results0.success)
% store solution w/o cuts
Pgen_base = results0.gen(:,2); % se activa solo para resultados parciales
% Qgen_base = results0.gen(:,3);
% Vbus_base = results0.bus(:,8);
ZLB = []; % stores Z_lo
ZUB=[];  % stores Z_ub
ZUB1=[]; % almacena el Z_ub sin cambios
% flow_P = results0.branch(:,14);


% mpopt = mpoption(mpopt,'opf.start',2); % warm start (is not neccesary but accelerates OPF convergence)
z_cut_obj = 0;
adaptive_cut = 'normal'; % Adaptive cut is to probing the actual cut
%% ******************** Contingency response ************
%  Create structure (cell) to store cut information
cut_x_contingency.f_obj =[]; % is a vector with variable size or fixed in cuts numbers
cut_x_contingency.mu_Pg =[]; % to dual multipliers for Pg
cut_x_contingency.constant =[]; % it is sum for constant in cut equation
cut_x_contingency.mu_volt = []; % dual multipliers for Vg
cut_x_contingency = repmat(cut_x_contingency,tot_con,1); % for every contingency there are one structure
% structure to store current contingency solution
contingency_sol.cost = 0; % cost of contingency
contingency_sol.delta = 0; % deltak value in contingency
contingency_sol.gen = zeros(size(results0.gen,1),2); % stores Pg, Qg
contingency_sol.Vbus = zeros(size(results0.bus,1),1); % stores Vb
contingency_sol.s_branch = zeros(size(results0.branch,1),1); % stores overload of branches
contingency_sol.s_balance = zeros(size(results0.bus,1),4); % stores balance slack variables for buses in P+,P-,Q+,Q-
contingency_sol = repmat(contingency_sol,tot_con,1); % replicate for every contingency
% **********************************************************
% 							BEGIN BENDERS PROCEDURE
% **********************************************************
% obj_func_prev = Inf(tot_con,1); % begin with inf value
% optionsipopt.alpha_for_y = 'safer-min-dual-infeas';
% mpopt = mpoption(mpopt,'ipopt.opts',optionsipopt);
mpopt1 = mpopt; % tomar los ajustes actuales
optionsipopt1 = optionsipopt;
% optionsipopt1.print_level = 5;
% optionsipopt1.mu_strategy = 'adaptive';
% optionsipopt1.mu_oracle = 'quality-function';
% optionsipopt1.fixed_mu_oracle ='average_compl';
% optionsipopt1.adaptive_mu_globalization = 'kkt-error';
% optionsipopt1.mu_linear_decrease_factor = 0.8;
% optionsipopt1.alpha_for_y = 'min-dual-infeas';
% optionsipopt1.accept_every_trial_step = 'yes';
% optionsipopt1.expect_infeasible_problem = 'yes';
mpopt1 = mpoption(mpopt1,'ipopt.opts',optionsipopt1);
upbound = Inf;
lobound = -Inf;
it_cuts = 1; % current iteration of benders cut procedure
cost_accept = 0.1*basecost;
no_mejora = 0; % if not improvement in benders procedure
cuts_tot = 30;
% parámetros bundle
delta_bundle = inf;
tol_bundle = 1e-5; % de prueba
m_bundle = 0.01;
t_bundle = 1e3; % valor de inicio del multiplicador (beta en doc tesis)
% fin parámetros bundle
while (it_cuts <=cuts_tot)&&(no_mejora <= 10)&& ((upbound-lobound) >= 	0.1*lobound) && (delta_bundle > tol_bundle) % outer loop, controls the process entirely
	disp(strcat('iteracion de benders No ',num2str(it_cuts)));
	% Foreach contingency
    Z_ub = basecost; % c(y_fixed) + contingency_cost (upper bound)
	muPgx = 0;
	muVbx = 0;
	parfor con=1:tot_con % maybe a parfor loop
		% additional options to IPOPT
		casoX_ext = results0; % copy results
		casoX_ext.Pgk.Pg0 = casoX_ext.gen(:,PG)/casoX_ext.baseMVA; % Active power base case (p.u.)
        casoX_ext.Vdev.V_i = casoX_ext.bus(:,8); % Voltage magnitude base case
		casoX_ext.Pgk.G_k = contingency(con).G_k; % generators in area
		% casoX_ext = rmfield(casoX_ext,{'var','nle','nli','qdc','f'}); % esta linea está rara que no funcione
		if con<=tot_gen
			casoX_ext.gen(contingency(con).out,GEN_STATUS) = 0; % simulates generator contingency
			casoX_ext.Pg_out = casoX_ext.gen(contingency(con).out,2)/casoX_ext.baseMVA; % Active power outage
		else
			casoX_ext.branch(contingency(con).out,11) = 0; % simulate branch out contingency
		end
		casoX_ext.fobj = casoX_ext.f;
        casoX_ext = rmfield(casoX_ext,{'nle','nli','qdc','f'});

		% ************************* Call toggles in contingency ***********************************
		casoX_ext.taoV = 1e-3;
		% casoX_ext = toggle_Pgk_polyhedral(casoX_ext);
		% casoX_ext = toggle_Pgk_newform(casoX_ext); % add active power redispatch equations
		casoX_ext = toggle_Pgk_newform2(casoX_ext); % add active power redispatch equations
		casoX_ext = toggle_voltages_newformulation(casoX_ext); % add voltages behavior equations in contingency for every generator or bus with generators.
        % fprintf('======================= contingencia %i ==================\n',con)
		% results2 = solvempec(casoX_ext,mpopt1);
		% results2 = solvempec_test(casoX_ext,mpopt1);
		% try
		% 	results2 = solvempec_test1(casoX_ext,mpopt1);
		% catch
		% 	error('problemas en la contingencia %i e iteracion %i',con,it_cuts);
		% end
        now2 = tic();
		results2 = solvempec_test1(casoX_ext,mpopt1);
        time2 = toc(now2);
        fprintf('\n el tiempo es %d s\n',time2);
        % results2 = solvempec1(casoX_ext,mpopt1);
		% results2 = runopf(casoX_ext,mpopt1);
		% if results2.success==0
			% fprintf('El OPF de la contingencia %i es %i\n',con,results2.success)
		% end

        % *************************   End toggles *****************************************

        % ************************* Creacion de cortes de Benders

        % [ConstantPg,muPg] = benderscutPg(results2);
		[ConstantPg,muPg] = benderscutPg_test(results2);
		[ConstantVb,muVb] = benderscutVb(results2);
		muPgx = muPgx + muPg;
		muVbx = muVbx + muVb;
		% [ConstantVb,muVb] = benderscutVb1(results2);
        % storeVB = [storeVB ConstantVb];
        % storemuVb = [storemuVb muVb];


        % introduccion de un cambio necesario
        % ConstantVb = 0;
        % muVb = zeros(size(results2.bus(:,1)));

		% ************************** Calculo de F(x) valor de la contingencia **************
		if isfield(results2.qdc,'Pgkcost')
            CostoPgk = results2.qdc.Pgkcost;
        else
            CostoPgk = 0;
        end

        if isfield(results2.qdc,'Vgkcost')
            CostoVgk = results2.qdc.Vgkcost;
        else
            CostoVgk = 0;
        end

        Obj_func = results2.f - CostoPgk - CostoVgk;
		% ********************************************************************************

		% *********************** asignacion de valores a la estructura de contingencia (cortes)  ***********************
        tao = 1.0;
		cut_x_contingency(con).f_obj = [cut_x_contingency(con).f_obj;tao*Obj_func]; % almacena funcion objetivo de cada contingencia
		cut_x_contingency(con).constant = [cut_x_contingency(con).constant;...
			tao*(ConstantPg+ConstantVb)]; % almacena valor escalar de las restricciones por contingencia
		cut_x_contingency(con).mu_Pg = [cut_x_contingency(con).mu_Pg,muPg]; % Pg multiplier
		cut_x_contingency(con).mu_volt = [cut_x_contingency(con).mu_volt,muVb]; % Vb multiplier
        Z_ub = Z_ub + Obj_func; % adding objective function of each contingency
        % ********************* Store contingency values ***************************
        contingency_sol(con).cost = Obj_func;
        contingency_sol(con).delta = results2.delta*results2.baseMVA;
        contingency_sol(con).gen = results2.gen(:,[PG;QG]);
        contingency_sol(con).Vbus = results2.bus(:,8);
		contingency_sol(con).s_branch = results2.softlims.RATE_A.overload;
		varnames = {'s_p_unbal_plus','s_p_unbal_minus',...
		's_q_unbal_plus','s_q_unbal_minus'};
		for i=1:numel(varnames)
			varname = varnames{i};
			contingency_sol(con).s_balance(:,i) = results2.var.val.(varname)*results2.baseMVA;
		end


		%adaptive_cut = 'normal';

		% else
			% % reduce original cost, attempting cut off higher costs
			% % calculate tao from overrun
			% mu_1 = zeros(length(results2.Pgk.Gin),1);
			% % represents the condition Pgmax <= Pg + alpha_g*deltak
			% mu_2 = zeros(length(results2.Pgk.Gin),1);
			% % represents the condition Pgmin >= Pg + alpha_g*deltak
			% mu_3 = zeros(length(results2.Pgk.Gin),1);
			% % represents the condition Pgk = Pg for generators out of area
			% mu_4 = results2.var.mu.u.Pg - results2.var.mu.l.Pg;
			% % Fill mu for generators in area
			% mu_1(results2.Pgk.Gin) = results2.lin.mu.u.res2_3 - results2.lin.mu.l.res2_3;
			% mu_2(results2.Pgk.Gin) = results2.lin.mu.u.res3_2;
			% mu_3(results2.Pgk.Gin) = results2.lin.mu.l.res4_2;
			% mu_4(results2.Pgk.Gin) = 0;
			% % shift to internal. Only for matching orders
			% results2.order.state = 'i';
			% aux1 = zeros(size(results2.gen,1),1);
			% % mapping to external. Both off line generators and contingency out generator must have mu = 0
			% mu_1 = i2e_data(results2,mu_1,aux1,'gen');
			% mu_2 = i2e_data(results2,mu_2,aux1,'gen');
			% mu_3 = i2e_data(results2,mu_3,aux1,'gen');
			% mu_4 = i2e_data(results2,mu_4,aux1,'gen');
			% results2.order.state = 'e';
			% % calculate cut as taylor first order approximation
			% % calculate mu*f(Pg)
			% mult_vector = mu_1.*results2.ygk(:,1)+mu_2.*results2.ygk(:,2)-mu_3.*results2.ygk(:,3) + mu_4;
			% const_term = mult_vector'*results0.gen(:,PG)/results2.baseMVA;
			% % assign values to cut struct of the k-th contingency
			% cut_x_contingency(con).f_obj = [cut_x_contingency(con).f_obj;Obj_func];
			% cut_x_contingency(con).mu_b = [cut_x_contingency(con).mu_b;const_term];
			% cut_x_contingency(con).mu_var = [cut_x_contingency(con).mu_var,mult_vector];
			% Z_ub = Z_ub + Obj_func;

			% % tao = 0.2;
			% % cut_x_contingency(con).f_obj(end)= tao*cut_x_contingency(con).f_obj(end);
			% % cut_x_contingency(con).mu_b(end)= tao*cut_x_contingency(con).mu_b(end);
			% % cut_x_contingency(con).mu_var(:,end) = tao*cut_x_contingency(con).mu_var(:,end);
			% % Z_ub = Z_ub + Obj_func;
			% % adaptive_cut = 'adaptive'; % selects how to change current cut
		% end % ends adaptive cut
        % casoX_ext = rmfield(casoX_ext,'userfcn');
	end % ends contingencies
	% Z_ub = Z_ub + Z_lo - z_cut_obj; % calculates Z upper bound solution
    auxP = zeros(size(results0.gen,1),1);
    auxV = zeros(size(results0.bus,1),1);
    results0.order.state = 'i'; 
    variables.Pgen = i2e_data(results0,variables.val.Pg,auxP,'gen',1);
    variables.Vm = i2e_data(results0,variables.val.Vm,auxV,'bus',1);
    results0.order.state = 'e';
	if it_cuts==1
		%revisar que hacer antes de pasar al otro
		f_xk = Z_ub; % como por llevar las cuentas
		best_var = variables;
		bestmuPg = muPgx;
		bestmuVb = muVbx;
	elseif (f_xk - Z_ub >= m_bundle*delta_bundle)
		epsilon_vector = [muPgx;muVbx]-[bestmuPg;bestmuVb];
		var_vector = [variables.Pgen;variables.Vm]-[best_var.Pgen;best_var.Vm];
		norma_eps = norm(epsilon_vector,2); % calculo la norma del vector
		% t_bundle = (1/t_bundle+var_vector'*epsilon_vector/norma_eps)\1;
%         t_bundle = 100*abs(Z_ub-f_xk);
        t_bundle = 1000*abs(Z_ub-f_xk);
		best_var = variables;
		bestmuVb = muVbx;
		bestmuPg = muPgx;
        f_xk = Z_ub;
    % else
    %     t_bundle = 100*abs(Z_ub-f_xk);
	end
	if Z_ub < upbound
		% update upper bound
		upbound = Z_ub;
		best_contingency = contingency_sol;
		best_base = results0;
        % *** bundle ****
        % variables = variables; % actualizo el mejor caso base
        % *** fin bundle ***
		no_mejora = 0;
	else
		no_mejora = no_mejora+1;
	end
	ZUB1 = [ZUB1;upbound];
	ZUB = [ZUB;Z_ub]; % instruccion para revisar como es realmente el upper bound
	% if it_cuts == 1
	% 	no_mejora = 0;
	% 	% ZUB = [ZUB;Z_ub]; % only for  iterative process
	% 	best_solution_base = results0; % almaceno la mejor solución
	% elseif Z_ub >= min(ZUB)
	% 	no_mejora = no_mejora+1;
	% 	% no_mejora = 0;
	% 	% ZUB = [ZUB;ZUB(end)];
	% else
	% 	no_mejora = 0;
	% 	best_solution = results0;
	% 	% ZUB = [ZUB;Z_ub];
	% end
    % ***** add cuts with proximal bundle ****
	% add cuts
	% resultsa = toggle_cuts_newform(results0,cut_x_contingency);
    resultsa = toggle_cuts_newform_bundle(results0,cut_x_contingency,best_var,t_bundle);
    % *** end add cuts with proximal bundle

	% run opf master problem
	mpopt.verbose = 0;
    % mpopt.opf.start = 0;
	results0 = runopf(resultsa,mpopt);
	
	% fprintf('tiempo con cortes es %f\n',results0.et)
	if results0.success==0 % doesn't converge master problem
		% fprintf('el problema maestro no converge\n')
		% change ipopt options
		mpopt2 = mpopt;
		mpopt2.verbose = 2;
        mpopt2.opf.start = 0;
        mpopt2.ipopt.opts.max_iter = 3000;
% 		mpopt2.ipopt.opts.expect_infeasible_problem = 'yes';
		mpopt2.ipopt.opts.print_level =3;
        mpopt2.ipopt.opts.start_with_resto = 'no';
		mpopt2.ipopt.opts.mu_strategy='adaptive';
        mpopt2.ipopt.opts.adaptive_mu_globalization='kkt-error';
        mpopt2.ipopt.opts.mu_oracle='quality-function';
        mpopt2.ipopt.opts.fixed_mu_oracle = 'quality-function';
%         mpopt2.ipopt.opts.alpha_for_y = 'min';
		mpopt2.ipopt.opts.acceptable_tol = 1e3;
		mpopt2.ipopt.opts.acceptable_constr_viol_tol = 1e-4;
		mpopt2.ipopt.opts.tol = 1e-4;
		% rerun
		results0 = runopf(resultsa,mpopt2);
		% fprintf('la nueva solucion es %i\n',results0.success)
		if ~results0.success
			error('mato proceso por no convergencia del maestro\n')
			% break;
		else
			% fprintf('funciona correctamente\n');
			% pause;
		end
	end
	

	mpopt.verbose = 0;
	Z_lo = results0.f; % cost solution of master problem
	lobound = Z_lo;
	ZLB = [ZLB;lobound];
	z_cut_obj = 0; % for substract contingency cuts cost
	for con = 1:numel(cut_x_contingency)
		z_cut_obj = z_cut_obj + results0.qdc.(strcat('z_cut_con_obj',num2str(con)));
	end
	% calculate basecost for upper bound and cost accept for while loop
	
	basecost = Z_lo - z_cut_obj - results0.qdc.BUNDLE; % creo que falta eliminar costo del bundle pero dejemoslo por el momento así
	cost_accept = 0.1*basecost;
	Pgen_base = [Pgen_base,results0.gen(:,2)];
	% Qgen_base = [Qgen_base,results0.gen(:,3)];
	% Vbus_base = [Vbus_base,results0.bus(:,8)];
	variables = results0.var;
	results0 = rmfield(results0,{'userfcn','var'}); % remove user equations and variables
	results0 = toggle_shuntgen(results0);
	results0 = toggle_softlims_ext1(results0,'on'); % add slack variables again
    it_cuts = it_cuts+1;
	delta_bundle = f_xk - Z_lo; % actualizo deltak de los parametros del bundle
end
toc
plot(1:length(ZLB),ZLB,1:length(ZUB),ZUB)
set(gca,'XTick',1:length(ZLB));
xlabel('Iteraciones')
ylabel('Valor [$]')
% impresion del valor
fprintf('porcentaje de solucion es %f%%\n',(Z_ub-Z_lo)/Z_lo*100);
% fprintf('la solucion de Pg es\n');
% disp(Pgen_base);