function out = solvempec_test1(mpc,mpopt)
% Funcion para resolver iterativamente la contingencia
% Se resuelve mediante los multiplicadores de lagrange, dejando fijas las restricciones y encontrando nuevos movimientos luego de una solución
% En caso de ser muy grande la red, que necesite mucho tiempo para resolver un caso se debe cortar de alguna manera
% mpc: caso de matpower a ser resuelto (viene con los toggles)
% out: salida con el caso resuelto
% mpopt: opciones del matpower se deben cargar
% se debe cumplir un criterio de parada.
it_total = 50;
% ipopt options
mpopt.ipopt.opts.max_cpu_time = 3600; % dejarlo en una hora
mpopt.ipopt.opts.acceptable_dual_inf_tol = 1;
mpopt.ipopt.opts.max_iter = 3000;
mpopt.verbose = 0;
mpopt.opf.start = 2;
% mpopt.ipopt.opts.slack_bound_frac = 0.5;
% mpopt.ipopt.opts.slack_bound_push = 0.5;
% mpopt.ipopt.opts.bound_mult_init_val = 1e3;
% mpopt.ipopt.opts.mehrotra_algorithm = 'no';
% mpopt.ipopt.opts.mu_strategy = 'adaptive';
% mpopt.ipopt.opts.mu_oracle = 'loqo';
% mpopt.ipopt.opts.fixed_mu_oracle = 'loqo';
% mpopt.ipopt.opts.max_iter = 10;
% mpopt.ipopt.opts.mu_linear_decrease_factor = 0.05;
% mpopt.ipopt.opts.corrector_type ='affine';
% mpopt.ipopt.opts.soc_method = 1;
% mpopt.ipopt.opts.ma57_automatic_scaling = 'yes';
% mpopt.ipopt.opts.alpha_for_y = 'primal-and-full';
mpopt.ipopt.opts.print_level = 5;
mpopt.ipopt.opts.acceptable_tol = 1e-3;
% mpopt.ipopt.opts.mu_target = 1;
mpopt.acceptable_obj_change_tol = 1;
% establecer un criterio con respecto al caso base
fobj_base = mpc.fobj;
deltafbase = 0.1*mpc.fobj; 
out = runopf(mpc,mpopt); % Solucion OPF contingencia
clear mpc; % elimino mpc para ahorrar espacio en memoria
out1 = out; % para saber si hay cambio (mucho almacenamiento no todo se necesita)
iter = 0; % inicio contador de iteraciones
beta_factor = 0.1; % con este valor controlo reduccion de la F.O. para entrar al ciclo
deltafmin = beta_factor * out.f; % calculo mejora  para algoritmo
% calculo mu (Tensiones y potencia) (u3-u2)
% mu = out.lin.mu.u.Pgeneral_3 - out.lin.mu.u.Pgeneral_2;
% calcular multiplicador (de las restricciones de igualdad)
delta2 = out.var.val.deltapgk_2;
% mu1 = out.lin.mu.u.Pgeneral_1 - out.lin.mu.u.Pgeneral_2;
% mu2 = out.lin.mu.u.Pgeneral_3 - out.lin.mu.u.Pgeneral_4;
mu = out.lin.mu.u.delta_form - out.lin.mu.l.delta_form;
lambda = out.lin.mu.u.Pgk - out.lin.mu.l.Pgk;
mu_p = mu + lambda;
% moving factor: enable changes in some cases, where is neccesary to avoid infeasible sets
sum_pos = 0;
sum_neg = 0;
sum_pos = sum(mu_p(delta2>(1-1e-4) & mu_p>0));
sum_neg = sum(mu_p(delta2<(1e-4) & mu_p<0));
dir_mov = sum_pos > (-1*sum_neg); % da la dirección
dir_mov = repmat(dir_mov,length(mu_p),1);
% movimiento de las tensiones
mu_v = out.lin.mu.u.qgk_newform - out.lin.mu.l.qgk_newform;
delta2_v = out.var.val.deltavgk_2;
% mu_v = out.lin.mu.u.voltgeneral_3 - out.lin.mu.u.voltgeneral_2;
stop_criteria = any(abs([mu_p;mu_v])>deltafmin); % criterio de parada
% fprintf('converge cont %i\n',out.success);
% out.success = 1; % solamente a modo de prueba
% extraer valores iniciales de las variables binarias
% El criterio de otra iteración se compone de:
% stop_criteria: multiplicadores para saber si un cambio mejora la F.O.
% out.success: que la contingencia converja
% iteraciones: que el numero de iteraciones no se supere
% que el valor obtenido sea mayor al 5% de la funcion objetivo del caso base
while stop_criteria && (iter < it_total) && (out.f>=0.05*fobj_base) && out.success
    % proceso de actualizacion de variables binarias (Potencia)
        if ((abs(sum_pos)>deltafmin)||(abs(sum_neg)>deltafmin))
            
%             fprintf('ingresa por heurística de Pgk\n');
%             pause();
            % extraigo los que participan
            ygk_update = out.auxvar.ygk_col(out.Pgk.Gin,:);
            % actualizo si se mueve a la derecha
            ygk_update(delta2>(1-1e-4)&dir_mov,1)=1;
            ygk_update(delta2>(1-1e-4)&dir_mov,2)=1;
            ygk_update(delta2<1e-4 & dir_mov,1)=0;
            ygk_update(delta2<1e-4 & dir_mov,2)=1;
            % actualizo si se mueve a la izquierda
            ygk_update(delta2<1e-4 & ~dir_mov,1)=0;
            ygk_update(delta2<1e-4 & ~dir_mov,2)=0;
            ygk_update(delta2>(1-1e-4) & ~dir_mov,1)=0;
            ygk_update(delta2>(1-1e-4) & ~dir_mov,2)=1;
            % xgk2zero = ygk_update(:,2)==0; % se puede igual porque es un entero
            % xgk1zero = ygk_update(:,1)==0; % se puede igual porque es un entero
            % % xgk2==0 y (u3-u2) > deltafmin cambian xgk2 por 1
            % ygk_update(xgk2zero & (mu > deltafmin) & dir_mov,2)=1;
            % % xgk2 != 0 & mu < - deltafmin & xgk1 ==1 cambian xgk1 a 0
            % ygk_update(~xgk2zero & (mu < -deltafmin) & ~xgk1zero & ~dir_mov,1)=0;
            % % mu < -deltafmin & xgk1 == 0 & xgk2 != 0 cambian xgk2 a cero
            % ygk_update((mu<-deltafmin) & xgk1zero & ~xgk2zero & ~dir_mov ,2) = 0;
            % % mu >deltafmin & xgk1 == 0 & xgk2 != 0 cambian xgk1 a 1
            % ygk_update((mu>deltafmin) & xgk1zero & ~xgk2zero & dir_mov,1) = 1;
            aux = repmat([0,1],size(out.gen,1),1);
            out.order.state = 'i'; % manejo interno por un momento
            out.auxvar.ygk_col = zeros(size(out.order.int.gen,1),2); % igualo tamanho en interno
            out.auxvar.ygk_col(out.Pgk.Gin,:) = ygk_update;
            out.auxvar.ygk_col = i2e_data(out,out.auxvar.ygk_col,aux,'gen',1);
            out.order.state = 'e';
    else
        % Por el momento no hacer nada  
        aux = repmat([0,1],size(out.gen,1),1);
        out.order.state = 'i';
        out.auxvar.ygk_col = i2e_data(out,out.auxvar.ygk_col,aux,'gen',1);
        out.order.state = 'e';
    end

    % proceso actualizacion de variables binarias (tension)
    if (any(abs(mu_v)>deltafmin))
%         fprintf('ingresa por heurística de Vgk\n');
        % vgk2zero = out.auxvar.vbin(:,2)==0; % se puede igualar porque es un entero
        % vgk1zero = out.auxvar.vbin(:,1)==0; % se puede igualar porque es un entero
        % vgk2==0 y (u3-u2) > deltafmin cambian vgk2 por 1
        % cambiar hacia la derecha
        out.auxvar.vbin(delta2_v>(1-1e-4) & mu_v>0,2)=1;
        out.auxvar.vbin(delta2_v>(1-1e-4) & mu_v>0,1)=1;
        out.auxvar.vbin(delta2_v<1e-4 & mu_v>0,2)=1;
        out.auxvar.vbin(delta2_v<1e-4 & mu_v>0,1)=0;
        % cambiar hacia la izquierda
        out.auxvar.vbin(delta2_v>(1-1e-4) & mu_v<0,2)=1;
        out.auxvar.vbin(delta2_v>(1-1e-4) & mu_v<0,1)=0;
        out.auxvar.vbin(delta2_v<1e-4 & mu_v<0,2)=0;
        out.auxvar.vbin(delta2_v<1e-4 & mu_v<0,1)=0;

        % out.auxvar.vbin(vgk2zero & (mu_v > deltafmin),2)=1;
        % % vgk2 != 0 & mu_v < - deltafmin & vgk1 ==1 cambian vgk1 a 0
        % out.auxvar.vbin(~vgk2zero & (mu_v < -deltafmin) & ~vgk1zero,1)=0;
        % % mu_v < -deltafmin & vgk1 == 0 & vgk2 != 0 cambian vgk2 a cero
        % out.auxvar.vbin((mu_v<-deltafmin) & vgk1zero & ~vgk2zero,2) = 0;
        % % mu_v >deltafmin & vgk1 == 0 & vgk2 != 0 cambian vgk1 a 1
        % out.auxvar.vbin((mu_v>deltafmin) & vgk1zero & ~vgk2zero,1) = 1;
        aux = repmat([0,1],size(out.gen,1),1);
        out.order.state = 'i';
        out.auxvar.vbin = i2e_data(out,out.auxvar.vbin,aux,'gen',1);
        out.order.state = 'e';
    else

        aux = repmat([0,1],size(out.gen,1),1);
        out.order.state = 'i';
        out.auxvar.vbin = i2e_data(out,out.auxvar.vbin,aux,'gen',1);
        out.order.state = 'e';
    end
    iter = iter+1; % actualizo # iteraciones
    % simulo modelo nuevamente
    fo_ant = out.f;
    % valor de la iteracion que cambia
%     xchange = sum(out.lin.mu.u.delta_form - out.lin.mu.l.delta_form);
    % fprintf('valor de la heurística nueva %f\n',xchange);
    out = runopf(out,mpopt); % Solucion opf contingencia 
    % fprintf('converge cont %i\n',out.success);
    % actualizar criterio de convergencia
    deltafmin = out.f * beta_factor; 
    mu = out.lin.mu.u.delta_form - out.lin.mu.l.delta_form;
    lambda = out.lin.mu.u.Pgk - out.lin.mu.l.Pgk;
    mu_p = mu + lambda;
    sum_pos = sum(mu_p(delta2>(1-1e-4) & mu_p>0));
    sum_neg = sum(mu_p(delta2<(1e-4) & mu_p<0));
    dir_mov = sum_pos > (-1*sum_neg); % da la dirección
    dir_mov = repmat(dir_mov,length(mu_p),1);
    % mu = out.lin.mu.u.Pgeneral_3 - out.lin.mu.u.Pgeneral_2;
    mu_v = out.lin.mu.u.qgk_newform - out.lin.mu.l.qgk_newform;
    % revision de cambio

    % verificar si existen diferencias entre inicial y calculado
    diff_p = any(any(out1.auxvar.ygk_col-out.auxvar.ygk_col));
    diff_v = any(any(out1.auxvar.vbin-out.auxvar.vbin));
    if ~out.success
        disp('no converge una vez')
        % disp('valores potencia')
        % out.auxvar.ygk_col
        % disp('valores tension')
        % pause;
        % out.auxvar.vbin
        mpopt2=mpopt;
        mpopt2.opf.start = 0; % arranque sin tener en cuenta la sol actual
        mpopt2.ipopt.opts.least_square_init_primal = 'yes';
        mpopt2.ipopt.opts.max_iter = 3000;
        mpopt2.ipopt.opts.start_with_resto = 'no';
        mpopt2.ipopt.opts.mu_strategy='adaptive';
        mpopt2.ipopt.opts.adaptive_mu_globalization='kkt-error';
        mpopt2.ipopt.opts.mu_oracle='quality-function';
        mpopt2.ipopt.opts.fixed_mu_oracle = 'quality-function';
        mpopt2.ipopt.opts.acceptable_tol = 1e3;
		mpopt2.ipopt.opts.acceptable_constr_viol_tol = 1e-4;
		mpopt2.ipopt.opts.tol = 1e-4;
        mpopt2.ipopt.opts.corrector_type='primal-dual';
        % modificar variables auxiliares para verificar que funcione
        % aux = repmat([0,1],size(out.gen,1),1);
        % out.order.state = 'i';
        % out.auxvar.vbin = i2e_data(out,out.auxvar.vbin,aux,'gen',1);
        % out.auxvar.ygk_col = i2e_data(out,out.auxvar.ygk_col,aux,'gen',1);
        % out.order.state = 'e';
        out.auxvar.vbin = repmat([0,1],size(out.gen,1),1);
        out.auxvar.ygk_col = repmat([0,1],size(out.gen,1),1);
        out = runopf(out,mpopt2);
        if ~out.success
            out = out1;
            disp('no converge contingencia')
            break; % sale si no converge
        end
    elseif out.f <= out1.f
        out1 = out;
    end
    stop_criteria = (any(abs([mu_p;mu_v])>deltafmin)|(out.f<((1-beta_factor)*fo_ant))) & diff_p & diff_v;
end
% out = rmfield(out,'auxvar'); % borro campo que no se necesita en el proceso posterior y reduzco la carga computacional
% out = rmfield(out,'holgura'); % borro variables de holgura porque no se necesita para caso base ni en otras contingencias
out = rmfield(out1,'auxvar'); % borro campo que no se necesita en el proceso posterior y reduzco la carga computacional
out = rmfield(out1,'holgura'); % borro variables de holgura porque no se necesita para caso base ni en otras contingencias
if out.f<0.05*fobj_base
    out.f = 0; % tal vez no se necesite pero lo dejo por el momento
    out.lin.mu.u.delta_form = zeros(size(out.lin.mu.u.delta_form));
    out.lin.mu.l.delta_form = zeros(size(out.lin.mu.l.delta_form));
    out.var.mu.u.Pg = zeros(size(out.var.mu.u.Pg));
    out.var.mu.l.Pg = zeros(size(out.var.mu.l.Pg));
    out.lin.mu.u.Volt_newform = zeros(size(out.lin.mu.u.Volt_newform));
    out.lin.mu.l.Volt_newform = zeros(size(out.lin.mu.l.Volt_newform));
end
% fprintf('el num iteraciones  es %f\n',iter);
