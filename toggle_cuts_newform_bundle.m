
function mpc = toggle_cuts_newform(mpc,cuts,variables,tao)
%TOGGLE_CUTS Agrega los cortes al modelo
% Recibe como argumentos la estructura del caso base y la estructura del
% corte.
% bundle: recibe los mejores valores de las variables hasta el momento
% La idea es crear un corte del tipo mu^T(X*-x)
% X*: valores de las variables de la contingencia
% x: variables del caso base
% Verificacion de parametros de entrada
if nargin<3
    error('Se necesitan 3 argumentos de entrada');
end
%% Chequeo de los parametros de entrada
if ~isstruct(cuts)
    error('Para generar el corte debe ingresar como estructura el modelo');
end
mpc.var1 = variables; % reemplazar y mantener (hace parte del bundle)
mpc.tao = tao;
%% Llamada a los callback para modificar el modelo del OPF
mpc = add_userfcn(mpc,'ext2int',@userfcn_addcut2_ext2int,cuts); % conversion de interna a externa
mpc = add_userfcn(mpc,'formulation',@userfcn_addcut_formulation); % Ingreso del corte
mpc = add_userfcn(mpc,'int2ext',@userfcn_addcut_int2ext); % Cambio de interno a externo en lo que se necesite
mpc = add_userfcn(mpc,'printpf',@userfcn_addcut_printpf); % impresion de resultados
mpc = add_userfcn(mpc,'savecase',@userfcn_addcut_savecase); % guardar resultados
%% Cuerpos de las funciones a usar en el modelo
%----------------- ext2int
function mpc=userfcn_addcut2_ext2int(mpc,mopt,args)
% Funcion para preparar los datos de entrada al modelo
% De acuerdo al manual de matpower se debe implementar
% Genera los vectores con la numeracion interna de acuerdo al conjunto sobre el que debe actuar
% Inicializacion
mpc.cuts = args;

for i = 1:numel(args)
    mpc.cuts(i).mu_Pg = e2i_data(mpc,args(i).mu_Pg,'gen',1); % change to internal order Pg
    mpc.cuts(i).mu_volt = e2i_data(mpc,args(i).mu_volt,'bus',1); % change to internal order Vb
end

%----------------- formulacion
function om = userfcn_addcut_formulation(om,opt,args)
% En esta funcion se hace cambios en la formulacion del modelo
% Agrega los cortes al caso base
mpc = om.get_mpc(); % obtener la estructura de datos del objeto opf
cuts = mpc.cuts; % obtengo la totalidad de cortes
%-------INCLUSION DE LOS CORTES (ITERATIVO) DE CADA CONTINGENCIA----------
for i = 1:numel(cuts) % para cada contingencia se agregan los cortes que se llevan hasta el momento
    % nb = size(mpc.bus,1); % numero de barras
    % ng = size(mpc.gen,1); % numero de generadores
    % [~,col_index] = ismember(mpc.gen(:,1),mpc.bus(:,1));
    % row_index = [1:length(mpc.gen(:,1))]'; % numero de fila
    % val_row_col = ones(ng,1);
    % Vbusgen = sparse(row_index,col_index,val_row_col,ng,nb);
    namecut_var = strcat('z_cut_','con_',num2str(i));
    namecut_restr = strcat('restr_con_',num2str(i));
    namecut_obj = strcat('z_cut_con_obj',num2str(i));
    om.add_var(namecut_var,1,0,0,Inf); % add cut variable
    % om.add_var(namecut_var,1,cuts(i).f_obj(end),0,Inf); % add cut variable
    Q = 0;
    C = 1;
    K = 0;
    om.add_quad_cost(namecut_obj,Q,C,K,{namecut_var});
    tao = 1.00;
    % mu_Pg = mu_Pg + cuts(i).mu_Pg;
    % mu_Vb = mu_Vb + cuts(i).mu_volt;
    % f_obj = f_obj + cuts(i).f_obj;
    % constant = cuts(i).constant + constant;
    W = [ones(size(cuts(i).f_obj)), cuts(i).mu_Pg',cuts(i).mu_volt'];
    om.add_lin_constraint(namecut_restr,W,tao*(cuts(i).f_obj+cuts(i).constant),...
        [],{namecut_var,'Pg','Vm'});
    % W = [ones(size(cuts(i).f_obj)), cuts(i).mu_var'];
    % om.add_lin_constraint(namecut_restr,W,cuts(i).f_obj+cuts(i).mu_b,...
    %     [],{namecut_var,'Pg'});
    % end
    % add cut alpha > = w* + lambda_p(f(Pg)-Pg) +lambda_v(f(Vg)-Vg)
    
end
% penalized bundle method (only Pg + Vm)
try
    t = mpc.tao; % falta crear rutina de actualizacion
    Q = 2*t*speye(length(mpc.var1.val.Pg)+length(mpc.var1.val.Vm));
    C = -2*t*[mpc.var1.val.Pg;mpc.var1.val.Vm];
    K = t*sum([mpc.var1.val.Pg;mpc.var1.val.Vm].^2);
    om.add_quad_cost('BUNDLE',Q,C,K,{'Pg','Vm'});
catch
    warning('ya existe funcion objetivo');
end


%---------------- int2ext
function results = userfcn_addcut_int2ext(results,  mpopt, args)
% no es necesario realizar llamados externos e internos
% disp('\n no hacer nada\n');
results = rmfield(results,'cuts'); % elimino cosas que no se necesitan.

function results = userfcn_addcut_printpf(results, fd, mpopt, args)
% fprintf(fd,'\n===========================================================')
% fprintf(fd,'\n                 ***FUnciona***                    ')
% fprintf(fd,'\n===========================================================\n')

function mpc = userfcn_addcut_savecase(mpc,fd,prefix,args)
% disp('no se si sirva');
