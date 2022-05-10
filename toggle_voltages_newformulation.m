function mpc = toggle_voltages_newformulation(mpc)
% TOGGLE para insertar redespacho de reactiva en contingencia
% 
% Vdev: estructura que contiene algunos campos que se deben tener en cuenta
% Vi_base: vector con todas las tensiones nodales de los resultados base
% G_k: conjunto de generadores que pueden redespachar reactivos, todos los generadores.
%
% ----------------------------------------------------------
%                       Rutina para conocer la direccion
% ----------------------------------------------------------




%% Chequeo de los parametros de entrada
if ~(isfield(mpc,'Vdev') && isstruct(mpc.Vdev) ...
        && isfield(mpc.Vdev,'V_i'))
    error('Falta informacion para integrar el modelo');
end
% estructura auxiliar para las variables
if ~(isfield(mpc,'auxvar') && isfield(mpc.auxvar,'vbin'))
    mpc.auxvar.vbin = repmat([0,1],size(mpc.gen,1),1);
    deltavgk2 = (mpc.gen(:,3)-mpc.gen(:,4))./(mpc.gen(:,5)-mpc.gen(:,4));
    deltavgk2((deltavgk2<0) | isnan(deltavgk2)) = 0;
    mpc.auxvar.deltavgk = [ones(size(mpc.gen(:,1))),deltavgk2,zeros(size(mpc.gen(:,1)))];
end
%% Llamada a los callback para modificar el modelo del OPF
mpc = add_userfcn(mpc,'ext2int',@userfcn_vdev_ext2int); % conversion
mpc = add_userfcn(mpc,'formulation',@userfcn_vdev_formulation); % formulacion
mpc = add_userfcn(mpc,'int2ext',@userfcn_vdev_int2ext); % devolver �ndices
mpc = add_userfcn(mpc,'printpf',@userfcn_vdev_printpf); % impresion de resultados
mpc = add_userfcn(mpc,'savecase',@userfcn_vdev_savecase); % guardar resultados
%% Cuerpos de las funciones a usar en el modelo
%----------------- ext2int
function mpc=userfcn_vdev_ext2int(mpc,mopt,args)
% Funcion para preparar los datos de entrada al modelo
% De acuerdo al manual de matpower se debe implementar
% Genera los vectores con la numeraci�n adecuada
% Inicializacion
o = mpc.order; % esta estructura la genera la funci�n ext2int es para
% M= 1e3; % valor de constante para evitar que se presenten cambios en la F.O.
% mapear entre una numeraci�n y otra por cambio en salida de l�neas y en
% generadores aislados y todo lo que genere no activaci�n de elementos
%% Verificaciones iniciales
% sirve para comprobar tamanhos de vectores y esas cosas
nb0 = size(o.ext.bus,1); % longitud de tensiones original
ng = size(o.ext.gen,1); % generadores original, aunque necesito el orden interno
if(length(mpc.Vdev.V_i)~=nb0) % si los vectores de tensiones difieren
    error('tamanho de vectores no es el mismo en ext2int');
end

% if(~any(ismember(mpc.Vdev.G_k,mpc.gen(:,1))))
    % error('no se encuentran generadores en la base de datos');
% end
%% Generar vector de cambios
% en esta seccion se genera el vector que tiene que ver con los nodos a los
% cuales se les puede modificar la tension.
% Falta implementar la comprobaci�n de generadores no activos y por ende
% sin control
% Paso 0: cambiar los costos
% mpc.gencost(:,1)=2; % conversion a polinomial de la forma de la F.O.
% mpc.gencost(:,4)=1; % considera un modelo de costos fijos
% mpc.gencost(:,5)=0; % los costos de los generadores son 0
% Paso 1: generar vector logico
% rowindex = (1:ng)'; % almacena posicion en las filas
% colindex = zeros(ng,1); % almacena posicion en las columnas
% val_row_col = ones(ng,1); % significa que el gen i est� conectado al nodo j
% for i=1:ng % recorre los generadores ya ordenados
% 	colindex(i) = find(o.ext.bus(:,1)==o.ext.gen(i,1));
% end
% mpc.Vdev.genbusmatrix = sparse(rowindex,colindex,val_row_col,ng,nb0);
% mpc.Vdev.vlog = ismember(mpc.bus(:,1),mpc.gen(mpc.gen(:,8)~=0,1)); % unicamente generadores activos
mpc.Vdev.qmin = mpc.gen(:,5)/mpc.baseMVA; % no creo que se necesite.
mpc.Vdev.qmax = mpc.gen(:,4)/mpc.baseMVA; % no creo que se necesite.
% Paso 2: convertir a orden interno
mpc.Vdev.Vlog = ismember(o.ext.bus(:,1),o.ext.gen(:,1)); % para tramites con la tension
mpc = e2i_field(mpc,{'Vdev','V_i'},'bus'); % orden de tensiones
mpc = e2i_field(mpc,{'Vdev','Vlog'},'bus'); % orden logico en interno
% mpc = e2i_field(mpc,{'Vdev','genbusmatrix'},'gen',1);
% mpc = e2i_field(mpc,{'Vdev','genbusmatrix'},'bus',2);
% mpc = e2i_field(mpc,{'Vdev','qmin'},'gen');
% mpc = e2i_field(mpc,{'Vdev','qmax'},'gen');

% mpc = e2i_field(mpc,{'Vdev','vlog'},'bus'); % nodos activos
% se necesita generar una matriz de incidencia gen -Vbus
% conversión de la estructura auxiliar a numeracion interna
mpc.auxvar.vbin = e2i_data(mpc,mpc.auxvar.vbin,'gen',1);
mpc.auxvar.deltavgk = e2i_data(mpc,mpc.auxvar.deltavgk,'gen',1);


%----------------- formulacion
function om = userfcn_vdev_formulation(om,opt,args)
% En esta funcion se hace cambios en la formulacion del modelo
% Unicamente se considera inclusion de restricciones
% si hay generadores que pertenecen a la misma barra deben cumplir la misma potencia
% puede que eso esté causando problemas en la convergencia del otro
mpc = om.get_mpc(); % obtener la estructura de datos del objeto opf
% implementacion de dos restricciones
V0 = mpc.Vdev.V_i; % tensiones iniciales de todas las barras (ya en orden interno)
qmin = mpc.Vdev.qmin;
qmax = mpc.Vdev.qmax;
nvar = size(mpc.gen,1); % devuelve nro de filas de generadores en interno
% Esta formulacion implementa el siguiente modelo
% Vgk = Vgk_min +(Vg-Vgk_min)deltavgk1 + (vgk_max-Vg)deltavgk3
% Qgk = Qgk_max + (Qgk_min-Qgk_max)*deltvgk2
% deltavgk3 <= vbingk1 <= deltavgk2 <= vbingk2 <= deltavgk1 
% 0 <= delgvgki <= 1
% vbingki in {0,1} variable binaria

% ----------------------- Variables binarias --------------------
% son solo dos variables binarias
% variable binaria vbingk1
% om.add_var('vbin_1',length(mpc.gen(:,1)),mpc.auxvar.vbin(:,1),zeros(nvar,1),ones(nvar,1));
% variable binaria vbingk2
% om.add_var('vbin_2',length(mpc.gen(:,1)),mpc.auxvar.vbin(:,2),zeros(nvar,1),ones(nvar,1));
% *****************************

% ------------- variables continuas representacion estados -------
% variables para cumplimiento de tensiones
% pienso que dos es un buen límite pero para mejorar es tomar el límite de la variable.
om.add_var('deltavgk_1',length(mpc.gen(:,1)),mpc.auxvar.deltavgk(:,1),zeros(nvar,1),ones(nvar,1)); % delta para cuando Vgk>=Vg
om.add_var('deltavgk_2',length(mpc.gen(:,1)),mpc.auxvar.deltavgk(:,2),zeros(nvar,1),ones(nvar,1)); % delta para cuando Vgk = Vg
om.add_var('deltavgk_3',length(mpc.gen(:,1)),mpc.auxvar.deltavgk(:,3),zeros(nvar,1),ones(nvar,1)); % delta para cuando Vgk <= Vg
% *********************************
% ---------------- Matriz dispersa para ir de Vb a Vg -----------
nb = size(mpc.bus,1); % numero de barras
ng = size(mpc.gen,1); % numero de generadores
[~,col_index] = ismember(mpc.gen(:,1),mpc.bus(:,1)); % posiciones de tensiones de generadores en vector de nodos
% crear matriz de posiciones
% voy a suponer que todos los generadores tienen nodo asociado
% esto lo da la explicación del ext2int
row_index = [1:ng]'; % numero de fila
val_row_col = ones(ng,1);
Vbusgen = sparse(row_index,col_index,val_row_col,ng,nb);
V0_matrix = sparse(diag(V0(col_index))); % tensiones de nodos matricial
gengen = speye(ng); % matriz identidad dispersa igual al numero de generadores
Busvmin = mpc.bus(col_index,13); % tension minima del nodo (revisar si es en numeración interna)
Busvmax = mpc.bus(col_index,12); % tension maxima del nodo (revisar si esta en numeracion interna)
%% ------------------ ec para la tension
% Vgk = Vgk_min +(Vg-Vgk_min)*deltavgk1 +(Vgk_max-Vg)*deltavgk3
% matriz de las variables
W = [Vbusgen -(V0_matrix-diag(Busvmin)) -(diag(Busvmax)-V0_matrix)];
varset = {'Vm','deltavgk_1','deltavgk_3'};
om.add_lin_constraint('Volt_newform',W,Busvmin,Busvmin,varset);
%% ----------------- ec para la potencia reactiva
% Qgk = Qgk_max +(Qgk_min-Qgk_max)*deltavgk_2
W = [gengen -diag(qmin-qmax)];
om.add_lin_constraint('qgk_newform',W,qmax,qmax,{'Qg','deltavgk_2'});
%% ----------------- ec general
% deltavgk3 <= vbingk1 <= deltavgk2 <= vbingk2 <= deltavgk1 
% deltavgk3 <= vbingk1
om.add_lin_constraint('voltgeneral_1',speye(nvar),[],mpc.auxvar.vbin(:,1),{'deltavgk_3'});
% - deltavgk2 <= -vbingk1
om.add_lin_constraint('voltgeneral_2',-speye(nvar),[],-mpc.auxvar.vbin(:,1),{'deltavgk_2'});
% deltavgk2 <= vbingk2
om.add_lin_constraint('voltgeneral_3',speye(nvar),[],mpc.auxvar.vbin(:,2),{'deltavgk_2'});
%  - deltvgk1 <= -vbingk2
om.add_lin_constraint('voltgeneral_4',-speye(nvar),[],-mpc.auxvar.vbin(:,2),{'deltavgk_1'});
%% ----------------------- Ecuaciones MPEC
% fcn = @(x)bin_v_restraint(x,mpc.taoV); % g y dg
% hess = @(x,lam)bin_v_hessrestraint(x,lam); % defincion de hessiana
% om.add_nln_constraint('bin_v',2*size(mpc.gen,1),0,...
% fcn,hess,{'vbin_1','vbin_2'}); % agregar restriccion al modelo
%% ----------------------- Modelo en funcion objetivo
% se propone un modelo en funcion objetivo para resolver el sistema
% resolver el sistema M*\sum x(1-x)
% se modela de manera cuadrática porque es mucho más fácil
% Q = -2*speye(2*length(mpc.gen(:,1)))/mpc.taoV;
% C = ones(2*length(mpc.gen(:,1)),1)/mpc.taoV;
% K = 0;
% om.add_quad_cost('Vgkcost',Q,C,K,{'vbin_1','vbin_2'});
%----------------------- int2ext  --------------------------------
function results = userfcn_vdev_int2ext(results,  mpopt, args)
% De aca salen los vectores para conocer la solucion
% results.vbin = zeros(length(results.gen(:,1)),3); % matriz para almacenar vbin_i

%results.qgk = [results.var.val.Qgk_1,results.var.val.Qgk_2,results.var.val.Qgk_3];
% results.vbin = [results.var.val.vbin_1,results.var.val.vbin_2]; % lleva las variables binarias
results.deltavgk = [results.var.val.deltavgk_1,results.var.val.deltavgk_2,results.var.val.deltavgk_3];
% results.vgk = zeros(length(results.bus(:,1)),3); % matriz para almacenar Vgk
% results.vgk = [results.var.val.Vgk_1,results.var.val.Vgk_2,...
%     results.var.val.Vgk_3];
aux1 = zeros(size(results.order.ext.gen,1),2); % tamaño de generadores original
% results.qgk = i2e_data(results,results.qgk,aux1,'gen',1);
% results.vbin = i2e_data(results,results.vbin,aux1,'gen',1);
aux2 = zeros(size(results.order.ext.gen,1),3);
results.deltavgk = i2e_data(results,results.deltavgk,aux2,'gen',1);
% convertir nuevamente a orden externo las tensiones que participan
% salida estructura auxiliar
results.auxvar.deltavgk = results.deltavgk;
% results.auxvar.vbin = results.vbin;

% -------------------- printf 
%no hace nada
function results = userfcn_vdev_printpf(results, fd, mpopt, args)
% print results
% disp('no imprimir nada');

% ----------------- savecase
% no hago nada tampoco
% ---------------- Funciones hiperesfera de radio 1
function mpc = userfcn_vdev_savecase(mpc,fd,prefix,args)
    % disp('no se si sirva');

% function [g, dg] = voltpoly_fcnhiper(x)
% % x: vector con las variables ygk
% [vbin1, vbin2, vbin3] = deal(x{:}); % extraer las variables
% m = length(vbin1);
% g = sum([vbin1 vbin2 vbin3].^2,2)-ones(m,1);
% a = 2*[vbin1;vbin2;vbin3];
% posrows = repmat((1:m)',3,1);
% aux = (1:3:3*m)';
% poscols = [aux;aux+1;aux+2];
% dg = sparse(posrows,poscols,a,m,3*m);
% % funcion para la matriz hessiana
% function d2g = voltpoly_hesshiper(x,lam)
% nr = length(x{1});
% a = 2*repelem(lam,length(x));
% posrows=(1:length(a))';
% poscols = (1:length(a))';
% d2g = sparse(posrows,poscols,a,length(a),length(a));
% function [g,dg] = bin_v_restraint(x,taoV)
%     % x: vector con las variables binarias
%     [vbin1, vbin2] = deal(x{:}); % extraer las variables
%     BinVvar = [vbin1;vbin2]; % almacena las tres variables en un solo vector
%     m = length(vbin1);
%     % generar restriccion cuadratica
%     % matH = speye(3*m); % generar matriz H
%     % cvector = ones(3*m,1); % genera vector
%     g = BinVvar - BinVvar.^2-taoV*ones(size(BinVvar)); % restriccion y*(y-1)<=tao
%     dg = sparse(diag(ones(length(BinVvar),1)-2*BinVvar)); % calculo de la matriz jacobiana
% % funcion para calcular la matriz hessiana
% function d2g = bin_v_hessrestraint(x,lam)
%     % basicamente es una matriz constante de 2
%     m = length(x{1}); % calculo la longitud del vector inicial
%     % lam posee las dimensiones de las restricciones
%     a = -2*lam; % vector a ubicar en la matriz dispersa
%     posicion = (1:2*m)'; % posicion de los multiplicadores
%     d2g = sparse(posicion,posicion,a,length(a),length(a));
% **************FIN