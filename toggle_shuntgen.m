function mpc = toggle_shuntgen(mpc)
% funcion para introducir las variables de los shunt variables en el balance nodal
% Modela los shunt como fuentes de potencia reactiva de acuerdo a la siguiente expresion
% \underline{b_i^sc}*v_i^2 <= qshunt_i <= \overline{b_i^sc}*v_i^2

% Comprobar que exista la estructura que contiene la información
% por el momento no se necesita nada más de ninguna conversión

if ~ (isfield(mpc,'comp')&& (size(mpc.comp,1)>=1) && any(mpc.comp(:,5)~=0))
    mpc = mpc;
    mpc.isshunt = false; % no tiene shunts
else
    % agregar las funciones que cumplen con el toggle
    mpc = add_userfcn(mpc,'ext2int',@userfcn_shuntgen_ext2int); % conversion
    mpc = add_userfcn(mpc,'formulation',@userfcn_shuntgen_formulation); % formulacion
    mpc = add_userfcn(mpc,'int2ext',@userfcn_shuntgen_int2ext); % devolver �ndices
    mpc = add_userfcn(mpc,'printpf',@userfcn_shuntgen_printpf); % impresion de resultados
    mpc = add_userfcn(mpc,'savecase',@userfcn_shuntgen_savecase); % guardar resultados
    mpc.isshunt = true;
end

function mpc = userfcn_shuntgen_ext2int(mpc,mopt,args)
    % funcion para convertir de externo a interno los campos y asociar las variables respectivas 
    % se realiza el mapa para saber cuales son los nodos donde se encuentran los shunts.
    % primero se sacan los que están fuera de operación
    % segundo los que no pertenecen a ningún nodo en especifico
    o = mpc.order; % depronto se necesite
    [busshunt]=ismember(mpc.comp(:,1),o.bus.i2e);
    mpc.order.ext.busshunt = busshunt & (mpc.comp(:,5)~=0); % shunts que entran al despacho
    [~, busextloc]=ismember(mpc.comp(:,1),o.ext.bus(:,1)); % se sabe a que nodo pertenece cada shunt en externa y es más fácil calcular todo
    mpc.order.ext.busextloc = busextloc;

function om=userfcn_shuntgen_formulation(om,opt,args)
    % formulacion del modelo
    % crea un vector de variables para la potencia reactiva inyectada por los shunts
    % primero: realizar el mapa 
    mpc = om.get_mpc(); % estructura en numeracion interna
    o = mpc.order; % se requiere por cambio de numeracion
    [busshunt shuntloc] = ismember(mpc.comp(:,1),o.bus.i2e);
    nshunt = size(mpc.comp(busshunt & mpc.comp(:,5)~=0,1),1); % numero shunts a tener en cuenta
    vinicial = (mpc.comp(busshunt & mpc.comp(:,5)~=0,3)+mpc.comp(busshunt & mpc.comp(:,5)~=0,4))/(2*mpc.baseMVA);
    om.add_var('qshunt',nshunt,vinicial);
    % cambiar potencias a valores en pu para procesamiento
    mpc.comp(:,3:4) = mpc.comp(:,3:4)/mpc.baseMVA;
    % introducir restricciones involucrando la tensión
    % toca crear f, df, d2f
    fcn = @(x)limqshunt(x,mpc.comp,busshunt,shuntloc); % calcula f y df
    hess = @(x,lam)limqshunt_hess(x,lam,busshunt,shuntloc,mpc.comp); % computa la matriz hessiana 
    om.add_nln_constraint('qshuntlims',2*nshunt,0,fcn,hess,{'qshunt','Vm'});
    % no tiene costos asociados
function results = userfcn_shuntgen_int2ext(results,mpopt,args)
    % generación y b de los shunts
    results.shuntgen = zeros(size(results.comp,1),2); % crea un vector para almacenamiento de resultados Q, b
    results.shuntgen(results.order.ext.busshunt,1) = results.var.val.qshunt*results.baseMVA;
    results.shuntgen(:,2) = results.shuntgen(:,1)./results.bus(results.order.ext.busextloc,8).^2;
    results.shuntgen(results.shuntgen(:,2)==Inf,2)=0; % por si sufre alguna indeterminacion

function results = userfcn_shuntgen_printpf(results,fd,mpopt,args)
    % tampoco hace nada
function mpc = userfcn_shuntgen_savecase(mpc,fd,prefix,args)
    % tampoco hace nada porque no quiero guardar el caso con nada adicional
% ******************************************************
%                   Funciones adicionales para la formulacion
% ******************************************************
function [g,dg]=limqshunt(x,comp,busshunt,shuntloc)
    % funcion para cálculo de g y de la matriz jacobiana
    % ingresan dos conjuntos en x, las variables Qshunt y Vm de los nodos
    [qshunt,vm]=deal(x{:});
    shuntloc1 = shuntloc(shuntloc~=0 & comp(:,5)~=0);
    ng = length(qshunt);
    nb = length(vm);
    g = [qshunt - comp(busshunt & comp(:,5)~=0,3).*vm(shuntloc1).^2;...
        comp(busshunt & comp(:,5)~=0,4).*vm(shuntloc1).^2-qshunt];
    % crear matrices
    Bmax = sparse([1:length(qshunt)],shuntloc1,comp(comp(:,5)~=0,3),ng,nb);
    Bmin = sparse([1:length(qshunt)],shuntloc1,comp(comp(:,5)~=0,4),ng,nb);
    dg = sparse([speye(ng) -2*Bmax*diag(vm); -speye(ng) 2*Bmin*diag(vm)]);
%     dg = sparse(nb+ng,2*ng);
% funcion para el calculo de la hessiana
function d2g = limqshunt_hess(x,lam,busshunt,shuntloc,comp)
    [qshunt,vm] = deal(x{:}); % extraigo las variables de las restricciones
    ng = length(qshunt);
    nb = length(vm);
    shuntloc1 = shuntloc(shuntloc~=0 & comp(:,5)~=0);
    Bmax = sparse([1:length(qshunt)],shuntloc1,comp(comp(:,5)~=0,3),ng,nb);
    Bmin = sparse([1:length(qshunt)],shuntloc1,comp(comp(:,5)~=0,4),ng,nb);
    lam1 = sparse(shuntloc1,[1:length(qshunt)],lam(1:ng),nb,ng);
    lam2 = sparse(shuntloc1,[1:length(qshunt)],lam(ng+1:2*ng),nb,ng);
    d2g = [sparse(ng,ng) sparse(ng,nb); sparse(nb,ng) -2*lam1*Bmax+2*lam2*Bmin];



    







