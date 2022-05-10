function mpc = toggle_Pgk_newform2(mpc)
%TOGGLE_TEST Extiende el OPF con las ecuaciones de redespacho de potencia activa
% Este modelo usa una nueva representacion del sistema no lineal
% Debe ser de mayor robustez
% En el  MPC debe adicionarse un vector (alfa) con los factores de
% participacion, aunque por la estructura del problema este se encuentra
% como un parametro dentro de los generadores.
% delta: valor del delta, en algun punto se usa
% lowlimit: variable que almacena el limite inferior de ygk_1, se debe
% pasar por ext2int
% if nargin>1
%     mpc.ModifDelta = delta;
%     mpc.lowlimit = lowlimit;
% else
%     mpc.ModifDelta = 0;
%     mpc.lowlimit = zeros(length(mpc.gen(:,1)),1); 
% end
%% Chequeo de los parametros de entrada
% 17-08-2020: cambio en el modelo para incluir un esquema heurístico con la generación de columnas
% se reemplazan las variables binarias por un esquema de columnas y variables de selección de cada columna
% para almacenar la columna se crea una matriz (auxvar.ygk_col) de tamanho gen x 2#iteraciones. 
% Es mucho mas facil mantener (hasta se puede con matrices dispersas) por memoria
if ~(isfield(mpc,'Pgk') && isstruct(mpc.Pgk) ...
        && isfield(mpc.Pgk,'Pg0') && isfield(mpc.Pgk,'G_k'))
    error('Falta informacion de entrada al modelo');
end
% iniciación de las variables 
if ~(isfield(mpc,'auxvar'))
    % si no existe el campo se deben iniciar las variables
    % creacion de estructura dentro de los resultados para almacenar los valores iniciales
    if ~(isfield(mpc,'Pg_out')) % falta evaluar si es cero        
        mpc.auxvar.deltak = 0;
        deltagk2 = (mpc.gen(:,2)-mpc.gen(:,10))./(mpc.gen(:,9)-mpc.gen(:,10));
        deltagk2((deltagk2>1) | isnan(deltagk2))=1;
        mpc.auxvar.deltapgk = [ones(size(mpc.gen(:,1))),...
            deltagk2,zeros(size(mpc.gen(:,1)))];
        mpc.auxvar.ygk_col = repmat([0,1],size(mpc.gen,1),1); % realizar pruebas con otras configuraciones (buscar mapeo de matriz)
        % mpc.auxvar.ygk_col(5,:)=[1,1]; % cambios en un gen
        % mpc.auxvar.colvarP = 1; % tamanho #col x 1
    else
        % obtener los que participan y que se encuentren activos
        Gk_in = ismember(mpc.gen(:,1),mpc.Pgk.G_k);
        Gk_in = Gk_in & (mpc.gen(:,8)>0);
        % estimacion inicial de delta
        delta_ini = (mpc.Pg_out)/sum(mpc.gen(Gk_in,21));
        % conocer los que llegan a su límite
        maxlim = (mpc.gen(Gk_in,9)/mpc.baseMVA - mpc.Pgk.Pg0(Gk_in))./mpc.gen(Gk_in,21) < (delta_ini*ones(sum(Gk_in),1));
        mpc.auxvar.ygk_col = repmat([0,1],size(mpc.gen,1),1); % realizar pruebas con otras configuraciones (buscar mapeo de matriz)
        mpc.auxvar.ygk_col(Gk_in,1) = maxlim;
        Pgkk = mpc.Pgk.Pg0 + mpc.gen(:,21)*delta_ini;
        aux_log = Pgkk*mpc.baseMVA > mpc.gen(:,9);
        Pgkk(aux_log) = mpc.gen(aux_log,9)/mpc.baseMVA;
        deltagk2 = Pgkk./(mpc.gen(:,9)/mpc.baseMVA);
%         deltagk2 = ((Pgkk*mpc.baseMVA - mpc.gen(:,10))./(mpc.gen(:,9)-mpc.gen(:,10))).*mpc.gen(:,21);
        % deltagk2 = (mpc.gen(:,2)-mpc.gen(:,10))./(mpc.gen(:,9)-mpc.gen(:,10));
        deltagk2((deltagk2>1) | isnan(deltagk2))=1;
        % deltagk2(Gk_in) = deltagk2(Gk_in).*(1-maxlim)+maxlim;
        mpc.auxvar.deltapgk = [ones(size(mpc.gen(:,1))),...
            deltagk2,zeros(size(mpc.gen(:,1)))];
        % suma de los que alcanzan el límite más los que no
%         sumpot = sum((mpc.gen(Gk_in,9)+mpc.gen(Gk_in,2)).*maxlim + mpc.gen(Gk_in,2).*(1-maxlim))/mpc.baseMVA;
        sumpot = sum((mpc.gen(Gk_in,9)-mpc.gen(Gk_in,2)).*maxlim)/mpc.baseMVA;
        mpc.auxvar.deltak = (mpc.Pg_out - sumpot)/sum(mpc.gen(Gk_in,21).*(1-maxlim));
        % mpc.auxvar.deltak = max((mpc.gen(Gk_in,9)/mpc.baseMVA - mpc.Pgk.Pg0(Gk_in))./mpc.gen(Gk_in,21));
    end
end
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
 MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
 QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
%% Fijar generaciones que no pertenecen al area
% Se fijan Pmin = Pmax, asi no se mueve
% IPOPT tiene un control para este tipo de variables, pero por el momento no está activo
% Revisar opciones de IPOPT para manejar mejor este problema aunque necesito el multiplicador.
aux1 = ~ismember(mpc.gen(:,GEN_BUS),mpc.Pgk.G_k); % Generadores que no pertenecen al area
aux1 = aux1 | (mpc.gen(:,APF)==0); % Agrego los que el factor de participacion es cero 
mpc.gen(aux1,PMIN)=mpc.Pgk.Pg0(aux1)*mpc.baseMVA; % Cambio de PMIN a Pg0
mpc.gen(aux1,PMAX)=mpc.Pgk.Pg0(aux1)*mpc.baseMVA; % Cambio de PMAX a Pg0
mpc.gencost(:,1)=2; % conversion a polinomial de la forma de la F.O.
mpc.gencost(:,4)=1; % considera un modelo de costos fijos
mpc.gencost(:,5)=0; % los costos de los generadores son 0 (evitar que aparezcan terminos en la optimizacion que no son relevantes en 
% la contingencia)
%% Llamada a los callback para modificar el modelo del OPF
mpc = add_userfcn(mpc,'ext2int',@userfcn_testGP_ext2int); % conversion
mpc = add_userfcn(mpc,'formulation',@userfcn_testGP_formulation); % formulacion
mpc = add_userfcn(mpc,'int2ext',@userfcn_testGP_int2ext); % devolver �ndices
mpc = add_userfcn(mpc,'printpf',@userfcn_testGP_printpf); % impresion de resultados
mpc = add_userfcn(mpc,'savecase',@userfcn_testGP_savecase); % guardar resultados
%% Cuerpos de las funciones a usar en el modelo
%----------------- ext2int
function mpc=userfcn_testGP_ext2int(mpc,mopt,args)
% Funcion para preparar los datos de entrada al modelo
% De acuerdo al manual de matpower se debe implementar
% Genera los vectores con la numeraci�n adecuada
% Inicializacion
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
 MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
 QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
orden = mpc.order; % trae la numeracion de las variables base de matpower (bus,gen,branch)
%M= 1e3; % valor de constante para evitar que se presenten cambios en la F.O.
%% Modificar los vectores a la numeracion interna
% Paso 1: vector l�gico de los generadores que entran a participar
mpc.Pgk.Gin = ismember(orden.ext.gen(:,1),mpc.Pgk.G_k); % con la estructura externa
mpc.Pgk.Gin = mpc.Pgk.Gin & (orden.ext.gen(:,21)~=0); % deja unicamente cuyo factor de participacion sea diferente de cero
% Paso 3: Guardar posiciones de los generadores que participan en el
% despacho
mpc.order.ext.Pgk.G_kpos = find(mpc.Pgk.Gin);
% guardar posiciones del vector de potencias a despachar
mpc.order.ext.Pgk.Pg0 = mpc.Pgk.Pg0;
% Paso 4: Modificar para que coincida con la numeracion interna
% mpc = e2i_field(mpc,{'lowlimit'},'gen');
mpc = e2i_field(mpc,{'Pgk','Gin'},'gen'); % vector de pertenecen o no al area
mpc = e2i_field(mpc,{'Pgk','Pg0'},'gen'); % potencias iniciales
% conversion de la estructura auxiliar a orden ingerno
mpc.auxvar.deltapgk = e2i_data(mpc,mpc.auxvar.deltapgk,'gen',1);
mpc.auxvar.ygk_col = e2i_data(mpc,mpc.auxvar.ygk_col,'gen',1);


% ----------------- formulacion
function om = userfcn_testGP_formulation(om,opt,args)
% En esta funcion se hace cambios en la formulacion del modelo
% con esta etapa se pretende modificar que la funcion objetivo obedezca a
% minimizar las desviaciones de tensi�n sin tener en cuenta las funciones
% de activaci�n sino simplemente inyectar reactiva sin que cambien las

% tensiones
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
 MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
 QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
%-----------------------
mpc = om.get_mpc(); % obtener la estructura de datos del objeto opf
%-------COSTOS----------
% mpc.gencost(:,1)=2; % conversion a polinomial de la forma de la F.O.
% mpc.gencost(:,4)=1; % considera un modelo de costos fijos
% mpc.gencost(:,5)=0; % los costos de los generadores son 0
%% creacion de variables
% los vectores de variables son del tamanho de los generadores a redespachar
% find(mpc.lowlimit.*mpc.Pgk.Gin)
limitsup = ones(sum(mpc.Pgk.Gin),1); % limite superior de las variables binarias
limitslo = zeros(sum(mpc.Pgk.Gin),1); % limite inferior de las variables binarias
deltamax = max(mpc.gen(mpc.Pgk.Gin,APF).\(mpc.gen(mpc.Pgk.Gin,PMAX)-mpc.gen(mpc.Pgk.Gin,PG)))/mpc.baseMVA; % abertura maxima (debe modificarse este l�mite)
deltamin = min(mpc.gen(mpc.Pgk.Gin,APF).\(mpc.gen(mpc.Pgk.Gin,PMIN)-mpc.gen(mpc.Pgk.Gin,PG)))/mpc.baseMVA; % abertura m�nima (debe modificarse este l�mite)
om.userdata.deltamax = deltamax; % almacenar en el modelo de optimizacion la informacion de deltamax
om.userdata.deltamin = deltamin; % almacenar en el modelo de optimizacion la informacion de deltamin
% Modelo para deltak (una sola variable)
om.add_var('DeltaK',1,mpc.auxvar.deltak,deltamin,deltamax); % Deltak que es la variable de salida del modelo
%% Matrices y vectores auxiliares
alphag = mpc.gen(:,APF); % vector de factores de participacion
alphag(~mpc.Pgk.Gin) = 0; % hago cero los que no son del area
% alphag = alphag/sum(alphag); % normalizo los factores de participacion
Pgup = mpc.gen(:,PMAX)/mpc.baseMVA; % limites maximos de generadores en numeracion interna
Pgup(~mpc.Pgk.Gin)=0; % hago cero los que no son del area
Pglo = mpc.gen(:,PMIN)/mpc.baseMVA; % limites minimos de generadores en numeracion interna
Pglo(~mpc.Pgk.Gin)=0; % hago cero los que no son del area
Pg0 = mpc.Pgk.Pg0; % potencias caso base en numeracion interna
Pg0(~mpc.Pgk.Gin) = 0; %hago cero los que no son del area
A1=diag(mpc.Pgk.Gin); % creo una matriz diagonal a partir del vector logico
A2=A1(:,mpc.Pgk.Gin); % elimino las columnas que no entran a participar
Pg0_1 = Pg0(mpc.Pgk.Gin); % unicamente los que entran a participar
alphag_1 = alphag(mpc.Pgk.Gin);
ngen_k = sum(mpc.Pgk.Gin); % numero de generadores que participan en el modelo
% creacion de variables auxiliares del modelo
% delta_ini = (Pglo(mpc.Pgk.Gin)-Pg0_1)./(Pgup(mpc.Pgk.Gin)-Pglo(mpc.Pgk.Gin));
om.add_var('deltapgk_1',ngen_k,mpc.auxvar.deltapgk(mpc.Pgk.Gin,1),limitslo,limitsup);
om.add_var('deltapgk_2',ngen_k,mpc.auxvar.deltapgk(mpc.Pgk.Gin,2),limitslo,limitsup);
om.add_var('deltapgk_3',ngen_k,mpc.auxvar.deltapgk(mpc.Pgk.Gin,3),limitslo,limitsup);
% Modelo de ecuaciones 
% Pgk = Pgmin + (Pgmax-Pgmin)*deltapgk_2
% involucra unicamente los generadores que participan
W1 = diag(Pgup(mpc.Pgk.Gin)-Pglo(mpc.Pgk.Gin)); % diagonalizo para tener todos los generadores
% W1 = W1(:,mpc.Pgk.Gin); % elimino columnas de los que no entran a participar
A3 = speye(size(mpc.gen,1));
A3 = A3(mpc.Pgk.Gin,:);
W = sparse([A3 -W1]);
varset = {'Pg','deltapgk_2'};
om.add_lin_constraint('Pgk',W,Pglo(mpc.Pgk.Gin),Pglo(mpc.Pgk.Gin),varset);
% Deltak = Deltamin + ((Pgmin-Pg)/fpg-Deltamin)*deltapgk_1 + (Pgup-Pglo)deltapgk_2/fpg +(Deltamax - (Pgup-Pg)/fpg)deltapgk_3
varset = {'DeltaK','deltapgk_1','deltapgk_2','deltapgk_3'};
wpgk1 = sparse(diag((Pglo(mpc.Pgk.Gin)-Pg0_1)./alphag(mpc.Pgk.Gin))-deltamin*eye(ngen_k));
wpgk2 = sparse(diag((Pgup(mpc.Pgk.Gin)-Pglo(mpc.Pgk.Gin))./alphag(mpc.Pgk.Gin)));
wpgk3 =sparse(deltamax*eye(ngen_k)-diag((Pgup(mpc.Pgk.Gin)-Pg0_1)./alphag(mpc.Pgk.Gin)));
W = [ones(ngen_k,1) -wpgk1 -wpgk2 -wpgk3];
om.add_lin_constraint('delta_form',W,deltamin*ones(ngen_k,1),deltamin*ones(ngen_k,1),varset);
% Restricciones adicionales del modelo 
% \sum_{q}^{Q} lambdaP_q = 1
% om.add_lin_constraint('sumcolvarP',ones(1,length(mpc.auxvar.colvarP)),1,1,{'colvarP'});
% restricciones de las columnas
% D = [-1,0;1,0;0,-1;0,1]
% primera restriccion deltapgk_3 <= xpgk_1 (fijo)
om.add_lin_constraint('Pgeneral_1',eye(ngen_k),[],mpc.auxvar.ygk_col(mpc.Pgk.Gin,1),{'deltapgk_3'});
% segunda restriccion  -deltapgk_2 <= -xpgk_1 (fijo)
om.add_lin_constraint('Pgeneral_2',-eye(ngen_k),[],-mpc.auxvar.ygk_col(mpc.Pgk.Gin,1),{'deltapgk_2'});
% tercera restriccion deltapgk_2 <= xpgk_2 (fijo)
om.add_lin_constraint('Pgeneral_3',eye(ngen_k),[],mpc.auxvar.ygk_col(mpc.Pgk.Gin,2),{'deltapgk_2'});
% cuarta restriccion -deltapgk_1 <= -xpgk_2 (fijo)
om.add_lin_constraint('Pgeneral_4',-eye(ngen_k),[],-mpc.auxvar.ygk_col(mpc.Pgk.Gin,2),{'deltapgk_1'});
%---------------- int2ext
function results = userfcn_testGP_int2ext(results,  mpopt, args)
    results.deltapgk = zeros(length(results.gen(:,1)),3);
    results.deltapgk(results.Pgk.Gin,:) = [results.var.val.deltapgk_1,results.var.val.deltapgk_2, results.var.val.deltapgk_3];
    % results.ygk(results.Pgk.Gin,1)=results.var.val.ygk_1;
    % results.ygk(results.Pgk.Gin,2)=results.var.val.ygk_2;
    results.Pgk.Pg0 = results.order.ext.Pgk.Pg0;
    aux1 = zeros(size(results.order.ext.gen,1),2);
    aux2 = zeros(size(results.order.ext.gen,1),3);
    % results.ygk = i2e_data(results,results.ygk,aux1,'gen',1);
    results.deltapgk = i2e_data(results,results.deltapgk,aux2,'gen',1);
    results.delta = results.var.val.DeltaK;
    % para enviar los resultados de las variables auxiliares

    results.auxvar.deltak = results.delta;
    results.auxvar.deltapgk = results.deltapgk;
function results = userfcn_testGP_printpf(results, fd, mpopt, args)
    % print results
    % fprintf(fd,'\n===========================================================')
    % fprintf(fd,'\n                 ***Valor del DeltaK***                    ')
    % fprintf(fd,'\n===========================================================')
    % fprintf(fd,'\n%3.10f\n',results.var.val.DeltaK)
    %----------------------------------------------------------------------

function mpc = userfcn_testGP_savecase(mpc,fd,prefix,args)
    % disp('no se si sirva');
    %% ********************************************************************