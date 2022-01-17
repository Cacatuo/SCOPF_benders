function [const,vector] = benderscutVb1(mpc)
% funcion para calcular el corte de Benders para las tensiones
% ingresa el mpc de la contingencia
% const: parte constante del corte
% vector: termino que multiplica a las tensiones de las barras (de generacion)
% mpc: caso resuelto en la contingencia
% obtener conjunto de participacion
%% modelo para formulacion de new form
mpc.order.state = 'i'; % cambio a orden interno
mu_x = mpc.lin.mu.u.Volt_newform - mpc.lin.mu.l.Volt_newform; % tamanho de generacion interna
aux1 = zeros(size(mpc.gen(:,1)));
aux2 = zeros(size(mpc.bus(:,1)));
%V0 = i2e_data(mpc,mpc.Vdev.V_i,aux2,'bus'); % extrae valor de tensiones iniciales
mu_x = i2e_data(mpc,mu_x,aux1,'gen'); % cambio a generacion externa
mpc.order.state ='e'; % cambio a orden externo
% mapear buses a generación
[LIA,GeninBus1]= ismember(mpc.gen(:,1),mpc.bus(:,1)); % buscar que generadores están en buses
GeninBus = GeninBus1(GeninBus1~=0); %elimino los que son cero
const = mpc.bus(GeninBus,8)+mpc.bus(GeninBus,13).*mpc.deltavgk(LIA,1)-mpc.bus(GeninBus,12).*mpc.deltavgk(LIA,3)-mpc.bus(GeninBus,13);
const = mu_x(LIA)'*const;
vec_aux = mu_x(LIA).*(mpc.deltavgk(LIA,1)-mpc.deltavgk(LIA,3));
nb = size(mpc.bus,1); % numero de buses
ng = size(GeninBus,1); % generadores que se encuentran activos
mat_aux = sparse(GeninBus,[1:ng]',ones(ng,1),nb,ng);
vector = mat_aux*vec_aux; % vector que multiplica las tensiones para el corte
% % mapear de generacion a bus
% se obtienen dos vectores LIA: que es logico de buses que están en generación; y 
% BusinGen: que me dice en qué posiciones de gen están los buses de barras
% [LIA,BusinGen]= ismember(mpc.bus(:,1),mpc.gen(:,1)); % buscar buses en el vector de generacion
% BusinGen = BusinGen(BusinGen~=0); % eliminar los buses que no son de generación
% mu_vb = zeros(size(mpc.bus(:,1))); % vector con ceros del tamanho de buses para hacer el mapa
% mu_vb(LIA) = mu_x(BusinGen); % almaceno en las posiciones los multiplicadores de cada tension de generacion
% % % traigo los deltavgk que se encuentran en generacion
% deltavgk = zeros(size(mpc.bus,1),3); % inicio con vector lleno de zeros
% deltavgk(LIA,:) = mpc.deltavgk(BusinGen,:); % convierto el vector deltavgk de gen a bus

% % % evaluacion de la constante de tension
% const = mpc.bus(:,8)+mpc.bus(:,13).*deltavgk(:,1)-mpc.bus(:,12).*deltavgk(:,3)-mpc.bus(:,13);
% const = mu_vb'*const; % multiplica por el multiplicador para extraer los valores 

% % evaluacion del vector para multiplicar por el vector de variable de tensiones
% vector = mu_vb.*(deltavgk(:,1)-deltavgk(:,3));

% %% fin de evaluacion new form

% %% evaluacion por medio de combinaciones lineales

% % mpc.order.state = 'i'; % cambio a modo interno
% % mu_x = mpc.lin.mu.u.sum_vgk - mpc.lin.mu.l.sum_vgk; % multiplicador en generacion interna
% % aux1 = zeros(size(mpc.gen(:,1)));
% % aux2 = zeros(size(mpc.bus(:,1)));
% % mu_x = i2e_data(mpc,mu_x,aux1,'gen'); % cambio a generacion externa
% % mpc.order.state = 'e';
% % % vectores para hacer mapeo de gen a bus
% % [LIA,BusinGen]= ismember(mpc.bus(:,1),mpc.gen(:,1));
% % BusinGen = BusinGen(BusinGen~=0); % eliminar los buses que no tienen asociados generadores
% % mu_vb = zeros(size(mpc.bus(:,1))); % vector con ceros del tamanho de buses para hacer el mapa
% % mu_vb(LIA) = mu_x(BusinGen); % almaceno en las posiciones los multiplicadores de cada tension de generacion
% % deltavgk = zeros(size(mpc.bus,1),4); % inicio con vector lleno de zeros, es 4 por los numeros que toca almacenar
% % deltavgk(LIA,:) = mpc.deltavgk(BusinGen,:); % convierto el vector deltavgk de gen a bus
% % % parte constante del corte para la tension
% % const = mpc.bus(:,8) - deltavgk(:,1).*mpc.bus(:,13) -deltavgk(:,4).*mpc.bus(:,12);
% % % multiplicacion por el vector cte
% % const = mu_vb'*const;
% % % vector que multiplica por las tensiones
% % vector = mu_vb.*(deltavgk(:,3)+deltavgk(:,2));

% %% evaluacion polihedral
% % El corte es mu1(Vgk1-Vg*ygk1) - mu2(Vgk2-Vg*ygk2) + mu3(Vgk3-Vg*ygk3)
% % cambiar a modo interno
% % mpc.order.state = 'i';
% % aux1 = zeros(size(mpc.gen(:,1))); % vector auxiliar
% % aux2 = repmat(aux1,1,3); % matriz auxiliar
% % mu_x = [mpc.lin.mu.u.Volt1 - mpc.lin.mu.l.Volt1,...
% %     mpc.lin.mu.l.Volt2, mpc.lin.mu.u.Volt3];
% % mu_x = i2e_data(mpc,mu_x,aux2,'gen');
% % mpc.order.state = 'e'; % cambio a modo externo
% % % vectores para hacer mapeo de gen a bus
% % [LIA,BusinGen]= ismember(mpc.bus(:,1),mpc.gen(:,1));
% % BusinGen = BusinGen(BusinGen~=0); % elimino los buses que no tienen asociados generadores

% % % parte constante del corte
% % const = mu_x(:,1)'*mpc.vgk(:,1) - mu_x(:,2)'*mpc.vgk(:,2) + mu_x(:,3)'*mpc.vgk(:,3);
% % % multiplicador de la variable Vbus
% % mu_y = mu_x.*mpc.vbin;
% % mu_y(:,2) = -mu_y(:,2);
% % mu_y = sum(mu_y,2);
% % % mapear de generadores a buses
% % mu_vb = zeros(size(mpc.bus(:,1)));
% % mu_vb(LIA) = mu_y(BusinGen);
% % vector = mu_vb;