function [const,vector] = benderscutVb(mpc)
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
% % mapear de generacion a bus
% se obtienen dos vectores LIA: que es logico de buses que están en generación; y 
% BusinGen: que me dice en qué posiciones de gen están los buses de barras
[LIA,BusinGen]= ismember(mpc.bus(:,1),mpc.gen(:,1)); % buscar buses en el vector de generacion
BusinGen = BusinGen(BusinGen~=0); % eliminar los buses que no son de generación
mu_vb = zeros(size(mpc.bus(:,1))); % vector con ceros del tamanho de buses para hacer el mapa
mu_vb(LIA) = mu_x(BusinGen); % almaceno en las posiciones los multiplicadores de cada tension de generacion
% % traigo los deltavgk que se encuentran en generacion
deltavgk = zeros(size(mpc.bus,1),3); % inicio con vector lleno de zeros
deltavgk(LIA,:) = mpc.deltavgk(BusinGen,:); % convierto el vector deltavgk de gen a bus

% evaluacion del vector para multiplicar por el vector de variable de tensiones
vector = mu_vb.*(deltavgk(:,1)-deltavgk(:,3));
const = vector'*mpc.Vdev.V_i;

%% fin de evaluacion new form