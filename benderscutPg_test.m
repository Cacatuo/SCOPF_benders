function [const,vector] = benderscutPg_test(mpc)
% funcion para calcular el corte de Benders por generacion
% ingresa el mpc de la contingencia
% const: parte constante del corte
% vector: termino de potencias que se multiplica por las variables potencia activa
% mpc: caso resuelto en la contingencia
% obtener conjunto de participacion

mpc.order.state = 'i'; % cambio a orden interno
aux1 = zeros(size(mpc.gen(:,1))); % se encuentra en orden externo
Pgk_Gin_ext = i2e_data(mpc,mpc.Pgk.Gin,aux1,'gen');  % trae el vector de los que participan a externo
Pgk_Gin_ext = logical(Pgk_Gin_ext); % para manejar los indices lógicos (se aplica porque no los deja como índices lógicos)
mu = zeros(size(mpc.Pgk.Gin)); % multiplicador en orden interno
mu(mpc.Pgk.Gin) = mpc.lin.mu.u.delta_form - mpc.lin.mu.l.delta_form; % multiplicador de los del área únicamente.
aux2 = mpc.var.mu.u.Pg - mpc.var.mu.l.Pg; % se encuentra en orden interno (Pgk = Pg)
mu(~mpc.Pgk.Gin) = aux2(~mpc.Pgk.Gin);
mu = i2e_data(mpc,mu,aux1,'gen'); % almacena en orden interno los multiplicadores que tienen que ver con Pg
mpc.order.state ='e'; % cambio a orden externo

% % obtener valores en pu
deltak = mpc.delta; % asignacion para llamarlo mucho más fácil (se encuentra en pu)
Pgmin = mpc.gen(:,10)/mpc.baseMVA; % convertir Pmin a pu
Pgmax = mpc.gen(:,9)/mpc.baseMVA; % convertir Pmax a pu
% invalphag = mpc.gen(:,21).\1; % alpha invertido
alphag = mpc.gen(:,21);
% cambiar los multiplicadores a cero de los generadores que se encuentran por fuera de operación
mu(mpc.gen(:,8)==0)=0;
% termino constante para la ecuacion Pgk = Pg (gen que no pertenecen al area) 
const = mu(~Pgk_Gin_ext)'*mpc.gen(~Pgk_Gin_ext,2)/mpc.baseMVA;
% termino constante para los que participan
term1 = (Pgmin./alphag - mpc.om.get_userdata('deltamin')).*mpc.deltapgk(:,1);
term2 = ((Pgmax-Pgmin).*mpc.deltapgk(:,2))./alphag;
term3 = (mpc.om.get_userdata('deltamax')-Pgmax./alphag).*mpc.deltapgk(:,3);
term4 = (deltak-mpc.om.get_userdata('deltamin'))*Pgk_Gin_ext; % 
% suma de términos constantes
% const = const + mu'*(-term1-term2-term3+term4);
% calculo de lo que multiplica a (-Pg) (para los que participan en el area)
term_vector = mu.*(-mpc.deltapgk(:,1)+mpc.deltapgk(:,3))./alphag;
% cambio del vector mu por los nuevos multiplicadores 
mu(Pgk_Gin_ext)= term_vector(Pgk_Gin_ext);
% const = const + mu(Pgk_Gin_ext)'*mpc.Pgk.Pg0(Pgk_Gin_ext);
const = mu'*mpc.Pgk.Pg0;
vector = mu;