function [const,vector] = benderscutPg(mpc)
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
const = const + mu'*(-term1-term2-term3+term4);
% calculo de lo que multiplica a (-Pg) (para los que participan en el area)
term_vector = mu.*(-mpc.deltapgk(:,1)+mpc.deltapgk(:,3))./alphag;
% cambio del vector mu por los nuevos multiplicadores 
mu(Pgk_Gin_ext)= term_vector(Pgk_Gin_ext);
vector = mu;

% para modelo polihedral
% obtener conjunto de participacion
% mpc.order.state = 'i';
% aux1 = zeros(size(mpc.gen(:,1))); % vector auxiliar
% aux2 = repmat(aux1,1,3); % matriz auxiliar
% mu_area = zeros(size(mpc.Pgk.Gin,1),3); % almacenar valores de multiplicador aca estan en orden interno
% % rellenar de acuerdo a las ecuaciones que cumplen
% % primera ecuacion Pgk1 = Pg*ygk1 + alphag*deltagk1
% mu_area(mpc.Pgk.Gin,1) = mpc.lin.mu.u.res2_3 - mpc.lin.mu.l.res2_3;
% % segunda ecuacion Pgmax*ygk2 <= Pg*ygk2 + alphag*deltagk2
% mu_area(mpc.Pgk.Gin,2) = mpc.lin.mu.u.res3_2;
% % tercera ecuacion Pgmin*ygk3 >= Pg*ygk3 + alphag*deltagk3
% mu_area(mpc.Pgk.Gin,3)= mpc.lin.mu.l.res4_2;
% % creacion del vector del multiplicador para la potencia activa (orden interno)
% mu = zeros(size(mpc.Pgk.Gin));
% % rellenar con los que generadores que no participan en el despacho
% mu(~mpc.Pgk.Gin) = mpc.var.mu.u.Pg(~mpc.Pgk.Gin) - mpc.var.mu.l.Pg(~mpc.Pgk.Gin);
% mu = i2e_data(mpc,mu,aux1,'gen'); % cambiar vector a orden externo
% mu_area = i2e_data(mpc,mu_area,aux2,'gen',1); % cambiar vector a orden externo
% Glogic = i2e_data(mpc,mpc.Pgk.Gin,aux1,'gen'); % cambiar vector logico de areas a externo
% Glogic = logical(Glogic); % conversion a logico para mantenerlo asi
% mpc.order.state = 'e'; % devuelve a orden externo

% % valores en pu
% Pgmin = mpc.gen(:,10)/mpc.baseMVA; % convertir Pmin a pu
% Pgmax = mpc.gen(:,9)/mpc.baseMVA; % convertir Pmax a pu
% % hallar el valor del multiplicador de acuerdo a las ecuaciones de participacion
% mux = mu_area.*mpc.ygk;
% mu(Glogic) = mux(Glogic,1)+mux(Glogic,2)-mux(Glogic,3);
% % termino constante del corte
% % constante de los que quedan fijos
% const = mu(~Glogic)'*mpc.gen(~Glogic,2)/mpc.baseMVA;
% % termino de los que participan en el area
% term1 = mu_area(:,1)'*(mpc.pgk_i(:,1)-mpc.gen(:,21).*mpc.deltai(:,1)); % mu1(Pgk1-alphag*deltagk1)
% term2 = mu_area(:,2)'*(Pgmax.*mpc.ygk(:,2)-mpc.gen(:,21).*mpc.deltai(:,2)); % mu2(Pgmax*ygk2 - alphag*deltagk2)
% term3 = mu_area(:,3)'*(Pgmin.*mpc.ygk(:,3)-mpc.gen(:,21).*mpc.deltai(:,3)); % mut3(Pgmin*ygk3 - alphag*deltagk3)
% const = const + term1 + term2 - term3;
% % volver cero la parte del vector que no participa
% mu(mpc.gen(:,8)==0)=0;
% vector = mu;