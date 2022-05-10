# SCOPF_benders
Security Constrained Optimal Power Flow
Repositorio para almacenar las rutinas del SCOPF resuelto por Benders
El desarrollo está en MATLAB usando MATPOWER y usando las funciones para la expansión del modelo
Creado por Camilo Acosta como proyecto de maestría en ingeniería eléctrica
UTP
2021
Requerimientos:
- MATLAB
- MATPOWER: https://matpower.org/
- IPOPT: https://github.com/ebertolazzi/mexIPOPT
- IPOPT (repositorio usado): https://github.com/jonathancurrie/OPTI

Para su funcionamiento se requiere comentar la siguiente línea:
om.add_nln_constraint(mis_cons, [nb;nb], 1, fcn_mis, hess_mis, nodal_balance_vars)
en el archivo opf_setup.m de MATPOWER

En la carpeta 'Data_sets' se encuentran los siguientes conjuntos de prueba:
    -   sistema 500 nodos
    -   sistema 2742 nodos
    -   sistema de 10k nodos

La tesis asociada se encuentra en: https://repositorio.utp.edu.co/items/1ea95b25-d495-455d-8330-3c7aafdd1c22
    En este documento se encuentran los parámetros de control del algoritmo
Para preguntas contactarse con: caanacosta@utp.edu.co


