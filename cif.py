# coding=utf-8
import iotbx.cif
from iotbx.cif import *
import math
import cmath
from cctbx import *
arq = open("/home/diego/PycharmProjects/untitled2/lab6.cif", "r")
cristal_as_cif = arq.read()
#import iotbx.cif
cristal_structure = iotbx.cif.reader(
  input_string=cristal_as_cif).build_crystal_structures()["lab6"]

#PEGANDO OS VALORES DA MULTIPLICIDADE ////////////////////////////////////// NUMERO DE ATOMOS DA CELA UNITARIA
cristal_structure.show_summary().show_scatterers()
print(len(cristal_structure.scatterers())) #NUMERO DE ATOMOS NA CELA UNITARIA
print("QUAL A RADIAcaO QUE ESTA SENDO UTILIZADA")
print("\n")
print("DIGITE 1 PARA Co")
print("DIGITE 2 PARA Cu")
print("DIGITE 3 PARA Mo")
print("DIGITE 4 PARA Fe")
print("DIGITE 5 PARA W")
#comprimento_de_onda = str(input("QUAL A RADIAcaO QUE ESTA SENDO UTILIZADA"))
#print cristal_structure.sites_cart()

#print B.multiplicity()
#print B.electron_count()
#print B.weight()
#a = int(B.multiplicity())
raio_do_eletron = 2.818*10**(-13) #raio classico do eletron

pi = 3.14159265359 #pi
h = 4.135*10**(-15) #constante de plank em eV.s
velocidade_da_luz = 2.99792*10**(8) #velocidade da luz em m/s
#comprimento_de_onda = int(input("qual o comprimento de onda utilizado")) #comprimento de onda da medida
comprimento_de_onda = 1.5405*10**(-8)
lamb = 1.5405
theta_maximo = 120
d_hkl = 0.889408
#energia = h*velocidade_da_luz/comprimento_de_onda #aqui vai ser a energia utilizada
energia = "8046.000" #teste com energia fixa
i = int(0)
y = int(0)
b = int(0)
c = range(len(cristal_structure.scatterers()))
k = range(len(cristal_structure.scatterers()))
xy = int(0)
atomos = range(len(cristal_structure.scatterers()))
atomos_variaveis = range(len(cristal_structure.scatterers()))
multiplicidades = range(len(cristal_structure.scatterers()))
coordenadas = range(len(cristal_structure.scatterers()))
F_0_2 = int(0)

print(atomos)
#PEGANDO ATOMOS DA CELA UNITARIA ///////////////////////////////////////
for scatterer in cristal_structure.scatterers():
  print "%s:" % scatterer.label, "%8.4f %8.4f %8.4f" % scatterer.site
  site_symmetry = cristal_structure.site_symmetry(scatterer.site)
  print "  point group type:", site_symmetry.point_group_type()
  print "  special position operator:", site_symmetry.special_op_simplified()
  coordenadas[xy] = scatterer.site #guardando as coordenadas em uma tabela
  atomos[xy] = scatterer.label #guardando atomos das cela unitaria
  #PEGANDO OS FATORES DE ESPALHAMENTO CORRIGIDOS ///////////////////////////////
  atomo = open("/home/diego/PycharmProjects/untitled2/fator de espalhamento/CorFEsp/%s.txt" %atomos[xy],"r")
  atomos_variaveis[xy] = cristal_structure.scatterers()[xy]
  multiplicidades[xy] = atomos_variaveis[xy].multiplicity()
  multiplicidades[xy] = int(multiplicidades[xy])
  print(multiplicidades[xy])
  print(xy)
  #CALCULO DAS CORRECOES DO FATOR DE ESPALHAMENTO f_1 e f_2 //////////////////////////////////////////////////////
  linhas = atomo.readlines()
  for i in range(len(linhas)):
      linha = linhas[i].strip()
      a = linha.split(" ")
      b = a[0]
      if (b == energia):
          c[xy] = float(a[-1])
          del(a[-1])
          del(a[0])

          for y in range(len(a)):
              if (a[y] != ""):
                  k[xy] = float(a[y])
                  print(k[xy])

          y = y + 1
  i = i +1
  atomo.close()
  xy = xy + 1
#CALCULO DO FATOR DE ESPALHAMENTO f_0  ///////////////////////////////////////////////
tabela = open("/home/diego/PycharmProjects/untitled2/table_WaasKirf.dat", "r")
linha = tabela.readlines()
i_2 = int(0)
i_3 = int(0)
b_2 = list(range(xy))
b_3 = list(range(11))
y_2 = int(0)
print((b_2[0]))
for i_2 in range(xy):
    numero = (input("qual o numero? "))
    elemento = (input("qual o elemento? "))
    for i_3 in range(len(linha)):
        if (linha[i_3] == "#S  %s  %s\n" %(numero,elemento)):
            aa = linha[i_3+3].split(" ")
            b_3 = list(range(11))
            x = int(0)
            y_2 = int(0)

            for x in range(len(aa)):
                if(aa[x] != "" and aa[x] != "\n"):
                    b_3[y_2] = float(aa[x])
                    y_2 = y_2 + 1

                b_2[i_2] = b_3


tabela.close()

print(b_2) #TABELA COM TODOS OS VALORES DOS COEFICIENTES DE WAASMAIER KIRFEL
print("\n")
print("atomos da cela unitaria f_1 e f_2 multiplicidades coodenadas")
print(atomos)
print(k) #f'
print(c) #f''
print(multiplicidades)
print(coordenadas)
print("\n")
#////FATORES DE ESTRUTURA CORRIGIDOS???
cristal_structure.scattering_type_registry(table="it1992")
f_calc = cristal_structure.structure_factors(d_min=d_hkl).f_calc()
f_calc_sq = f_calc.as_intensity_array()
#f_calc_sq.show_summary().show_array()
#INDICES DE MILLER ///////////////////////////////////////////////////////////////////
print(f_calc_sq.size()) #DIZ QUANTOS PLANOS EXISTEM
print list(f_calc_sq.indices()) #DIZ QUAIS SaO PS PLANOS
'''
print("TESTE INDICES DE MILLER ESPECIFICOS")
print(list(f_calc_sq.indices()[1])[1])
'''
print("\n")
f_calc_sq.as_cif_simple(
  array_type="calc", data_name="lab6", out=open("/home/diego/PycharmProjects/untitled2/lab6.hkl", "wb"))
#f_calc_sq.show_array()



from iotbx.cif import model

cif = model.cif()
cif_block = model.block()

unit_cell = cristal_structure.unit_cell()

print type(unit_cell)

params = unit_cell.parameters()
cif_block["_cell_length_a"] = params[0]
cif_block["_cell_length_b"] = params[1]
cif_block["_cell_length_c"] = params[2]
cif_block["_cell_angle_alpha"] = params[3]
cif_block["_cell_angle_beta"] = params[4]
cif_block["_cell_angle_gamma"] = params[5]
cif_block["_cell_volume"] = unit_cell.volume()
cell_a = float(cif_block["_cell_length_a"])
cell_b = float(cif_block["_cell_length_b"])
cell_c = float(cif_block["_cell_length_c"])
cell_alpha = float(cif_block["_cell_angle_alpha"])
cell_beta = float(cif_block["_cell_angle_beta"])
cell_gamma = float(cif_block["_cell_angle_gamma"])

#CALCULO DO FATOR DE ESPALHAMENTO /////////////////////////////////////////////
#CALCULO DO d (distancia interplanar) ////////////////////////////////////////
hkl = int(0)
d_1 = list(range(f_calc_sq.size()))
d_2 = list(range(f_calc_sq.size()))
for hkl in range(f_calc_sq.size()):
    d_1[hkl] = ((1)/(1 + 2*math.cos(cell_alpha*((pi)/180))*math.cos(cell_beta*((pi)/180))*math.cos(cell_gamma*((pi)/180)) - (math.cos(cell_alpha*((pi)/180)))**(2) - (math.cos(cell_beta*((pi)/180)))**(2) - (math.cos(cell_gamma*((pi)/180)))**(2)))*((((list(f_calc_sq.indices()[hkl])[0])**(2)*(math.sin(cell_alpha*((pi)/180)))**(2))/(cell_a**2)) + (((list(f_calc_sq.indices()[hkl])[1])**(2)*(math.sin(cell_beta*((pi)/180)))**(2))/(cell_b**2)) + (((list(f_calc_sq.indices()[hkl])[2])**(2)*(math.sin(cell_gamma*((pi)/180)))**(2))/(cell_c**2)) + ((2*(list(f_calc_sq.indices()[hkl])[0])*(list(f_calc_sq.indices()[hkl])[1]))/(cell_a*cell_b))*(math.cos(cell_alpha*((pi)/180))*math.cos(cell_beta*((pi)/180)) - math.cos(cell_gamma*((pi)/180))) + ((2*(list(f_calc_sq.indices()[hkl])[1])*(list(f_calc_sq.indices()[hkl])[2]))/(cell_b*cell_c))*(math.cos(cell_beta*((pi)/180))*math.cos(cell_gamma*((pi)/180)) - math.cos(cell_alpha*((pi)/180))) + ((2*(list(f_calc_sq.indices()[hkl])[2])*(list(f_calc_sq.indices()[hkl])[0]))/(cell_a*cell_c))*(math.cos(cell_gamma*((pi)/180))*math.cos(cell_alpha*((pi)/180)) - math.cos(cell_beta*((pi)/180))))
    d_2[hkl] = math.sqrt(d_1[hkl]) #ISSO AQUI JA E O 1/dhkl

print("ESSE E O VALOR DE 1/d")
print(d_2)
print("\n")
#CALCULO PROPRIAMENTE DITO DO FATOR DE ESPALHAMENTO

i_4 = int(0)
i_5 = int(0)
i_6 = int(0)
f_0_tudo = list(range(xy))
i_7 = int(0)

print(xy)
for i_4 in range(xy):
    f_0 = list(range(f_calc_sq.size()))
    for i_6 in range(f_calc_sq.size()):
        for i_5 in range(5):
            print(b_2[i_4][i_5])
            print(b_2[i_4][i_5+6])
            f_0[i_6] = f_0[i_6] + ((b_2[i_4][i_5])*cmath.exp(-b_2[i_4][i_5+6]*(d_2[i_6]/2)**(2)))
            print(f_0[i_6])
        f_0[i_6] = f_0[i_6] + b_2[i_4][5]
        print(f_0[i_6])

    print(f_0)
    print("\n")
    print(i_4)
    f_0_tudo[i_4] = f_0
    print("\n")
    print(f_0_tudo)

print("\n")
print(f_0_tudo) #FATORES DE ESPALHAMENTO PARA TODOS OS PLANOS E TODOS OS ATOMOS


#FATOR DE ESPALHAMENTO SOMADO COM A CORRECAO///////////////////////////////
i_4 = int(0)
i_5 = int(0)
i_6 = int(0)
for i_4 in range(xy):
    for i_6 in range(f_calc_sq.size()):
        print("\n")
        print(k[i_4])
        print(c[i_4])
        f_0_tudo[i_4][i_6] = f_0_tudo[i_4][i_6] + k[i_4] + 1j*c[i_4]
        print("\n")

print(f_0_tudo)

#CALCULO DO FATOR DE ABSORCAO /////////////////////////////////////////////////
#CALCULO DO GAMMA/////////////////////////////////////////////////////////////
print("GAMMA")
volume_da_cela_unitaria = float(cif_block["_cell_volume"])
gamma = (raio_do_eletron*(comprimento_de_onda)**2)/(pi*volume_da_cela_unitaria*10**(-24))
print(gamma)
print("\n")
#CALCULO DO F_0_2/////////////////////////////////////////////////////////////////
for xy in (range(len(cristal_structure.scatterers()))):
    F_0_2 = F_0_2 + multiplicidades[xy] * c[xy]

print("FATOR F_0_2")
print(F_0_2)
print("\n")

#CALCULO DO FATOR DE ABSORCAO////////////////////////////////////////////////////
fator_de_absorcao = ((2*pi)/(comprimento_de_onda))*gamma * F_0_2
print("FATOR DE ABSORCAO")
print(fator_de_absorcao)
print("\n")

#CALCULO DO FATOR DE ESTRUTURA CORRIGIDO
f_de_estrutura = range(f_calc_sq.size())
print("FATOR DE ESTRUTURA")
i = int(0)
y = int(0)
f_espalhamento = [1, 1]
print(len(f_de_estrutura))
print("\n")
for i in range(f_calc_sq.size()):
    f_de_estrutura[i] = 0 + 0j
    for y in range(xy+1):
        f_de_estrutura[i] = f_de_estrutura[i] + (multiplicidades[y]*f_0_tudo[y][i]*cmath.exp((2*pi*1j)*(coordenadas[y][0]*list(f_calc_sq.indices()[i])[0] + coordenadas[y][1]*list(f_calc_sq.indices()[i])[1] + coordenadas[y][2]*list(f_calc_sq.indices()[i])[2])))

        print(multiplicidades[y])
        print(f_0_tudo[y][i])
        print("\n")
        print(f_de_estrutura[i])


print(f_de_estrutura)
print("\n")
i_8 = int(0)
i_9 = int(0)
tabelaa = list(range(xy + 1))
print(tabelaa)
for i_8 in range(xy + 1):
    for i_9 in range(1):
        print(i_8)
        # Esse objeto contém as informações sobre o grupo espacial 221
        grupo_espacial_info = sgtbx.space_group_info(221)
        # a partir dele, obtemos a representação em python do grupo em si, que permite fazer operações
        grupo_espacial = grupo_espacial_info.group()
        # geramos uma representação genérica das tabelas wyckoff
        tabela_wyckoff = grupo_espacial_info.wyckoff_table()
        # obtemos a simetria do atomo 2 a partir das informações do grupo espacial
        simetria = cristal_structure.site_symmetry(coordenadas[i_8])
        # geramos agora um mapeamento das matrizes de transformação para as posições, a partir da simetria
        mapeamento = tabela_wyckoff.mapping(simetria)
        # obtemos a posição a partir do mapeamento
        posicao_wyckoff = mapeamento.position()
        # a partir da posição, obtemos as matrizes de transformação
        ops = posicao_wyckoff.unique_ops(grupo_espacial)
        # aplicamos cada matriz ao ponto original e obtemos, assim, os pontos seguintes
        for op in ops:
            print(op * coordenadas[i_8])
            tabelab = op * coordenadas[i_8]
            print(tabelab)
            tabelac = list(tabelab)
            print(type(tabelac))
            #tabelaa[i_8] = tabelaa[i_8] + tabelac



print(tabelaa)


'''

space_group = cristal_structure.space_group()
symop_loop = model.loop(header=("_space_group_symop_id",
                                "_space_group_symop_operation_xyz"))
for symop_id, symop in enumerate(space_group):
  symop_loop.add_row((symop_id + 1, symop.as_xyz()))

space_group_type = cristal_structure.space_group_info().type()
cif_block["_space_group_crystal_system"] = space_group.crystal_system().lower()
cif_block["_space_group_IT_number"] = space_group_type.number()
cif_block["_space_group_name_H-M_alt"] = space_group_type.lookup_symbol()
cif_block["_space_group_name_Hall"] = space_group_type.hall_symbol()
cif_block.add_loop(symop_loop)

cif["lab6"] = cif_block
print cif
'''


arq.close()