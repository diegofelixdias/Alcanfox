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
print("\n")
print("BEM VINDO AO SSK MARK III")
print("\n")
print("\n")
#ESSE PROGRAMA FOI DESENVOLVIDO NO LABORATORIO DE RAIOS-X DA UNIVERSIDADE FEDERAL DO CEARA
#VOCE TEM TOTAL LIBERDADE PARA MODIFICAR ESSE PROGRAMA DO JEITO QUE BEM ENTENDER, MAS POR FAVOR CITAR
#PEGANDO OS VALORES DA MULTIPLICIDADE ////////////////////////////////////// NUMERO DE ATOMOS DA CELA UNITARIA
cristal_structure.show_summary().show_scatterers()
print(len(cristal_structure.scatterers())) #NUMERO DE ATOMOS NA CELA UNITARIA
print("QUAL A RADIACAO QUE ESTA SENDO UTILIZADA")
print("\n")
print("DIGITE 1 PARA Co")
print("DIGITE 2 PARA Cu")
print("DIGITE 3 PARA Mo")
print("DIGITE 4 PARA Fe")
print("DIGITE 5 PARA W")
print("DIGITE 6 PARA DIGITAR UM COMPRIMENTO DE ONDA")
print("\n")
numero = int(input())
if numero == 1:
    print("VOCE ESCOLHEU COBALTO")
    comprimento_de_onda = 1.7890*10**(-8)

if numero == 2:
    print("VOCE ESCOLHEU COBRE")
    comprimento_de_onda = 1.5405*10**(-8)

if numero == 3:
    print("VOCE ESCOLHEU MOLIBDENIO")
    comprimento_de_onda = 0.7098*10**(-8)

if numero == 4:
    print("VOCE ESCOLHEU FERRO")
    comprimento_de_onda = 1.936*10**(-8)

if numero == 5:
    print("VOCE ESCOLHEU TUNGSTENIO")
    comprimento_de_onda = 0.2091*10**(-8)

if numero == 6:
    print("VOCE VAI DIGITAR O COMPRIMENTO DE ONDA")
    print("POR FAVOR DIGITE O COMPRIMENTO DE ONDA EM ANGSTRONS")
    comprimento_de_onda = float(input())
    comprimento_de_onda = comprimento_de_onda*10**(-8)


raio_do_eletron = 2.818*10**(-13) #raio classico do eletron
pi = 3.14159265359 #pi
h = 4.135*10**(-15) #constante de plank em eV.s
velocidade_da_luz = 2.99792*10**(8) #velocidade da luz em m/s
energia2 = h*velocidade_da_luz*10**(2)/comprimento_de_onda #aqui vai ser a energia utilizada
energia = "8046.000" #teste com energia fixa
print("QUAL O RANGE DA MEDIDA, DIGITE EM GRAUS")
angulo = float(input())
d_hkl = comprimento_de_onda*10**(8)/(2*math.sin(angulo*(pi/180))) #usado para gerar os planos da difracao

#print cristal_structure.sites_cart()

#print B.multiplicity()
#print B.electron_count()
#print B.weight()
#a = int(B.multiplicity())


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
print("\n")
#////////////////////////////////////////////////////////////////////////////////
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
#CORRECOES DO FATOR DE ESPALHAMENTO f_' e f_'' //////////////////////////////////////////////////////
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


          y = y + 1
  i = i +1
  atomo.close()
  xy = xy + 1

#////////////////////////////////////////////////////////////////////////////////
#CALCULO DO FATOR DE ESPALHAMENTO f_0////////////////////////////////////////////////////
tabela = open("/home/diego/PycharmProjects/untitled2/table_WaasKirf.dat", "r")
linha = tabela.readlines()
i_2 = int(0)
i_3 = int(0)
b_2 = list(range(xy))
b_3 = list(range(11))
y_2 = int(0)
print("\n")
for i_2 in range(xy):
    numero = (input("qual o numero atomico do elemento? "))
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
print("\n")
print("TABELA COM TODOS OS VALORES DOS COEFICIENTES DE WAASMAIER KIRFEL PARA O CALCULO DO FATOR DE ESPALHAMENTO")
print("a1, a2, a3, a4, a5, c, b1, b2, b3, b4, b5")
print(b_2) #TABELA COM TODOS OS VALORES DOS COEFICIENTES DE WAASMAIER KIRFEL
print("\n")
print("ATOMOS DA CELA UNITARIA")
print(atomos)
print("\n")
print("VALOR DE F' PARA A CORREÇÃO DO FATOR DE ESPALHAMENTO")
print(k) #f'
print("\n")
print("VALOR DE F'' PARA A CORREÇÃO DO FATOR DE ESPALHAMENTO")
print(c) #f''
print("\n")
print("MULTIPLICIDADE DOS ATOMOS")
print(multiplicidades)
print("\n")
print("COORDENADAS DOS ATOMOS NA CELA UNITARIA")
print(coordenadas)
print("\n")
#////////////////////////////////////////////////////////////////////////////////

#////FATORES DE ESTRUTURA CORRIGIDOS???
cristal_structure.scattering_type_registry(table="it1992")
f_calc = cristal_structure.structure_factors(d_min=d_hkl).f_calc()
f_calc_sq = f_calc.as_intensity_array()
#f_calc_sq.show_summary().show_array()

#////////////////////////////////////////////////////////////////////////////////

#INDICES DE MILLER ///////////////////////////////////////////////////////////////////
print("NUMERO DE PLANOS DA MEDIDA")
print(f_calc_sq.size()) #DIZ QUANTOS PLANOS EXISTEM
print("PLANOS DA MEDIDA")
print list(f_calc_sq.indices()) #DIZ QUAIS SaO PS PLANOS
'''
print("TESTE INDICES DE MILLER ESPECIFICOS")
print(list(f_calc_sq.indices()[1])[1])
'''
print("\n")
f_calc_sq.as_cif_simple(
  array_type="calc", data_name="lab6", out=open("/home/diego/PycharmProjects/untitled2/lab6.hkl", "wb"))
#f_calc_sq.show_array()
#////////////////////////////////////////////////////////////////////////////////


from iotbx.cif import model

cif = model.cif()
cif_block = model.block()

unit_cell = cristal_structure.unit_cell()

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

print("ESSE E O VALOR DE 1/d PARA O CALCULO DO FATOR DE ESPALHAMENTO")
print(d_2)
print("\n")
#////////////////////////////////////////////////////////////////////////////////

#CALCULO PROPRIAMENTE DITO DO FATOR DE ESPALHAMENTO//////////////////////////////

i_4 = int(0)
i_5 = int(0)
i_6 = int(0)
f_0_tudo = list(range(xy))
i_7 = int(0)

for i_4 in range(xy):
    f_0 = list(range(f_calc_sq.size()))
    for i_6 in range(f_calc_sq.size()):
        for i_5 in range(5):
            f_0[i_6] = f_0[i_6] + ((b_2[i_4][i_5])*cmath.exp(-b_2[i_4][i_5+6]*(d_2[i_6]/2)**(2)))
        f_0[i_6] = f_0[i_6] + b_2[i_4][5]

    f_0_tudo[i_4] = f_0


print("\n")
print("FATOR DE ESPALHAMENTO SEM AS CORREÇÕES, CALCULADOS A PARTIR DE 5 TERMOS A PARTIR DA TABELA DE WASSKIRF PARA CADA UM DOS PLANOS")
print(f_0_tudo) #FATORES DE ESPALHAMENTO PARA TODOS OS PLANOS E TODOS OS ATOMOS
#////////////////////////////////////////////////////////////////////////////////

#FATOR DE ESPALHAMENTO SOMADO COM A CORRECAO///////////////////////////////
i_4 = int(0)
i_5 = int(0)
i_6 = int(0)
for i_4 in range(xy):
    for i_6 in range(f_calc_sq.size()):
        f_0_tudo[i_4][i_6] = f_0_tudo[i_4][i_6] + k[i_4] + 1j*c[i_4]

print("\n")
print("FATOR DE ESPALHAMENTO CORRIGIDO COM f' & f'' PARA TODOS OS ATOMOS PARA CADA UM DOS PLANOS")
print(f_0_tudo)
#////////////////////////////////////////////////////////////////////////////////

#CALCULO DO FATOR DE ABSORCAO /////////////////////////////////////////////////
#CALCULO DO GAMMA/////////////////////////////////////////////////////////////
print("\n")
print("GAMMA")
volume_da_cela_unitaria = float(cif_block["_cell_volume"])
gamma = (raio_do_eletron*(comprimento_de_onda)**2)/(pi*volume_da_cela_unitaria*10**(-24))
print(gamma)
print("\n")
#////////////////////////////////////////////////////////////////////////////////

#CALCULO DO F_0_2linhas/////////////////////////////////////////////////////////////////
for xy in (range(len(cristal_structure.scatterers()))):
    F_0_2 = F_0_2 + multiplicidades[xy] * c[xy]

print("FATOR F_0_2")
print(F_0_2)
print("\n")
#////////////////////////////////////////////////////////////////////////////////

#CALCULO DO FATOR DE ABSORCAO////////////////////////////////////////////////////
fator_de_absorcao = ((2*pi)/(comprimento_de_onda))*gamma * F_0_2
print("FATOR DE ABSORCAO")
print(fator_de_absorcao)
print("\n")
#////////////////////////////////////////////////////////////////////////////////

#POSICOES EQUIVALENTES PARA O CALCULO DO FATOR DE ESTRUTURA////////////////////////////
i_8 = int(0)
i_9 = int(0)
tabelaa = (range(xy + 1))
tabelab = (range(xy + 1))
for i_8 in range(xy + 1):
    for i_9 in range(1):
        tabelab[i_8] = int(0)
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
        tabelab[i_8] = []
        for op in ops:
            tabelaa = list(op * coordenadas[i_8])
            tabelab[i_8] = tabelab[i_8] + tabelaa
'''
print("POSICOES EQUIVALENTES")
print(tabelab) #ja e float
print("\n")
'''
#//////////////////////////////////////////////////////////////////

#CALCULO DO FATOR DE ESTRUTURA CORRIGIDO/////////////////////////////////////////
f_de_estrutura = range(f_calc_sq.size())
print("FATOR DE ESTRUTURA")
i = int(0)
y = int(0)
z = int(0)
f_espalhamento = [1, 1]
for i in range(f_calc_sq.size()):
    f_de_estrutura[i] = 0 + 0j
    for y in range(xy+1):
        for z in range(multiplicidades[y]):
            atomo_a = z + z*2
            f_de_estrutura[i] = f_de_estrutura[i] + (f_0_tudo[y][i]*cmath.exp((2*pi*1j)*(tabelab[y][atomo_a]*list(f_calc_sq.indices()[i])[0] + tabelab[y][atomo_a + 1]*list(f_calc_sq.indices()[i])[1] + tabelab[y][atomo_a + 2]*list(f_calc_sq.indices()[i])[2])))



print(f_de_estrutura)
print("\n")
#//////////////////////////////////////////////////////////////////////////////////////

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