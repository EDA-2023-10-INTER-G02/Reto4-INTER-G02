"""
 * Copyright 2020, Departamento de sistemas y Computación, Universidad
 * de Los Andes
 *
 *
 * Desarrolado para el curso ISIS1225 - Estructuras de Datos y Algoritmos
 *
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along withthis program.  If not, see <http://www.gnu.org/licenses/>.
 """

import config as cf
import sys
import controller
from DISClib.ADT import list as lt
from DISClib.ADT import stack as st
from DISClib.ADT import queue as qu
from DISClib.ADT import map as mp
from DISClib.DataStructures import mapentry as me
assert cf
from tabulate import tabulate
import traceback

"""
La vista se encarga de la interacción con el usuario
Presenta el menu de opciones y por cada seleccion
se hace la solicitud al controlador para ejecutar la
operación solicitada
"""


def new_controller():
    """
        Se crea una instancia del controlador
    """
    #TO DO: Llamar la función del controlador donde se crean las estructuras de datos
    return controller.new_controller()


def print_menu():
    print("Bienvenido")
    print("1- Cargar información")
    print("2- Ejecutar Requerimiento 1")
    print("3- Ejecutar Requerimiento 2")
    print("4- Ejecutar Requerimiento 3")
    print("5- Ejecutar Requerimiento 4")
    print("6- Ejecutar Requerimiento 5")
    print("7- Ejecutar Requerimiento 6")
    print("8- Ejecutar Requerimiento 7")
    print("9- Ejecutar Requerimiento 8")
    print("0- Salir")

def opciones_tamaño():
    tamano = int(input("Elija el tamaño del archivo:\n1.Small (1%)\n2.5%\n3.10%\n4.20%\n5.30%\n6.50%\n7.80%\n8.Large (100%)\n"))
    tamanos=["small.csv","5pct.csv","10pct.csv","20pct.csv","30pct.csv","50pct.csv","80pct.csv","large.csv"]
    return tamanos[tamano-1]

def printLoadData(control):
    info = controller.printLoadData(control)
    print(tabulate(info,headers=["Location Long","Location Lat","node-id","individual-id","adjacents nodes"],tablefmt='grid'))
    
def load_data(control,filename):
    """
    Carga los datos
    """
    #TODO: Realizar la carga de datos
    info = controller.load_data(control,"wolfs/BA-Grey-Wolf-tracks-utf8-" +filename)
    return info


def print_Tabulate_req1(lst):
    """
        Función que imprime el tabulate del req 1
    """
    #TO DO: Realizar la función para imprimir un elemento
    lstOflst = []
    for vertex in lt.iterator(lst):
        lstOflst.append(vertex['elements'])
    print(tabulate(lstOflst,headers=["Identificador del punto de encuentro","Location Long","Location Lat","Número de lobos que transitan por el punto",
    "3 primeros y últimos identificadores de lobos","Siguiente nodo","Distancia al siguiente vértice"],tablefmt='grid',maxheadercolwidths=15,maxcolwidths=15))

def printTabulate_req3(lst):
    print(tabulate(lst,headers=['SCCID','Node IDs','SCC size','Min Lat','Max Lat','Min Lon','Max Lon','Wolf Count','Wolf Details'],tablefmt='grid',maxheadercolwidths=[2,28,3,5,5,5,5,5],
                   maxcolwidths=[2,28,3,5,5,5,5,5],numalign="right"))

def printTabulateWolfInfo(lst):
    print(tabulate([lst],headers=['individual-id','animal-taxon','animal-life-stage','animal-sex','study-site','Total-travel-dist','deployment-comments'],tablefmt='grid',maxheadercolwidths=20,maxcolwidths=20))

def printTabulatePathInfo(lst):
    print(tabulate(lst,headers=['node-id','location long','location lat','individual-id','wolf count'],tablefmt='grid',maxheadercolwidths=30,maxcolwidths=30))
    
    
def print_req_1(control):
    """
        Función que imprime la solución del Requerimiento 1 en consola
    """
    # TODO: Imprimir el resultado del requerimiento 1
    pointI = input('Ingrese el Identificador del punto de encuentro de origen: ')
    pointF = input('Ingrese el Identificador del punto de encuentro de destino: ')
    return controller.req_1(control,pointI,pointF)

def print_req_2(control):
    """
        Función que imprime la solución del Requerimiento 2 en consola
    """
    # TODO: Imprimir el resultado del requerimiento 2
    ini= input('Ingrese el Identificador del punto de encuentro de origen: ')
    fin= input('Ingrese el Identificador del punto de encuentro de destino: ')
    lista_camino2, num_gathering, num_vert, num_track, edges, tot_dist=controller.req_2(control, ini, fin)
    print('El total de nodos en el camino es ' + str(num_vert))
    print('El total de puntos de encuentro es ' + str(num_gathering))
    print( 'El total de puntos de seguimiento es '+ str(num_track) )
    print('La distancia total es de ' + str(tot_dist)+ (' km'))
    print('Cantidad nodos camino BFS ' + str(num_vert))
    print('Cantidad arcos camino BFS ' + str(edges))
    table= []
    headers= 'location-log-aprox','location-lat-aprox', 'node-id', 'individual-id', 'individual-count','edge-to', 'edge-distance-km'
    for coor in lt.iterator(lista_camino2):
        table.append([coor[0], coor[1], coor[2], coor[3],coor[4],coor[5], coor[6]])
    print(tabulate(table,headers=headers,tablefmt='grid',maxcolwidths= 11, maxheadercolwidths= 11))




def print_req_3(control):
    """
        Función que imprime la solución del Requerimiento 3 en consola
    """
    # TODO: Imprimir el resultado del requerimiento 3
    return controller.req_3(control)


def print_req_4(control):
    """
        Función que imprime la solución del Requerimiento 4 en consola
    """
    # TODO: Imprimir el resultado del requerimiento 4
    ini= input('Nodo inicial')
    fin= input('Nodo final')
    lista_camino2, lista_gathering, num_arc, num_vert, size_gath, dist_tot, dist_menor_ini, dist_menor_dest, tot_dist, ino, longino, latino,individual_idino, dest, longdest, latdest, individual_id_dest= controller.req_4(control, ini, fin)
    print('La distancia de la coordenada ingresada como punto de origen al punto de encuentro mas cercano es: ' + str(dist_menor_ini))
    table0_1= []
    table0_1.append([ino, longino, latino,individual_idino])
    headers01= 'node-id', 'location-log-aprox','location-lat-aprox', 'individual-id'
    print(tabulate(table0_1,headers=headers01,tablefmt='grid',maxcolwidths= 11, maxheadercolwidths= 11))
    print('La distancia de la coordenada ingresada como punto de destino al punto de encuentro mas cercano es: ' + str(dist_menor_dest))
    table0_2=[]
    table0_2.append([dest, longdest, latdest, individual_id_dest])
    print(tabulate(table0_2,headers=headers01,tablefmt='grid',maxcolwidths= 11, maxheadercolwidths= 11))
    print('La distancia del punto de encuentro de origen y el punto de encuentro de destino es: '+ str(dist_tot))
    print('La distancia total del recorrido es : '+ str(tot_dist))
    print('La cantidad de arcos en el camino son:'+ str(num_arc))
    print('La cantidad de arcos en el camino son:'+ str(num_vert))
    print('La cantidad de puntos de encuentro en el recorrido son: ' + str(size_gath))
    table= []
    headers= 'src-node-id','location-lat-src', 'location-long-src', 'tgt-node-id','location-lat-tgt', 'location-long-tgt','individual-id' 'distance-km'
    for coor in lt.iterator(lista_camino2):
        table.append([coor[0], coor[1], coor[2], coor[3],coor[4],coor[5], coor[6],coor[7]])
    print(tabulate(table,headers=headers,tablefmt='grid',maxcolwidths= 11, maxheadercolwidths= 11))
    table2=[]
    headers2= 'node-id','location-long-aprox', 'location-lat-aprox','individual-id','individual-count'
    for coor in lt.iterator(lista_gathering):
        lista=[]
        for cada in lt.iterator(coor[3]):
            lista.append(cada)
        table2.append([coor[0], coor[1], coor[2], lista,coor[4]])
    print(tabulate(table2,headers=headers2,tablefmt='grid',maxcolwidths= 11, maxheadercolwidths= 11))

def print_req_5(control, identificador, maxDis, minPuntos):
    """
        Función que imprime la solución del Requerimiento 5 en consola
    """
    # TODO: Imprimir el resultado del requerimiento 5
    #print(control)
    tabla = controller.req_5(control, identificador, maxDis, minPuntos)
    print("=============== Req No. 5 Inputs ===============")
    print("Min of Gathering/Tracking Points to visit:", minPuntos)
    print("Max roundtrip traveling distance:", maxDis, "[km]")
    print("Origin gathering point: '" + identificador + "'\n")
    if lt.isEmpty(tabla):
        print("The origin gathering point is not conected to other points.")
    else:
        print("The origin gathering point is conected to other points.")
        print("Preparing to find the longest path...\n")
        print("=============== Req No. 5 Answer ===============")
        print("The are", lt.size(tabla), "possible paths from point '" + identificador + "'")
        print("The minimun number of gathering points to visit is:", minPuntos)
        print("The maximun roundtrip distance to travel is:", maxDis, "[km]")
        print(tabulate([tabla['elements'][1]], headers='keys', tablefmt='grid',maxcolwidths= 20, maxheadercolwidths= 20))


def print_req_6(control):
    """
        Función que imprime la solución del Requerimiento 6 en consola
    """
    # TODO: Imprimir el resultado del requerimiento 6
    animal_sex= input('El sexo registrado del animal: ')
    ini= input('Fecha inicial: ')
    fin= input('Fecha final: ')
    wolfMaxDist, maxWolfInfo, pathInfoMax,wolfMinDist,minWolfInfo,pathInfoMin = controller.req_6(control, animal_sex, ini, fin)
    print('\n--------- Parte 1 ---------')
    print('El individuo con mayor recorrido fue: ' +wolfMaxDist)
    printTabulateWolfInfo(maxWolfInfo)
    totalNodes, totalEdges, nodesInfo = pathInfoMax
    print('El camino más largo para el lobo ' +wolfMaxDist+' tiene: \nTotal de nodos: '+str(totalNodes)+'\nTotal de arcos: '+str(totalEdges))
    print('Detalles del camino con los tres primeros y últimos nodos: ')
    printTabulatePathInfo(nodesInfo)
    print('\n--------- Parte 2 ---------')
    print('El individuo con menor recorrido fue: ' +wolfMinDist)
    printTabulateWolfInfo(minWolfInfo)
    totalNodes, totalEdges, nodesInfo = pathInfoMin
    print('El camino más largo para el lobo ' +wolfMinDist+' tiene: \nTotal de nodos: '+str(totalNodes)+'\nTotal de arcos: '+str(totalEdges))
    print('Detalles del camino con los tres primeros y últimos nodos: ')
    printTabulatePathInfo(nodesInfo)
    print('\n')


def print_req_7(control,dateI, dateF,tempMax,tempMin):
    """
        Función que imprime la solución del Requerimiento 7 en consola
    """
    # TODO: Imprimir el resultado del requerimiento 7
    numScc, FivemaxManInfo = controller.req_7(control,dateI, dateF,tempMax,tempMin)
    print('\nHay ' +str(numScc)+' componentes fuertemente conectados en el grafo creado')
    print('Las tres manadas con mayor dominio sobre el territorio son:\n')
    printTabulate_req3(FivemaxManInfo)


def print_req_8(control):
    """
        Función que imprime la solución del Requerimiento 8 en consola
    """
    # TODO: Imprimir el resultado del requerimiento 8
    pass


# Se crea el controlador asociado a la vista
control = new_controller()

# main del reto
if __name__ == "__main__":
    """
    Menu principal
    """
    working = True
    #ciclo del menu
    while working:
        print_menu()
        inputs = input('Seleccione una opción para continuar\n')
        try:
            if int(inputs) == 1:
                filename = opciones_tamaño()
                print("Cargando información de los archivos ....\n")
                wolfsNum, rtas,tracksNum,wolfIndividualVertex,graphSize,mayorlat,menorlat,mayorlon,menorlon,totalTime = load_data(control,filename)
                totalMTPs,totalWolfsMTPs,WeightZeroEdges,vertexNum = rtas
                totalVertices,totalEdges = graphSize
                print("Total de lobos reconocidos en el estudio: " +str(wolfsNum))
                print("Total de puntos de encuentro reconocidos (MTPs): " +str(totalMTPs))
                print("Total de lobos presentes en los puntos de encuentro (MTPs): " +str(totalWolfsMTPs))
                print("Total de eventos cargados durante el estudio: "+str(tracksNum))
                print('Total de vértices creados para representar el movimiento de los individuos: '+str(wolfIndividualVertex))
                print("Total de arcos para unir nodos de encuentro y puntos de seguimiento: "+str(WeightZeroEdges))
                print('\nTotal de vértices en el grafo: ' +str(totalVertices))
                print('Total de arcos en el grafo: '+str(totalEdges))
                print('Rango del área rectangular que ocupan los lobos grises de Boutin Alberta en Canadá:')
                print('Latitudes: desde '+str(menorlat)+' hasta '+str(mayorlat))
                print('Longitudes: desde '+str(menorlon)+' hasta '+str(mayorlon))
                print('-> Tiempo de ejecución: '+str(totalTime)+"\n"+"\n"+"Primeros y últimos 5 nodos cargados en el grafo dirigido")
                printLoadData(control)
                #controller.imprimir_nodo_prueba(control)
                print("\n")
                
            elif int(inputs) == 2:
                totalNodes,trackingPoints,totalDist,Mtps,rta = print_req_1(control)
                if rta == []:
                    print('No hay camino')
                else:
                    print('\nTotal de nodos del camino: '+str(totalNodes))
                    print("Total de puntos de encuentro: " + str(Mtps))
                    print('Número de tracking points: ' +str(trackingPoints))
                    print("Distancia total entre punto de origen y destino: "+str(round(totalDist,4))+str(' km')+"\n")
                    print("Cinco primeros y últimos vértices de la ruta: ")
                    print_Tabulate_req1(rta)
                    print("\n")

            elif int(inputs) == 3:
                print_req_2(control)

            elif int(inputs) == 4:
                numScc, FivemaxManInfo =print_req_3(control)
                print('\nHay ' +str(numScc)+' componentes fuertemente conectados en el grafo')
                print('Las cinco manadas con mayor dominio sobre el territorio son:\n')
                printTabulate_req3(FivemaxManInfo)

            elif int(inputs) == 5:
                print_req_4(control)

            elif int(inputs) == 6:
                identificador = input("Escriba el identificador del punto de encuentro de origen: ")
                maxDis = input("Digite la máxima distancia que el guardabosques puede recorrer en kilometros: ")
                minPuntos = input("Digite la minima cantidad de puntos de encuentros que el guardabosques desea inspeccionar: ")
                print_req_5(control, identificador, maxDis, minPuntos)

            elif int(inputs) == 7:
                print_req_6(control)

            elif int(inputs) == 8:
                dateI = input('Ingrese fecha inicial: ')
                dateF = input('Ingrese fecha final: ')
                tempMin = float(input('Ingrese temperatura mínima: '))
                tempMax = float(input('Ingrese temperatura máxima: '))
                print_req_7(control,dateI, dateF,tempMax,tempMin)

            elif int(inputs) == 9:
                print_req_8(control)

            elif int(inputs) == 0:
                working = False
                print("\nGracias por utilizar el programa")
                
            else:
                print("Opción errónea, vuelva a elegir.\n")
        except Exception as exp:
            print("ERR:", exp)
            traceback.print_exc()
    sys.exit(0)

