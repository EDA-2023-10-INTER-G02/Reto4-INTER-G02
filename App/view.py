﻿"""
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


def print_data(control, id):
    """
        Función que imprime un dato dado su ID
    """
    #TODO: Realizar la función para imprimir un elemento
    pass

def print_req_1(control):
    """
        Función que imprime la solución del Requerimiento 1 en consola
    """
    # TODO: Imprimir el resultado del requerimiento 1
    pointI = input('Ingrese el Identificador del punto de encuentro de origen: ')
    pointF = input('Ingrese el Identificador del punto de encuentro de destino: ')
    controller.req_1(control,pointI,pointF)


def print_req_2(control):
    """
        Función que imprime la solución del Requerimiento 2 en consola
    """
    # TODO: Imprimir el resultado del requerimiento 2
    ini= input('Nodo inicial')
    fin= input('Nodo final')
    path, longitud= controller.req_2(control, ini, fin)
    print(path)
    print(longitud)


def print_req_3(control):
    """
        Función que imprime la solución del Requerimiento 3 en consola
    """
    # TODO: Imprimir el resultado del requerimiento 3
    pass


def print_req_4(control):
    """
        Función que imprime la solución del Requerimiento 4 en consola
    """
    # TODO: Imprimir el resultado del requerimiento 4
    pass


def print_req_5(control):
    """
        Función que imprime la solución del Requerimiento 5 en consola
    """
    # TODO: Imprimir el resultado del requerimiento 5
    pass


def print_req_6(control):
    """
        Función que imprime la solución del Requerimiento 6 en consola
    """
    # TODO: Imprimir el resultado del requerimiento 6
    pass


def print_req_7(control):
    """
        Función que imprime la solución del Requerimiento 7 en consola
    """
    # TODO: Imprimir el resultado del requerimiento 7
    pass


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
                wolfsNum, rtas,tracksNum,wolfIndividualEdges,graphSize,mayorlat,menorlat,mayorlon,menorlon,totalTime = load_data(control,filename)
                totalMTPs,totalWolfsMTPs,WeightZeroEdges,vertexNum = rtas
                totalVertices,totalEdges = graphSize
                print("Total de lobos reconocidos en el estudio: " +str(wolfsNum))
                print("Total de puntos de encuentro reconocidos (MTPs): " +str(totalMTPs))
                print("Total de lobos presentes en los puntos de encuentro (MTPs): " +str(totalWolfsMTPs))
                print("Total de eventos cargados durante el estudio: "+str(tracksNum))
                print("Total de arcos para unir nodos de encuentro y puntos de seguimiento: "+str(WeightZeroEdges))
                print('Total de arcos creados para representar el movimiento de los individuos: '+str(wolfIndividualEdges))
                print('\nTotal de vértices en el grafo: ' +str(totalVertices))
                print('Total de arcos en el grafo: '+str(totalEdges))
                print('Rango del área rectangular que ocupan los lobos grises de Boutin Alberta en Canadá:')
                print('Latitudes: desde '+str(menorlat)+' hasta '+str(mayorlat))
                print('Longitudes: desde '+str(menorlon)+' hasta '+str(mayorlon))
                print('-> Tiempo de ejecución: '+str(totalTime)+"\n"+"\n"+"Primeros y últimos 5 nodos cargados en el grafo dirigido")
                printLoadData(control)
                controller.imprimir_nodo_prueba(control)
                print("\n")
                
            elif int(inputs) == 2:
                print_req_1(control)

            elif int(inputs) == 3:
                print_req_2(control)

            elif int(inputs) == 4:
                print_req_3(control)

            elif int(inputs) == 5:
                print_req_4(control)

            elif int(inputs) == 6:
                print_req_5(control)

            elif int(inputs) == 7:
                print_req_6(control)

            elif int(inputs) == 8:
                print_req_7(control)

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

