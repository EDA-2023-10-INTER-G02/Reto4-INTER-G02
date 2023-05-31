"""
 * Copyright 2020, Departamento de sistemas y Computación,
 * Universidad de Los Andes
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
import model
import time
import csv
import tracemalloc
from DISClib.ADT import graph as gr
import datetime
import tracemalloc
"""
El controlador se encarga de mediar entre la vista y el modelo.
"""


def new_controller():
    """
    Crea una instancia del modelo
    """
    #TO DO: Llamar la función del modelo que crea las estructuras de datos
    control = model.new_data_structs()
    return control


# Funciones para la carga de datos

def load_data(control, filename):
    """
    Carga los datos del reto
    """
    # TODO: Realizar la carga de datos
    tracksNum = 0
    tracksfile = cf.data_dir + filename
    input_file = csv.DictReader(open(tracksfile, encoding="utf-8"), delimiter=",")
    
    wolfsFile = cf.data_dir + 'wolfs/BA-Grey-Wolf-individuals-utf8-large.csv'
    input_fileWolfs = csv.DictReader(open(wolfsFile, encoding="utf-8"), delimiter=",")
    timeI = get_time()
    for wolf in input_fileWolfs:
        model.addWolfsData(control,wolf)
     
    for track in input_file:
        model.add_data(control,track)
        tracksNum += 1
        
    model.sortData(control)
    wolfIndividualVertex,mayorlat,menorlat,mayorlon,menorlon = model.addTrackConnection(control)
    rtas = model.addPositionConnection(control)
    timeF = get_time()
    wolfsNum = model.data_size(control['lobos'])
    graphSize = model.graphSize(control)
    totalTime = delta_time(timeI,timeF)
    
    return wolfsNum, rtas,tracksNum,wolfIndividualVertex,graphSize,mayorlat,menorlat,mayorlon,menorlon,totalTime
    
def printLoadData(control):
    tabulate = model.TabulateCD(control)
    return tabulate

def imprimir_nodo_prueba(control):
    nodo = 'm111p369_57p033_32261A_32261'
    model.imprimir(control,nodo)
    
    
# Funciones de ordenamiento

def sort(control):
    """
    Ordena los datos del modelo
    """
    #TODO: Llamar la función del modelo para ordenar los datos
    pass


# Funciones de consulta sobre el catálogo

def get_data(control, id):
    """
    Retorna un dato por su ID.
    """
    #TODO: Llamar la función del modelo para obtener un dato
    pass


def req_1(control,pointI,pointF):
    """
    Retorna el resultado del requerimiento 1
    """
    # TODO: Modificar el requerimiento 1
    return model.req_1(control,pointI,pointF)


def req_2(control, initialStation, destStation):
    """
    Retorna el resultado del requerimiento 2
    """
    # TODO: Modificar el requerimiento 2
    return model.req_2(control, initialStation, destStation)
     


def req_3(control):
    """
    Retorna el resultado del requerimiento 3
    """
    # TODO: Modificar el requerimiento 3
    return model.req_3(control)
    


def req_4(control, ini, fin):
    """
    Retorna el resultado del requerimiento 4
    """
    # TODO: Modificar el requerimiento 4
    ini= ini.replace('(', '')
    ini= ini.replace(')', '')
    ini= ini.replace(',', '')
    ini= ini.split()
    fin= fin.replace('(', '')
    fin= fin.replace(')', '')
    fin= fin.replace(',', '')
    fin= fin.split()
    return model.req_4(control, ini, fin)



def req_5(control, identificador, maxDis, minPuntos):
    """
    Retorna el resultado del requerimiento 5
    """
    # TODO: Modificar el requerimiento 5
    return model.req_5(control, identificador, maxDis, minPuntos)

def req_6(control,animal_sex, ini, fin):
    """
    Retorna el resultado del requerimiento 6
    """
    # TODO: Modificar el requerimiento 6
    ini= datetime.datetime.strptime(ini,'%Y-%m-%d %H:%M')
    fin= datetime.datetime.strptime(fin,'%Y-%m-%d %H:%M')
    #timeI = get_time()
    tracemalloc.start()
    memory_i = get_memory()
    x= model.req_6(control,animal_sex, ini, fin)
    memory_f = get_memory()
    tracemalloc.stop()
    #timeF = get_time()
    #time = delta_time(timeI,timeF)
    final_mem = delta_memory(memory_f,memory_i)
    return x
    


def req_7(control,dateI, dateF,tempMax,tempMin):
    """
    Retorna el resultado del requerimiento 7
    """
    # TODO: Modificar el requerimiento 7
    dateI = '2012-11-28 00:00'
    dateF = '2014-05-17 23:59'
    tempMax = 9.7
    tempMin = -17.3
    dateI= datetime.datetime.strptime(dateI,'%Y-%m-%d %H:%M')
    dateF= datetime.datetime.strptime(dateF,'%Y-%m-%d %H:%M')
    return model.req_7(control,dateI, dateF,tempMax,tempMin)


def req_8(control):
    """
    Retorna el resultado del requerimiento 8
    """
    # TODO: Modificar el requerimiento 8
    pass


# Funciones para medir tiempos de ejecucion

def get_time():
    """
    devuelve el instante tiempo de procesamiento en milisegundos
    """
    return float(time.perf_counter()*1000)


def delta_time(start, end):
    """
    devuelve la diferencia entre tiempos de procesamiento muestreados
    """
    elapsed = float(end - start)
    return elapsed

def get_memory():
    """
    toma una muestra de la memoria alocada en instante de tiempo
    """
    return tracemalloc.take_snapshot()


def delta_memory(stop_memory, start_memory):
    """
    calcula la diferencia en memoria alocada del programa entre dos
    instantes de tiempo y devuelve el resultado en bytes (ej.: 2100.0 B)
    """
    memory_diff = stop_memory.compare_to(start_memory, "filename")
    delta_memory = 0.0

    # suma de las diferencias en uso de memoria
    for stat in memory_diff:
        delta_memory = delta_memory + stat.size_diff
    # de Byte -> kByte
    delta_memory = delta_memory/1024.0
    return delta_memory

def get_memory():
    """
    toma una muestra de la memoria alocada en instante de tiempo
    """
    return tracemalloc.take_snapshot()


def delta_memory(stop_memory, start_memory):
    """
    calcula la diferencia en memoria alocada del programa entre dos
    instantes de tiempo y devuelve el resultado en bytes (ej.: 2100.0 B)
    """
    memory_diff = stop_memory.compare_to(start_memory, "filename")
    delta_memory = 0.0

    # suma de las diferencias en uso de memoria
    for stat in memory_diff:
        delta_memory = delta_memory + stat.size_diff
    # de Byte -> kByte
    delta_memory = delta_memory/1024.0
    return delta_memory

