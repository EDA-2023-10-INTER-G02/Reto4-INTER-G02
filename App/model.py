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
 *
 * Contribuciones:
 *
 * Dario Correal - Version inicial
 """


import config as cf
from DISClib.ADT import list as lt
from DISClib.ADT import stack as st
from DISClib.ADT import queue as qu
from DISClib.ADT import map as mp
from DISClib.ADT import minpq as mpq
from DISClib.ADT import indexminpq as impq
from DISClib.ADT import orderedmap as om
from DISClib.DataStructures import mapentry as me
from DISClib.ADT import graph as gr
from DISClib.Algorithms.Graphs import scc
from DISClib.Algorithms.Graphs import dijsktra as djk
from DISClib.Algorithms.Graphs import bellmanford as bf
from DISClib.Algorithms.Graphs import bfs
from DISClib.Algorithms.Graphs import dfs
from DISClib.Algorithms.Graphs import prim
from DISClib.Algorithms.Sorting import shellsort as sa
from DISClib.Algorithms.Sorting import insertionsort as ins
from DISClib.Algorithms.Sorting import selectionsort as se
from DISClib.Algorithms.Sorting import mergesort as merg
from DISClib.Algorithms.Sorting import quicksort as quk
import datetime
assert cf

"""
Se define la estructura de un catálogo de videos. El catálogo tendrá
dos listas, una para los videos, otra para las categorias de los mismos.
"""

# Construccion de modelos


def new_data_structs():
    """
    Inicializa las estructuras de datos del modelo. Las crea de
    manera vacía para posteriormente almacenar la información.
    """
    #TODO: Inicializar las estructuras de datos
    
    data_structs = {'ordered_data': None,'posiciones':None, 'grafo_dir':None}
    
    data_structs['ordered_data'] = {}
    data_structs['posiciones'] = mp.newMap(maptype='PROBING')
    data_structs['grafoDir'] = gr.newGraph(datastructure= "ADJ_LIST",directed=True)

    return data_structs

# Funciones para agregar informacion al modelo

def add_ordered_data(data_structs, data):
    """
    Función para agregar nuevos elementos a la lista
    """
    #TO DO: Crear la función para agregar elementos a una lista
    wolf_id = data['individual-local-identifier']
    date = data['timestamp']
    date_format = datetime.datetime.strptime(date, "%Y-%m-%d %H:%M")
    data['timestamp'] = date_format
    
    if wolf_id not in data_structs['ordered_data']:
        data_structs['ordered_data'][wolf_id] = [data]
    else:
        data_structs['ordered_data'][wolf_id].append(data)

def order_time(data_structs):
    for lobo in data_structs['ordered_data']:
        data_structs['ordered_data'][lobo].sort(key=lambda l: l['timestamp'])
    
    lst = []
    for lobo in data_structs['ordered_data']:
        for data in data_structs['ordered_data'][lobo]:
            lst.append(data)
            
    data_structs['ordered_data'] = lst
# Funciones para creacion de datos

def addTrackConnection(data_strcuts, track,lasttrack):
    origin = formatVertex(lasttrack)
    destination = formatVertex(track)
    addPosition(data_strcuts, origin)
    addPosition(data_strcuts, destination)


def addPosition(data_structs,positionId):
    """
    Adiciona una posición como un vertice del grafo
    """
    if not gr.containsVertex(data_structs['grafoDir'], positionId):
        gr.insertVertex(data_structs['grafoDir'], positionId)
    

def new_data(id, info):
    """
    Crea una nueva estructura para modelar los datos
    """
    #TODO: Crear la función para estructurar los datos
    pass


# Funciones de consulta

def get_data(data_structs, id):
    """
    Retorna un dato a partir de su ID
    """
    #TODO: Crear la función para obtener un dato de una lista
    pass


def data_size(data_structs):
    """
    Retorna el tamaño de la lista de datos
    """
    #TODO: Crear la función para obtener el tamaño de una lista
    pass


def req_1(data_structs):
    """
    Función que soluciona el requerimiento 1
    """
    # TODO: Realizar el requerimiento 1
    pass


def req_2(data_structs):
    """
    Función que soluciona el requerimiento 2
    """
    # TODO: Realizar el requerimiento 2
    pass


def req_3(data_structs):
    """
    Función que soluciona el requerimiento 3
    """
    # TODO: Realizar el requerimiento 3
    pass


def req_4(data_structs):
    """
    Función que soluciona el requerimiento 4
    """
    # TODO: Realizar el requerimiento 4
    pass


def req_5(data_structs):
    """
    Función que soluciona el requerimiento 5
    """
    # TODO: Realizar el requerimiento 5
    pass


def req_6(data_structs):
    """
    Función que soluciona el requerimiento 6
    """
    # TODO: Realizar el requerimiento 6
    pass


def req_7(data_structs):
    """
    Función que soluciona el requerimiento 7
    """
    # TODO: Realizar el requerimiento 7
    pass


def req_8(data_structs):
    """
    Función que soluciona el requerimiento 8
    """
    # TODO: Realizar el requerimiento 8
    pass


# Funciones utilizadas para comparar elementos dentro de una lista

def compare(data_1, data_2):
    """
    Función encargada de comparar dos datos
    """
    #TODO: Crear función comparadora de la lista
    pass

# Funciones de ordenamiento


def sort_criteria(data_1, data_2):
    """sortCriteria criterio de ordenamiento para las funciones de ordenamiento

    Args:
        data1 (_type_): _description_
        data2 (_type_): _description_

    Returns:
        _type_: _description_
    """
    #TODO: Crear función comparadora para ordenar
    pass


def sort(data_structs):
    """
    Función encargada de ordenar la lista con los datos
    """
    #TODO: Crear función de ordenamiento
    pass

# Funciones helper

def formatVertex(track):
    """
    Se formatea el nombrer del vertice con la longitud y latitud del event seguido del identificador del lobo.
    """
    latitud = str(round(float(track['location-lat']),4))
    nl = latitud.replace("-", "m")
    new_lat = nl.replace(".","p")
    
    longitud = str(round(float(track['location-lon']),4))
    nlo = longitud.replace("-", "m")
    new_lon = nlo.replace(".","p")
    
    name = new_lat + "_" + new_lon + "_" + track['individual-local-identifier']
    return name
    
position_id = "111p4534_56p1245_332"
print(position_id[0:4])
