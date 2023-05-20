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
import math
assert cf
import sys
default_limit = 1000
sys.setrecursionlimit(default_limit*10)

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
    #TO DO: Inicializar las estructuras de datos
    data_structs = {'orderedData': None,'posiciones':None, 'grafoDir':None}
    
    data_structs['orderedData'] = mp.newMap(maptype='PROBING')
    data_structs['posiciones'] = mp.newMap(maptype='PROBING')
    data_structs['grafoDir'] = gr.newGraph(datastructure= "ADJ_LIST",directed=True)
    data_structs['lobos'] = lt.newList(datastructure='ARRAY_LIST')

    return data_structs

# Funciones para agregar informacion al modelo
def addWolfsData(data_structs, data):
    lt.addLast(data_structs['lobos'],data)
    entry = mp.get(data_structs['orderedData'],data['animal-id'])
    if entry is None:
        lst = lt.newList(datastructure='ARRAY_LIST')
        mp.put(data_structs['orderedData'],data['animal-id'],lst)

def add_data(data_structs, data):
    """
    Función para agregar información
    """
    #TO DO: ordena la información en un mapa (key= id lobo, value: lista de dicc con los eventos del lobo)
    entry = mp.get(data_structs['orderedData'],data['individual-local-identifier'])
    lstEvents = me.getValue(entry)
    lt.addLast(lstEvents,data)
        
def sortData(data_structs):
    keySet = mp.keySet(data_structs['orderedData'])
    for wolf in lt.iterator(keySet):
        entry = mp.get(data_structs['orderedData'],wolf)
        lstEvents = me.getValue(entry)
        quk.sort(lstEvents,sortCriteriaTimeStamp)

# Funciones para creacion de datos

def addTrackConnection(data_structs):
    """
    Adiciona las posiciones al grafo como vertices y arcos entre las
    estaciones adyacentes.

    Los vertices tienen por nombre la posición longitud-latitud seguido del lobo al que sirven
    """
    #TODO: Crear la función para estructurar los datos
    wolfIndividualEdges = 0
    mayorlat = 0
    menorlat = 1000
    mayorlon = -1000
    menorlon = 1000
    keyWolfSet = mp.keySet(data_structs['orderedData'])
    for wolf in lt.iterator(keyWolfSet):
        entry = mp.get(data_structs['orderedData'],wolf)
        wolfEventsLst = me.getValue(entry)
        lasttrack = None
        for track in lt.iterator(wolfEventsLst):
            if lasttrack != None:
                origin = formatVertex(data_structs,lasttrack)
                destination = formatVertex(data_structs,track)
                addPosition(data_structs, origin)
                addPosition(data_structs, destination)
                distance = getDistance(track,lasttrack)
                wolfIndividualEdges = addConnection(data_structs, origin, destination, distance,wolfIndividualEdges)
                
                if round(float(track['location-lat']),3) > mayorlat:
                    mayorlat = round(float(track['location-lat']),3)
                elif round(float(track['location-lat']),3) < menorlat:
                    menorlat = round(float(track['location-lat']),3)
                if round(float(track['location-long']),3) > mayorlon:
                    mayorlon = round(float(track['location-long']),3)
                elif round(float(track['location-long']),3) < menorlon:
                    menorlon = round(float(track['location-long']),3)
            lasttrack = track
    return wolfIndividualEdges, mayorlat,menorlat,mayorlon,menorlon
        
def formatVertex(data_structs,track):
    """
    Se formatea el nombre del vertice con la longitud y latitud del 
    event seguido del identificador del lobo.
    """
        
    latitud = str(round(float(track['location-lat']),3))
    nl = latitud.replace("-", "m")
    new_lat = nl.replace(".","p")
    
    longitud = str(round(float(track['location-long']),3))
    nlo = longitud.replace("-", "m")
    new_lon = nlo.replace(".","p")
    
    name = new_lon + "_" + new_lat + "_" + track['individual-local-identifier']
    position = new_lon + "_" + new_lat
    
    entry = mp.get(data_structs['posiciones'],position)
    if entry is None:
        lstWolfsEvents = lt.newList(datastructure='ARRAY_LIST')
        mp.put(data_structs['posiciones'],position,lstWolfsEvents)
        lt.addLast(lstWolfsEvents,track)
    else:
        lstWolfsEvents = me.getValue(entry)
        lstIds = lt.newList('ARRAY_LIST')
        for event in lt.iterator(lstWolfsEvents):
            """if track['individual-local-identifier'] != event['individual-local-identifier']:
                lt.addLast(lstWolfsEvents,track)"""
            lt.addLast(lstIds,event['individual-local-identifier'])
        containsWolfEvent = lt.isPresent(lstIds,track['individual-local-identifier'])
        if containsWolfEvent == 0:
            lt.addLast(lstWolfsEvents,track)
        
    return name

def addPosition(data_structs,vertexName):
    """
    Adiciona una posición como un vertice del grafo
    """
    if not gr.containsVertex(data_structs['grafoDir'], vertexName):
        gr.insertVertex(data_structs['grafoDir'], vertexName)
    return data_structs
    
def getDistance(origin, destination):
    lat2 = math.radians(abs(float(origin['location-lat'])))
    lon2 = math.radians(abs(float(origin['location-long'])))
    lat1 = math.radians(abs(float(destination['location-lat'])))
    lon1 = math.radians(abs(float(destination['location-long'])))
    sin2 = (math.sin((lat2-lat1)/2))**2
    part2 = math.cos(lat1)*math.cos(lat2)*(math.sin((lon2-lon1)/2))**2
    distancia = 2 * math.asin(math.sqrt(sin2 + part2)) * 6371
    return round(distancia,3)


def addConnection(data_structs, origin, destination, distance,wolfIndividualEdges):
    """
    Adiciona un arco entre dos posiciones
    """
    edge = gr.getEdge(data_structs['grafoDir'], origin, destination)
    if edge is None:
        gr.addEdge(data_structs['grafoDir'], origin, destination, distance)
        wolfIndividualEdges += 1
    return wolfIndividualEdges

def addPositionConnection(data_structs):
    
    """Por cada vertice (cada posicion) se recorre la lista
    de lobos servidas en dicha posicion y se crean
    arcos entre ellos para representar el cambio de ruta
    que se puede realizar en una posicion."""
    
    lstPosiciones = mp.keySet(data_structs['posiciones'])
    totalMTPs = 0
    totalWolfsMTPs = 0
    WeightZeroEdges = 0

    for key in lt.iterator(lstPosiciones):
        lstWolfsEvents = mp.get(data_structs['posiciones'], key)['value']
        if lt.size(lstWolfsEvents) >= 2:
            totalWolfsMTPs += lt.size(lstWolfsEvents)
            for event in lt.iterator(lstWolfsEvents):
                vertexName = key + '_' + event['individual-local-identifier']
                containsMTP = gr.containsVertex(data_structs['grafoDir'],key)
                if not containsMTP:
                    gr.insertVertex(data_structs['grafoDir'],key)
                    totalMTPs += 1
                gr.addEdge(data_structs['grafoDir'],key,vertexName,3000)
                WeightZeroEdges += 1

    vertexNum = gr.numVertices(data_structs['grafoDir'])
    return totalMTPs,totalWolfsMTPs,WeightZeroEdges,vertexNum
    
def firstFiveMPTs(data_structs):
    lstPosiciones = mp.keySet(data_structs['posiciones'])
    ret = []
    counter = 0
    for MTP in lt.iterator(lstPosiciones):
        if counter < 5:
            lstEvents = mp.get(data_structs['posiciones'],MTP)['value']
            if lt.size(lstEvents) >= 2:
                track = lt.firstElement(lstEvents)
                coordenadas = "("+ track['location-long'] + ", " + track['location-lat']+ ")"
                lst = [MTP,coordenadas,lt.size(lstEvents)]
                ret.append(lst)
                counter += 1
        else:
            return ret[0:5]
            
# Funciones de consulta

def graphSize(data_structs):
    """
    Retorna el tamaño de un grafo (numero total de vértices y nodos)
    """
    #TODO: Crear la función para obtener un dato de una lista
    totalVertices = gr.numVertices(data_structs['grafoDir'])
    totalEdges = gr.numEdges(data_structs['grafoDir'])
    
    return totalVertices,totalEdges


def data_size(lst):
    """
    Retorna el tamaño de la lista de datos
    """
    #TODO: Crear la función para obtener el tamaño de una lista
    return lt.size(lst)

def imprimir(control,nodoPrueba):
    lstNodos = gr.adjacents(control['grafoDir'],nodoPrueba)
    for nodo in lt.iterator(lstNodos):
        print(gr.getEdge(control['grafoDir'],nodoPrueba,nodo))
    
    
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


def sortCriteriaTimeStamp(data1, data2):
    """sortCriteria criterio de ordenamiento para las funciones de ordenamiento
    """
    #TODO: Crear función comparadora para ordenar
    if datetime.datetime.strptime(data1['timestamp'],'%Y-%m-%d %H:%M')<datetime.datetime.strptime(data2['timestamp'],'%Y-%m-%d %H:%M'):
        return True
    else:
        return False
        

def sortTimeStamp(lst):
    """
    Función encargada de ordenar la lista con los datos
    """
    #TODO: Crear función de ordenamiento
    return quk.sort(lst,sortCriteriaTimeStamp)

def SortLat(track,mayorLat,menorLat):
    if round(float(track['location-lat']),3) > mayorLat:
        mayorLat = round(float(track['location-lat']),3)
    elif round(float(track['location-lat']),3) < menorLat:
        menorLat = round(float(track['location-lat']),3)
    return mayorLat,menorLat
