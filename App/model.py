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
from math import sin, cos, sqrt, atan2, radians
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
    data_structs['grafoNoDir'] = gr.newGraph(datastructure= "ADJ_LIST",directed=False)
    data_structs['lobos'] = lt.newList(datastructure='ARRAY_LIST')
    data_structs['MTPs'] = mp.newMap(maptype='PROBING')
    data_structs['individualPoints'] = mp.newMap(maptype='PROBING')
    

    return data_structs

# Funciones para agregar informacion al modelo
def addWolfsData(data_structs, data):
    data['individual-id'] = data['animal-id'] + "_" + data['tag-id']
    lt.addLast(data_structs['lobos'],data)
    entry = mp.get(data_structs['orderedData'],data['individual-id'])
    if entry is None:
        lst = lt.newList(datastructure='ARRAY_LIST')
        mp.put(data_structs['orderedData'],data['individual-id'],lst)

def add_data(data_structs, data):
    """
    Función para agregar información
    """
    #TO DO: ordena la información en un mapa (key= id lobo, value: lista de dicc con los eventos del lobo)
    data['individual-id'] = data['individual-local-identifier'] +"_" + data['tag-local-identifier']
    data['location-lat'] = round(float(data['location-lat']),3)
    data['location-long'] = round(float(data['location-long']),3)
    entry = mp.get(data_structs['orderedData'],data['individual-id'])
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
    data_structs['5Vertices'] = lt.newList(datastructure='ARRAY_LIST')
    counter = 0
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
                addConnection(data_structs, origin, destination, distance)
                if counter < 5 and track['node-id'] != lasttrack['node-id']:
                    lt.addLast(data_structs['5Vertices'],track)
                    counter += 1
                if track['location-lat'] > mayorlat:
                    mayorlat = track['location-lat']
                elif track['location-lat'] < menorlat:
                     menorlat = track['location-lat']
                if track['location-long'] > mayorlon:
                    mayorlon = track['location-long']
                elif track['location-long'] < menorlon:
                     menorlon = track['location-long']
            lasttrack = track
    wolfIndividualVertex = lt.size(mp.keySet(data_structs['individualPoints']))
    return wolfIndividualVertex,mayorlat,menorlat,mayorlon,menorlon
        
def formatVertex(data_structs,track):
    """
    Se formatea el nombre del vertice con la longitud y latitud del 
    event seguido del identificador del lobo.
    """
        
    latitud = str(track['location-lat'])
    nl = latitud.replace("-", "m")
    new_lat = nl.replace(".","p")
    
    longitud = str(track['location-long'])
    nlo = longitud.replace("-", "m")
    new_lon = nlo.replace(".","p")
    
    name = new_lon + "_" + new_lat + "_" + track['individual-id']
    position = new_lon + "_" + new_lat
    track['node-id'] = name
    
    entry = mp.get(data_structs['posiciones'],position)
    if entry is None:
        lstWolfsEvents = lt.newList(datastructure='ARRAY_LIST')
        lt.addLast(lstWolfsEvents,track)
        mp.put(data_structs['posiciones'],position,lstWolfsEvents)
        
    else:
        lstWolfsEvents = me.getValue(entry)
        lstIds = lt.newList('ARRAY_LIST')
        for event in lt.iterator(lstWolfsEvents):
            lt.addLast(lstIds,event['individual-id'])
        containsWolfEvent = lt.isPresent(lstIds,track['individual-id'])
        if containsWolfEvent == 0:
            lt.addLast(lstWolfsEvents,track)
            
    entry2 = mp.get(data_structs['individualPoints'],name)
    if entry2 is None:
        lstWolfsEventsInd = lt.newList(datastructure='ARRAY_LIST')
        mp.put(data_structs['individualPoints'],name,lstWolfsEventsInd)
        lt.addLast(lstWolfsEventsInd,track)
    else:
        lstWolfsEventsInd = me.getValue(entry2)
        lstTime = lt.newList('ARRAY_LIST')
        for event in lt.iterator(lstWolfsEventsInd):
            lt.addLast(lstTime,event['timestamp'])
        containsWolfEventInd = lt.isPresent(lstTime,track['timestamp'])
        if containsWolfEventInd == 0:
            lt.addLast(lstWolfsEventsInd,track)
        
    return name

def addPosition(data_structs,vertexName):
    """
    Adiciona una posición como un vertice del grafo
    """
    if not gr.containsVertex(data_structs['grafoDir'], vertexName):
        gr.insertVertex(data_structs['grafoDir'], vertexName)
    if not gr.containsVertex(data_structs['grafoNoDir'], vertexName):
        gr.insertVertex(data_structs['grafoNoDir'], vertexName)
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


def addConnection(data_structs, origin, destination, distance):
    """
    Adiciona un arco entre dos posiciones
    """
    edge1 = gr.getEdge(data_structs['grafoDir'], origin, destination)
    if edge1 is None:
        gr.addEdge(data_structs['grafoDir'], origin, destination, distance)
    edge2 = gr.getEdge(data_structs['grafoNoDir'], origin, destination)
    if edge2 is None:
        gr.addEdge(data_structs['grafoNoDir'], origin, destination, distance)
    

def addPositionConnection(data_structs):
    
    """Por cada vertice (cada posicion) se recorre la lista
    de lobos servidas en dicha posicion y se crean
    arcos entre ellos para representar el cambio de ruta
    que se puede realizar en una posicion."""
    
    lstPosiciones = mp.keySet(data_structs['posiciones'])
    totalWolfsMTPs = 0
    WeightZeroEdges = 0

    for key in lt.iterator(lstPosiciones):
        lstWolfsEvents = mp.get(data_structs['posiciones'], key)['value']
        if lt.size(lstWolfsEvents) >= 2:
            totalWolfsMTPs += lt.size(lstWolfsEvents)
            for event in lt.iterator(lstWolfsEvents):
                vertexName = key + '_' + event['individual-id']
                containsMTP = gr.containsVertex(data_structs['grafoDir'],key)
                if not containsMTP:
                    gr.insertVertex(data_structs['grafoDir'],key)
                    lstMTPs =lt.newList('ARRAY_LIST')
                    lt.addLast(lstMTPs,event)
                    mp.put(data_structs['MTPs'],key,lstMTPs)
                else:
                    entry = mp.get(data_structs['MTPs'],key)
                    lstMTPs = me.getValue(entry)
                    lt.addLast(lstMTPs,event)
                gr.addEdge(data_structs['grafoDir'],key,vertexName,0)
                gr.addEdge(data_structs['grafoDir'],vertexName,key,0)
                WeightZeroEdges += 1
                containsMTP2 = gr.containsVertex(data_structs['grafoNoDir'],key)
                if not containsMTP2:
                    gr.insertVertex(data_structs['grafoNoDir'],key)
                gr.addEdge(data_structs['grafoNoDir'],key,vertexName,0)
                
    totalMTPs = lt.size(mp.keySet(data_structs['MTPs']))
    vertexNum = gr.numVertices(data_structs['grafoDir'])
    return totalMTPs,totalWolfsMTPs,WeightZeroEdges,vertexNum
    
def TabulateCD(data_structs):
    lst = []
    for elem in lt.iterator(data_structs['5Vertices']):
        lstEvent = [] 
        lstEvent.append(elem['location-long'])
        lstEvent.append(elem['location-lat'])
        lstEvent.append(elem['node-id'])
        lstEvent.append(elem['individual-id'])
        lstadjacents = gr.adjacents(data_structs['grafoDir'],elem['node-id'])
        lstEvent.append(lt.size(lstadjacents))
        lst.append(lstEvent)
    
    MTPs = mp.keySet(data_structs['MTPs'])
    LastMTPs = lt.subList(MTPs,(lt.size(MTPs)-5),5)
    for mtp in lt.iterator(LastMTPs):
        lstEvent = []
        elem = mp.get(data_structs['MTPs'],mtp)['value']['elements'][0]
        lstEvent.append(elem['location-long'])
        lstEvent.append(elem['location-lat'])
        lstEvent.append(mtp)
        lstEvent.append(elem['individual-id'])
        lstadjacents = gr.adjacents(data_structs['grafoDir'],mtp)
        lstEvent.append(lt.size(lstadjacents))
        lst.append(lstEvent)
    return lst
        
# Funciones de consulta

def graphSize(data_structs):
    """
    Retorna el tamaño de un grafo (numero total de vértices y nodos)
    """
    #TO DO: Crear la función para obtener un dato de una lista
    totalVertices = gr.numVertices(data_structs['grafoDir'])
    totalEdges = gr.numEdges(data_structs['grafoDir'])
    
    return totalVertices,totalEdges


def data_size(lst):
    """
    Retorna el tamaño de la lista de datos
    """
    #TO DO: Crear la función para obtener el tamaño de una lista
    return lt.size(lst)

def imprimir(control,nodoPrueba):
    print(mp.get(control['MTPs'],'m112p107_56p895'))
    
    
def req_1(data_structs,initialPoint,destPoint):
    """
    Función que soluciona el requerimiento 1
    """
    # TODO: Realizar el requerimiento 1
    data_structs['search']= dfs.DepthFirstSearch(data_structs['grafoDir'],initialPoint) 
    camino = dfs.hasPathTo(data_structs['search'],destPoint)
    lstCamino = lt.newList('ARRAY_LIST')
    if camino:
        ruta = dfs.pathTo(data_structs['search'],destPoint)
        while (not st.isEmpty(ruta)):
            lt.addLast(lstCamino,st.pop(ruta))
            
    if camino == False: 
        return 0,0,[]
    else:
        lastElem = lt.lastElement(lstCamino)
        lt.addLast(lstCamino,lastElem)
        totalDist = 0
        lasttrack = None
        totalMtps = 0
        lstReturn = lt.newList('ARRAY_LIST')
        for nodo in lt.iterator(lstCamino):
            lstInfo = lt.newList('ARRAY_LIST')
            if lasttrack != None:
                try:
                    totalDist += gr.getEdge(data_structs['grafoDir'],lasttrack,nodo)['weight']
                except:
                    totalDist += 0
                entry = mp.get(data_structs['individualPoints'],lasttrack)
                if entry != None:
                    commonWolfs = 1
                    wolfsId = me.getValue(entry)['elements'][0]['individual-id']
                    nodeId = me.getValue(entry)['elements'][0]['node-id']
                else: 
                    entry = mp.get(data_structs['MTPs'],lasttrack)
                    totalMtps += 1
                    value = me.getValue(entry)
                    nodeId = lasttrack
                    commonWolfs = lt.size(value)
                    wolfsIds = lt.newList('ARRAY_LIST')
                    for event in lt.iterator(value):
                        lt.addLast(wolfsIds,event['individual-id'])
                    if lt.size(wolfsIds) > 6:
                        wolfsId = getIFirstandLast(wolfsIds,3)
                    else:
                        lis = []
                        for wolf in lt.iterator(wolfsIds):
                            lis.append(wolf)
                        wolfsId = lis
                value = me.getValue(entry)['elements'][0]
                lt.addLast(lstInfo,nodeId)
                lt.addLast(lstInfo,value['location-long'])
                lt.addLast(lstInfo,value['location-lat'])
                lt.addLast(lstInfo,commonWolfs)
                lt.addLast(lstInfo,wolfsId)
                try:
                    lt.addLast(lstInfo,(gr.getEdge(data_structs['grafoNoDir'],lasttrack,nodo)['weight']))
                except:
                    lt.addLast(lstInfo,0)
                lt.addLast(lstReturn,lstInfo)
            lasttrack = nodo
        
        if lt.size(lstReturn) > 10:
            rta = getiFirstandLast(lstReturn,5)
        else:
            rta = lstReturn
        return totalDist,totalMtps,rta

def req_2(data_structs, initialStation, destination):
    """
    Función que soluciona el requerimiento 2
    """
    # TODO: Realizar el requerimiento 2
    data_structs['search']= bfs.BreadhtFisrtSearch(data_structs['grafoDir'],initialStation)
    hay_camino= bfs.hasPathTo(data_structs['search'],destination)
    num_gathering=0
    tot_dist=0
    ult_nodo= None
    if hay_camino==True:
        lista_camino= lt.newList('ARRAY_LIST')
        lista_camino2= lt.newList('ARRAY_LIST')
        camino= bfs.pathTo(data_structs['search'],destination)
        i=1
        if camino is not None:
            num_vert = st.size(camino)
            while (not st.isEmpty(camino)):
                stop = st.pop(camino)
                lt.addLast(lista_camino, stop)
            for cada_nodo in lt.iterator(lista_camino):
                i+=1
                lista_stop= []
                node_id= cada_nodo
                if ult_nodo!= None and cada_nodo!= destination:
                    tot_dist+= gr.getEdge(data_structs['grafoDir'],cada_nodo, lt.getElement(lista_camino, i))['weight']
                    dist= gr.getEdge(data_structs['grafoDir'],cada_nodo,lt.getElement(lista_camino, i))['weight']      
                else:
                    dist= 0.0
                cant_= cada_nodo.count('_')
                if cant_==1:
                    num_gathering+=1
                    individual_count0=gr.adjacentEdges(data_structs['grafoDir'],cada_nodo)
                    individual_count= lt.size(individual_count0)
                    info_stop = mp.get(data_structs['MTPs'],cada_nodo)
                    info= me.getValue(info_stop)['elements'][1]
                    long= info['location-long']
                    lat=info['location-lat']
                    lista_ady=[]
                    for ady1 in lt.iterator(individual_count0):
                        ad= ady1['vertexB']
                        lista_ady.append(ad)
                    ady= lista_ady
                else:
                        info_stop=mp.get(data_structs['individualPoints'], cada_nodo)
                        info= me.getValue(info_stop)
                        tot_dist+= gr.getEdge(data_structs['grafoDir'],ult_nodo,cada_nodo)['weight']
                        individual_count= 1
                        long= info['elements'][0]['location-long']
                        lat=info['elements'][0]['location-lat']
                        ady=1
                if cada_nodo!=destination:
                        edge_to= lt.getElement(lista_camino, i)
                else:
                    dist= 'Unknown'
                    edge_to='Unknown'
                ult_nodo= cada_nodo
                lista_stop.append(long)
                lista_stop.append(lat)
                lista_stop.append(node_id)
                lista_stop.append(ady)
                lista_stop.append(individual_count)
                lista_stop.append(edge_to)
                lista_stop.append(dist)
                lt.addLast(lista_camino2,lista_stop)
    return lista_camino2


def req_3(data_structs):
    """
    Función que soluciona el requerimiento 3
    """
    # TODO: Realizar el requerimiento 3
    pass


def req_4(data_structs, ini, fin):
    """
    Función que soluciona el requerimiento 4
    """
    # TODO: Realizar el requerimiento 4
    longini= float(ini[0])
    latini= float(ini[1])
    longfin= float(fin[0])
    latfin= float(fin[1])
    dist_menor_ini=9999999999
    mapa= data_structs['MTPs']
    keys= mp.keySet(mapa)
    ino=''
    dest=''
    for llave in lt.iterator(keys):
        lista_mtp= mp.get(mapa, llave)
        lista_mtp= me.getValue(lista_mtp)
        info= lt.getElement(lista_mtp, 1)
        long2=info['location-long']
        lat2= info['location-lat']
        dist= haversine(latini, longini, lat2, long2)
        if dist==0:
            continue
        if dist<dist_menor_ini:
            dist_menor_ini=dist
            ino=llave
    dist_menor_dest=999999999
    for llave in lt.iterator(keys):
            print(llave)
            lista_mtp= mp.get(mapa, llave)
            lista_mtp= me.getValue(lista_mtp)
            info= lt.getElement(lista_mtp, 1)
            long2=info['location-long']
            lat2= info['location-lat']
            dist= haversine( lat2, long2, latfin, longfin)
            print(dist)
            if dist==0:
                continue
            if dist<dist_menor_dest:
                dist_menor_dest=dist
                dest= llave
    paths= djk.Dijkstra(data_structs['grafoDir'], ino)
    hay_path= djk.hasPathTo(paths, dest)
    path = djk.pathTo(paths, dest)
    print(path)
    '''hay= mp.get(mapa,'m111p866_57p451')
    print(hay)'''
    #print(dist_menor_dest, dist_menor_ini, dest, ino)
        

def haversine(lat1, lon1, lat2, lon2):
    lat1, lon1, lat2, lon2 = map(radians, [lat1, lon1, lat2, lon2])
    R = 6371
    dlat = lat2 - lat1
    dlon = lon2 - lon1
    a = sin(dlat / 2) ** 2 + cos(lat1) * cos(lat2) * sin(dlon / 2) ** 2
    c = 2 * atan2(sqrt(a), sqrt(1 - a))
    distance = R * c
    return distance

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

def getiFirstandLast(lst,i):
    newLst = lt.newList('ARRAY_LIST')
    lstAux = lt.newList('SINGLE_LINKED')
    for i in range(0,i):
        first = lt.firstElement(lst)
        lt.addLast(newLst,first)
        lt.removeFirst(lst)
    for j in range(0,i+1):
        last = lt.lastElement(lst)
        lt.addFirst(lstAux,last)
        lt.removeLast(lst)
    for elem in lt.iterator(lstAux):
        lt.addLast(newLst,elem)
    return newLst

def getIFirstandLast(lst,i):
    newLst = []
    lista = getiFirstandLast(lst,i)
    for elem in lt.iterator(lista):
        newLst.append(elem)
    return newLst
