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
from tabulate import tabulate
from math import sin, cos, sqrt, atan2, radians
default_limit = 1000
sys.setrecursionlimit(default_limit*10000)

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
    data_structs['lobos'] = mp.newMap(maptype='PROBING')
    data_structs['MTPs'] = mp.newMap(maptype='PROBING')
    data_structs['individualPoints'] = mp.newMap(maptype='PROBING')
    

    return data_structs

# Funciones para agregar informacion al modelo
def addWolfsData(data_structs, data):
    data['individual-id'] = data['animal-id'] + "_" + data['tag-id']
    if data['deploy-off-date']== '':
        data['deploy-off-date']= 'Unknown'
    if data['animal-death-comments']== '':
        data['animal-death-comments']= 'Unknown'
    if data['animal-life-stage']== '':
        data['animal-life-stage']='Unknown'
    if data['animal-sex']== '':
        data['animal-sex']= 'Unknown'
    if data['deployment-comments']== '':
        data['deployment-comments']= 'Unknown'
    if data['study-site']=='':
        data['study-site']= 'Unknown'
    if data['tag-beacon-frequency']=='':
        data['tag-beacon-frequency']= 'Unknown'
    if data['tag-mass']=='':
        data['tag-mass']= 'Unknown'
    if data['tag-model']=='':
        data['tag-model']= 'Unknown'
    mp.put(data_structs['lobos'],data['individual-id'],data)
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
    track['coordenada'] = position
    
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
    trackingPoints = 0
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
                    trackingPoints += 1
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
                        wolfsIds = getiFirstandLast(wolfsIds,3)
                    res = ""
                    for wolf in lt.iterator(wolfsIds):
                        res += wolf+", "
                        wolfsId = res
                value = me.getValue(entry)['elements'][0]
                lt.addLast(lstInfo,nodeId)
                lt.addLast(lstInfo,value['location-long'])
                lt.addLast(lstInfo,value['location-lat'])
                lt.addLast(lstInfo,commonWolfs)
                lt.addLast(lstInfo,wolfsId)
                if lasttrack != nodo:
                    lt.addLast(lstInfo,nodo)
                else:
                    lt.addLast(lstInfo,'Unknown')
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
        totalNodes = lt.size(lstCamino) 
        return totalNodes-1,trackingPoints,totalDist,totalMtps,rta

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
                    info= me.getValue(info_stop)['elements'][0]
                    long= info['location-long']
                    lat=info['location-lat']
                    lista_ady=[]
                    for ady1 in lt.iterator(individual_count0):
                        ad= ady1['vertexB']
                        info_stop=mp.get(data_structs['individualPoints'], ad)
                        info= me.getValue(info_stop)
                        ady2= info['elements'][0]['individual-id']
                        lista_ady.append(ady2)
                    if len(lista_ady) > 6:
                        lista_ady = ult3_prim3(lista_ady)
                    else:
                        lista_ady = lista_ady
                    ady= lista_ady
                else:
                        info_stop=mp.get(data_structs['individualPoints'], cada_nodo)
                        info= me.getValue(info_stop)
                        individual_count= 1
                        long= info['elements'][0]['location-long']
                        lat=info['elements'][0]['location-lat']
                        ady= info['elements'][0]['individual-id']
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
    num_track= num_vert - num_gathering
    edges= num_vert- 1
    tot_dist= round(tot_dist, 4)
    if lt.size(lista_camino2) > 10:
        lista_camino2 = getiFirstandLast(lista_camino2,5)
    else:
        lista_camino2 = lista_camino2
    return lista_camino2, num_gathering, num_vert, num_track, edges, tot_dist




def req_3(data_structs):
    """
    Función que soluciona el requerimiento 3
    """
    # TODO: Realizar el requerimiento 3
    sccStruct = scc.KosarajuSCC(data_structs['grafoDir'])
    numScc = scc.connectedComponents(sccStruct)
    idSccMap = mp.newMap(maptype='PROBING')
    for nodo in lt.iterator(mp.keySet(sccStruct['idscc'])):
        idScc = mp.get(sccStruct['idscc'],nodo)['value']
        contains = mp.get(idSccMap,idScc)
        if contains == None:
            lstNodes = lt.newList('ARRAY_LIST')
            lt.addLast(lstNodes,nodo)
            mp.put(idSccMap,idScc,lstNodes)
        else:
            lstNodes = me.getValue(contains)
            lt.addLast(lstNodes,nodo)
            
    omSize = om.newMap()
    for idScc in lt.iterator(mp.keySet(idSccMap)):
        lstNodes = me.getValue(mp.get(idSccMap,idScc))
        if lt.size(lstNodes) > 2:
            info = mp.newMap(maptype='PROBING')
            mp.put(info,'idScc',idScc)
            mp.put(info,'nodes',lstNodes)
            om.put(omSize,lt.size(lstNodes),info)
    c =0
    FivemaxManInfo = []
    while c < 5:
        minLat = 1000
        maxLat = 0
        minLong = 0
        maxLong = -1000
        infoScc = []        
        maxM= om.maxKey(omSize)
        info = me.getValue(om.get(omSize,maxM))
        sccId = me.getValue(mp.get(info,'idScc'))
        lstNodesId = me.getValue(mp.get(info,'nodes'))
        nodeIds = getIFirstandLast(lstNodesId,3)
        wolfs = lt.newList('ARRAY_LIST')
        for nodo in lt.iterator(lstNodesId):
            entry = mp.get(data_structs['MTPs'],nodo)
            if entry == None:
                entry = mp.get(data_structs['individualPoints'],nodo)
                wolf = me.getValue(entry)['elements'][0]['individual-id']
                if not lt.isPresent(wolfs,wolf):
                    lt.addLast(wolfs,wolf)
            event = me.getValue(entry)['elements'][0]
            if event['location-lat'] > maxLat:
                maxLat = event['location-lat']
            elif event['location-lat'] < minLat:
                minLat = event['location-lat']
            if event['location-long'] > maxLong:
                maxLong = event['location-long']
            elif event['location-long'] < minLong:
                minLong = event['location-long']    
        infoScc.append(sccId)
        res = []
        p =0
        for Id in nodeIds:
            lst = [Id]
            res.append(lst)
            p+= 1
            if p ==3:
                pt = ['...']
                res.append(pt) 
        infoScc.append(tabulate(res,tablefmt="plain"))
        infoScc.append(maxM)
        infoScc.append(minLat)
        infoScc.append(maxLat)
        infoScc.append(minLong)
        infoScc.append(maxLong)
        infoScc.append(lt.size(wolfs))
        lstOflst =[]
        if lt.size(wolfs) > 6:
            wolfs = getiFirstandLast(wolfs,3)
        for lobo in lt.iterator(wolfs):
            infoLobo = []
            details = me.getValue(mp.get(data_structs['lobos'],lobo))
            infoLobo.append(details['individual-id'])
            infoLobo.append(details['animal-taxon'])
            infoLobo.append(details['animal-sex'])
            infoLobo.append(details['animal-life-stage'])
            infoLobo.append(details['study-site'])
            lstOflst.append(infoLobo)
        wolftable = tabulate(lstOflst,headers=['indiv-id','wolf taxon','wolf sex','life-stage','study-site'],
                             tablefmt='grid',maxheadercolwidths=5,maxcolwidths=5,stralign="center")
        infoScc.append(wolftable)
        FivemaxManInfo.append(infoScc)
        om.deleteMax(omSize)
        c += 1
    return numScc, FivemaxManInfo
    
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
            lista_mtp= mp.get(mapa, llave)
            lista_mtp= me.getValue(lista_mtp)
            info= lt.getElement(lista_mtp, 1)
            long2=info['location-long']
            lat2= info['location-lat']
            dist= haversine( lat2, long2, latfin, longfin)
            if dist==0:
                continue
            if dist<dist_menor_dest:
                dist_menor_dest=dist
                dest= llave
    info_stop = mp.get(data_structs['MTPs'],ino)
    info= me.getValue(info_stop)
    longino= info['elements'][0]['location-long']
    latino=info['elements'][0]['location-lat']
    individual_idino=gr.adjacentEdges(data_structs['grafoDir'],ino)
    lista_adyini= []
    for ady1 in lt.iterator(individual_idino):
        ad= ady1['vertexB']
        info_stop=mp.get(data_structs['individualPoints'], ad)
        info= me.getValue(info_stop)
        ady2= info['elements'][0]['individual-id']
        lista_adyini.append(ady2)
    individual_idino= lista_adyini

    info_stop = mp.get(data_structs['MTPs'],dest)
    info= me.getValue(info_stop)
    longdest= info['elements'][0]['location-long']
    latdest=info['elements'][0]['location-lat']
    individual_id_dest=gr.adjacentEdges(data_structs['grafoDir'],dest)
    lista_adyini= []
    for ady1 in lt.iterator(individual_id_dest):
        ad= ady1['vertexB']
        info_stop=mp.get(data_structs['individualPoints'], ad)
        info= me.getValue(info_stop)
        ady2= info['elements'][0]['individual-id']
        lista_adyini.append(ady2)
    individual_id_dest= lista_adyini
    paths= djk.Dijkstra(data_structs['grafoDir'], ino)
    hay_path= djk.hasPathTo(paths, dest)
    if hay_path ==True:
        path = djk.pathTo(paths, dest)
        lista_camino= lt.newList('ARRAY_LIST')
        lista_gathering= lt.newList('ARRAY_LIST')
        lista_camino2= lt.newList('ARRAY_LIST')
        dist_tot=0
        if path is not None:
            num_vert = st.size(path)
            while (not st.isEmpty(path)):
                stop = st.pop(path)
                lt.addLast(lista_camino, stop)
        for info in lt.iterator(lista_camino):
            lista= []
            src_node_id= info['vertexA']
            tgt_node_id= info['vertexB']
            dist= info[ 'weight']
            dist_tot+= info[ 'weight']
            src_count=src_node_id.count('_')
            tgt_count=tgt_node_id.count('_')
            if src_count==1:
                lista_gath= []
                info_stop = mp.get(data_structs['MTPs'],src_node_id)
                info= me.getValue(info_stop)['elements'][0]
                longsrc= info['location-long']
                latsrc=info['location-lat']
                lt.addLast(lista_gathering, src_node_id)
            else:
                info_stop=mp.get(data_structs['individualPoints'], src_node_id)
                info= me.getValue(info_stop)
                longsrc= info['elements'][0]['location-long']
                latsrc=info['elements'][0]['location-lat']
            if tgt_count==1:
                info_stop = mp.get(data_structs['MTPs'],tgt_node_id)
                info= me.getValue(info_stop)['elements'][0]
                longtgt= info['location-long']
                lattgt=info['location-lat']
                individual_idtgt= info['individual-id']
                lt.addLast(lista_gathering, tgt_node_id)
            else:
                info_stop=mp.get(data_structs['individualPoints'], tgt_node_id)
                info= me.getValue(info_stop)
                longtgt= info['elements'][0]['location-long']
                lattgt=info['elements'][0]['location-lat']
                individual_idtgt= info['elements'][0]['individual-id']
            lista.append(src_node_id)
            lista.append(latsrc)
            lista.append(longsrc)
            lista.append(tgt_node_id)
            lista.append(lattgt) 
            lista.append(longtgt)
            lista.append(individual_idtgt) 
            lista.append(dist)
            lt.addLast(lista_camino2, lista)
    listagathering= lt.newList('ARRAY_LIST')
    for mtp in lt.iterator(lista_gathering):
        lista1=[]
        info_stop = mp.get(data_structs['MTPs'],ino)
        info= me.getValue(info_stop)
        longmtp= info['elements'][0]['location-long']
        latmtp=info['elements'][0]['location-lat']
        individual_idmtp=gr.adjacentEdges(data_structs['grafoDir'],mtp)
        lista_adymtp= lt.newList('ARRAY_LIST')
        for ady1 in lt.iterator(individual_idmtp):
            ad= ady1['vertexB']
            info_stop=mp.get(data_structs['individualPoints'], ad)
            info= me.getValue(info_stop)
            ady2= info['elements'][0]['individual-id']
            lt.addLast(lista_adymtp,ady2)
        tam= lt.size(individual_idmtp)
        individual_idmtp= lista_adymtp
        if lt.size(individual_idmtp) > 6:
            individual_idmtp = getiFirstandLast(individual_idmtp,3)
        else:
            individual_idmtp = individual_idmtp
        lista1.append(mtp)
        lista1.append(longmtp)
        lista1.append(latmtp)
        lista1.append(individual_idmtp)
        lista1.append(tam)
        lt.addLast(listagathering, lista1)
    num_arc= num_vert-1
    size_gath= lt.size(lista_gathering)
    tot_dist= dist_tot+ dist_menor_ini+ dist_menor_dest
    if lt.size(listagathering) > 6:
        listagathering = getiFirstandLast(listagathering,3)
    else:
        listagathering = listagathering
    return lista_camino2, listagathering, num_arc, num_vert, size_gath, dist_tot, dist_menor_ini, dist_menor_dest, tot_dist, ino, longino, latino,individual_idino, dest, longdest, latdest, individual_id_dest
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

def req_5(data_structs, identificador, maxDis, minPuntos):
    """
    Función que soluciona el requerimiento 5
    """
    # TODO: Realizar el requerimiento 5
    listaNodos = lt.newList(datastructure='ARRAY_LIST')
    resu_djk = djk.Dijkstra(data_structs["grafoNoDir"], identificador)
    for nodo in resu_djk['iminpq']['elements']['elements']:
            if nodo['index'] <= maxDis:
                lt.addLast(listaNodos, nodo)
    listaNodosCondiciones = lt.newList(datastructure='ARRAY_LIST')
    if minPuntos != 0:
        for final in lt.iterator(listaNodos):
            camino = djk.pathTo(resu_djk, final["key"])
            path = camino['first']
            cant = 1
            if path:
                listaPuntos = path['info']['vertexA']
                path = path['next']
                while path:
                    listaPuntos += ", " + path['info']['vertexA']
                    cant += 1
                    if 'last' in path.keys():
                        path = path['last']
                    else:
                        path = path['next']
            if cant >= minPuntos:
                lt.addLast(listaNodosCondiciones, {"Point Count" : cant, "Path distance [km]": final["index"], "Point List": listaPuntos})
    return listaNodosCondiciones


def req_6(data_structs,animal_sex, ini, fin):
    """
    Función que soluciona el requerimiento 6
    """
    # TODO: Realizar el requerimiento 6
    lista_lobos= lt.newList('ARRAY_LIST')
    lbs= data_structs['lobos']
    keys= mp.keySet(lbs)
    mapa_filt= mp.newMap(maptype='CHAINING')
    for lobo in lt.iterator(keys):
        entry= mp.get(lbs, lobo)
        info= me.getValue(entry)
        a_s= info['animal-sex']
        if a_s == animal_sex:
            lt.addLast(lista_lobos,lobo )
    mapa= data_structs['orderedData']
    for lobo in lt.iterator(lista_lobos):
        entry= mp.get(mapa, lobo)
        info= me.getValue(entry)
        for inf in lt.iterator(info):
            fecha= inf['timestamp']
            fecha= datetime.datetime.strptime(fecha,'%Y-%m-%d %H:%M')
            if fecha > ini and fecha< fin:
                if mp.contains(mapa_filt, lobo):
                    por_fecha=mp.get(mapa_filt, lobo)
                    por_fecha=me.getValue(por_fecha)
                    lt.addLast(por_fecha, inf)
                else:
                    por_anio= lt.newList('ARRAY_LIST')
                    lt.addLast(por_anio, inf)
                    mp.put(mapa_filt,lobo, por_anio)
                    
    omDist = om.newMap()
    for wolf in lt.iterator(mp.keySet(mapa_filt)):
        distance = 0
        lstEvents = me.getValue(mp.get(mapa_filt,wolf))
        lasttrack = None
        for event in lt.iterator(lstEvents):
            if lasttrack != None:
                startP = lasttrack['node-id']
                destP = event['node-id']
                distance += gr.getEdge(data_structs['grafoDir'],startP,destP)['weight']
            lasttrack = event
        om.put(omDist,distance,lstEvents)
     
    maxDist = om.maxKey(omDist)        
    minDist = om.minKey(omDist) 
    wolfMaxDist = me.getValue(om.get(omDist,maxDist))['elements'][0]['individual-id']
    wolfMinDist = me.getValue(om.get(omDist,minDist))['elements'][0]['individual-id']
    maxWolfInfo = ObtainWolfInfo(data_structs,wolfMaxDist,maxDist)
    minWolfInfo = ObtainWolfInfo(data_structs,wolfMinDist,minDist)
    #segunda parte
    pathInfoMax = pathInfo(data_structs,omDist,wolfMaxDist,maxDist)
    pathInfoMin = pathInfo(data_structs,omDist,wolfMinDist,minDist)
    
    return wolfMaxDist, maxWolfInfo, pathInfoMax,wolfMinDist,minWolfInfo,pathInfoMin

def req_7(data_structs,dateI, dateF,tempMax,tempMin):
    """
    Función que soluciona el requerimiento 7
    """
    # TO DO: Realizar el requerimiento 7
    
    mapaFilt = mp.newMap()
    for wolf in lt.iterator(mp.keySet(data_structs['orderedData'])):
        lstWolfs = lt.newList('ARRAY_LIST')
        for track in lt.iterator(me.getValue(mp.get(data_structs['orderedData'],wolf))):
            fecha = track['timestamp'] 
            fecha= datetime.datetime.strptime(fecha,'%Y-%m-%d %H:%M')
            if float(track['external-temperature']) >= tempMin and float(track['external-temperature']) <= tempMax and fecha >= dateI and fecha <= dateF:
                lt.addLast(lstWolfs,track)
        mp.put(mapaFilt,wolf,lstWolfs)
    
    mapaPosiciones = mp.newMap()
    graFilt = grafoFilt(mapaFilt,mapaPosiciones)
    for pos in lt.iterator(mp.keySet(mapaPosiciones)):
        if lt.size(me.getValue(mp.get(mapaPosiciones,pos))) >= 2:
            gr.insertVertex(graFilt,pos)
            for node in lt.iterator(me.getValue(mp.get(mapaPosiciones,pos))):
                gr.addEdge(graFilt,node['node-id'],pos,0)
                gr.addEdge(graFilt,pos,node['node-id'],0)
    #para los componentes
    sccStruct = scc.KosarajuSCC(graFilt)
    numScc = scc.connectedComponents(sccStruct)
    idSccMap = mp.newMap(maptype='PROBING')
    for nodo in lt.iterator(mp.keySet(sccStruct['idscc'])):
        idScc = mp.get(sccStruct['idscc'],nodo)['value']
        contains = mp.get(idSccMap,idScc)
        if contains == None:
            lstNodes = lt.newList('ARRAY_LIST')
            lt.addLast(lstNodes,nodo)
            mp.put(idSccMap,idScc,lstNodes)
        else:
            lstNodes = me.getValue(contains)
            lt.addLast(lstNodes,nodo)
            
    omSize = om.newMap()
    for idScc in lt.iterator(mp.keySet(idSccMap)):
        lstNodes = me.getValue(mp.get(idSccMap,idScc))
        if lt.size(lstNodes) > 2:
            info = mp.newMap(maptype='PROBING')
            mp.put(info,'idScc',idScc)
            mp.put(info,'nodes',lstNodes)
            om.put(omSize,lt.size(lstNodes),info)
            
    c =0
    FivemaxManInfo = []
    while c < 3:
        minLat = 1000
        maxLat = 0
        minLong = 0
        maxLong = -1000
        infoScc = []        
        maxM= om.maxKey(omSize)
        info = me.getValue(om.get(omSize,maxM))
        sccId = me.getValue(mp.get(info,'idScc'))
        lstNodesId = me.getValue(mp.get(info,'nodes'))
        nodeIds = getIFirstandLast(lstNodesId,3)
        wolfs = lt.newList('ARRAY_LIST')
        for nodo in lt.iterator(lstNodesId):
            entry = mp.get(data_structs['MTPs'],nodo)
            if entry == None:
                entry = mp.get(data_structs['individualPoints'],nodo)
                wolf = me.getValue(entry)['elements'][0]['individual-id']
                if not lt.isPresent(wolfs,wolf):
                    lt.addLast(wolfs,wolf)
            event = me.getValue(entry)['elements'][0]
            if event['location-lat'] > maxLat:
                maxLat = event['location-lat']
            elif event['location-lat'] < minLat:
                minLat = event['location-lat']
            if event['location-long'] > maxLong:
                maxLong = event['location-long']
            elif event['location-long'] < minLong:
                minLong = event['location-long']    
        infoScc.append(sccId)
        res = []
        p =0
        for Id in nodeIds:
            lst = [Id]
            res.append(lst)
            p+= 1
            if p ==3:
                pt = ['...']
                res.append(pt) 
        infoScc.append(tabulate(res,tablefmt="plain"))
        infoScc.append(maxM)
        infoScc.append(minLat)
        infoScc.append(maxLat)
        infoScc.append(minLong)
        infoScc.append(maxLong)
        infoScc.append(lt.size(wolfs))
        lstOflst =[]
        if lt.size(wolfs) > 6:
            wolfs = getiFirstandLast(wolfs,3)
        for lobo in lt.iterator(wolfs):
            infoLobo = []
            details = me.getValue(mp.get(data_structs['lobos'],lobo))
            infoLobo.append(details['individual-id'])
            infoLobo.append(details['animal-taxon'])
            infoLobo.append(details['animal-sex'])
            infoLobo.append(details['animal-life-stage'])
            infoLobo.append(details['study-site'])
            lstOflst.append(infoLobo)
        wolftable = tabulate(lstOflst,headers=['indiv-id','wolf taxon','wolf sex','life-stage','study-site'],
                             tablefmt='grid',maxheadercolwidths=5,maxcolwidths=5,stralign="center")
        infoScc.append(wolftable)
        FivemaxManInfo.append(infoScc)
        om.deleteMax(omSize)
        c += 1
    return numScc, FivemaxManInfo
    
def grafoFilt(mapaFilt,mapaPosiciones):
    grafoFilt = gr.newGraph(datastructure= "ADJ_LIST",directed=True)
    for wolf in lt.iterator(mp.keySet(mapaFilt)):
        entry = mp.get(mapaFilt,wolf)
        wolfEventsLst = me.getValue(entry)
        lasttrack = None
        for track in lt.iterator(wolfEventsLst):
            if lasttrack != None:
                origin = lasttrack['node-id']
                destination = track['node-id']
                if not gr.containsVertex(grafoFilt,origin):
                    gr.insertVertex(grafoFilt,origin)
                if not gr.containsVertex(grafoFilt,destination):
                    gr.insertVertex(grafoFilt,destination)
                distance = getDistance(track,lasttrack)
                if gr.getEdge(grafoFilt,origin,destination) ==  None:
                    gr.addEdge(grafoFilt,origin,destination,distance)
                try:
                    if mp.get(mapaPosiciones,lasttrack['coordenada']) == None:
                        lstTracks = lt.newList('ARRAY_LIST')
                        mp.put(mapaPosiciones,track['coordenada'],lstTracks)
                    lstTracks = me.getValue(mp.get(mapaPosiciones,track['coordenada']))
                    lstwolfs = lt.newList()
                    for ev in lt.iterator(lstTracks):
                        lt.addLast(lstwolfs,ev['individual-id'])
                    if not lt.isPresent(lstwolfs,lasttrack['individual-id']):
                        lt.addLast(lstTracks,track)
                except:
                    None
            lasttrack = track
    return grafoFilt

def req_8(data_structs):
    """
    Función que soluciona el requerimiento 8
    """
    # TODO: Realizar el requerimiento 8
    pass

def ult3_prim3(lista):
    lsta=[]
    l1= lista[0]
    l2= lista[1]
    l3= lista[2]
    ult1= lista[len(lista)-3]
    ult2= lista[len(lista)-2]
    ult3= lista[len(lista)-1]
    lsta.append(l1)
    lsta.append(l2)
    lsta.append(l3)
    lsta.append(ult1)
    lsta.append(ult2)
    lsta.append(ult3)
    return lsta
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
    """
    Retorna una lista DiscLib
    """
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
    """
    Retorna una lista python a partir de una lista DiscLib
    """
    newLst = []
    lista = getiFirstandLast(lst,i)
    for elem in lt.iterator(lista):
        newLst.append(elem)
    return newLst

def ObtainWolfInfo(data_structs, wolfId,dist):
    """
    Retorna una lista de listas con la info del lobo
    """
    details = me.getValue(mp.get(data_structs['lobos'],wolfId))
    dist =  round(dist,3)
    lstOfLst = [wolfId,details['animal-taxon'],details['animal-life-stage'],details['animal-sex'],details['study-site'],
                str(dist),details['deployment-comments']]
    return lstOfLst         
   
def pathInfo(data_structs,omMap,Id,dist):
    totalNodes =  lt.size(me.getValue(om.get(omMap,dist)))
    totalEdges = totalNodes-1
    if lt.size(me.getValue(om.get(omMap,dist))) > 6:
        nodes = getIFirstandLast(me.getValue(om.get(omMap,dist)),3)
    else:
        nodes = me.getValue(om.get(omMap,dist))['elements']
    nodeInfo = nodesInfo(nodes,data_structs)
    return totalNodes, totalEdges, nodeInfo
    
def nodesInfo(lstTracks,data_struct):
    lstOflst = []
    for track in lstTracks:
        entry = mp.get(data_struct['individualPoints'],track['node-id'])
        wolfs = me.getValue(entry)['elements'][0]['individual-id']
        node_id = track['node-id']
        wolfCount = 1
        if entry == None:
            entry = mp.get(data_struct['MTPs'],track['node-id'])
            print('punto de encuentro!')
            wolfs = []
            wolfCount = 0
            node_id = track['node-id']
            for ev in lt.iterator(me.getValue(entry)):
                wolfs.append(ev['individual-id'])
                wolfCount += 1
            
        event = me.getValue(entry)['elements'][0]
        lstInfo = [node_id,event['location-long'],event['location-lat'],wolfs,wolfCount]
        lstOflst.append(lstInfo)
    return lstOflst