#!/usr/bin/env python
#! python/bin/user
import sys
import sets
import random
import collections
import igraph
import signal
import time
import igraph
import itertools

########## Start : Expected by OPTIL #############
class Killer:
    exit_now = False
    def __init__(self):
        signal.signal(signal.SIGINT, self.exit)
        signal.signal(signal.SIGTERM, self.exit)

    def exit(self,signum, frame):
        try:
            #print "Exited due to kill signal"
            #print "Number of time Exp Ran", COUNTER
            print "VALUE", Min_wt_st_random
            for x, y in Min_tree_edges_random:
                print x, y
        except ValueError:
            pass
    
        sys.exit()
#self.exit_now = True
########## End : Expected by OPTIL #############


########## Start: Graph Algorithms ############
def merge_vertices(G, v1_name, v2_name):
    """
    Input: Graph G and two of its vertices v1_name, v2_name
    Output: Nothing. Modifies the input grpah

    This function will add an edge between every vertex in neighborhood of v2_name with v1_name. It then simplify the graph and delete v2_name.
    """ 
    # Checking whether v1_name and v_name are in G
    try:
        v1 = G.vs.find(name=v1_name)
        v2 = G.vs.find(name=v2_name)
    except:
        raise ValueError("Vertices not found in graph")

    for u in G.neighbors(v2):
        G.add_edge(v1.index, u)

    G.delete_vertices([v2])
    G.simplify()

########## End: Graph Algorithms ##############

########## Start: Heuristics ##################    
    
def approx_mst(G):
    """
    Input: Graph G
    Output: Minimum weight spanning tree of G
    """
    Edge_Weights = [e['weight'] for e in G.es()]
    Tree_edge_ids = G.spanning_tree(weights=Edge_Weights, return_tree=False)
    wt_steiner_tree = 0
    for edge_id in Tree_edge_ids:
        wt_steiner_tree += G.es()[edge_id]['weight']
    return wt_steiner_tree, Tree_edge_ids


def approx_mst_pruning(G, Terminals):
    """
    Input: Graph G
    Output: Weight and edges of Minimum weight spanning tree of G.

    Delete remaining edges and prone the tree to make sure that only leaves are terminals
    """
    Edge_Weights = [e['weight'] for e in G.es()]
    Edge_IDs = [e.index for e in G.es()]
    Tree_edge_ids = G.spanning_tree(weights=Edge_Weights, return_tree=False)
    Set_tree_edge_ids = set(Tree_edge_ids)
    Edge_IDs_to_delete = [edge_id for edge_id in Edge_IDs if edge_id not in Set_tree_edge_ids]
    G.delete_edges(Edge_IDs_to_delete)
    Terminal_names = set([v["name"] for v in Terminals])
    Pendent_vertices = [v for v in G.vs() if (G.degree(v) == 1 and not v["name"] in Terminal_names)]
    while Pendent_vertices:
        G.delete_vertices(Pendent_vertices)
        Pendent_vertices = [v for v in G.vs() if (v.degree() == 1 and not v["name"] in Terminal_names)]

    wt_steiner_tree = 0
    Tree_edges = []
    for edge in G.es():
        wt_steiner_tree += edge['weight']
        u = G.vs()[edge.source]["name"]
        v = G.vs()[edge.target]["name"]
        Tree_edges.append((u, v))
        
    return wt_steiner_tree, Tree_edges


def random_select(G, Terminals):
    """
    Input: Graph G, set of Terminals T
    Output: Weight and edges of Steiner Tree

    Randomly select two terminals. Find shortest path connecting them. Contract this path.
    """
    Terminal_names = set([v["name"] for v in Terminals])
    wt_steiner_tree = 0
    Tree_edges = set([])
    while len(Terminal_names) > 1:
        t1_name, t2_name = random.sample(Terminal_names, 2)
        Terminal_names.remove(t2_name)
        t1 = G.vs().find(name=t1_name)
        t2 = G.vs().find(name=t2_name)
        shortest_epath = G.get_shortest_paths(t1, t2, weights="weight", output="epath")[0]
        for edge_id in shortest_epath:
            edge = G.es()[edge_id]
            u = G.vs()[edge.source]
            v = G.vs()[edge.target]
            # If edge (u, v) or (v, u) is present in Tree_edges then do nothing.
            if (v["name"], u["name"]) in Tree_edges or (u["name"], v["name"]) in Tree_edges:
                continue
            Tree_edges.add((u["name"], v["name"]))
            wt_steiner_tree += edge["weight"]
            edge["weight"] = 0
            
    return wt_steiner_tree, Tree_edges

########## End: Heuristics ####################    


########## Start: Kernelization ###############

def rr_deg1(G, Terminals):
    """
    Input : Graph G
    Output : Nothing (Modifies the input graph)

    Applies reduction rule.
    Reduction Rule : Delete all non-terminal vertices of deg at most 1.
    """ 
    small_degree = [v for v in G.vs() if (G.degree(v) <= 1 and not v in Terminals)]
    G.delete_vertices(small_degree)

def rr_deg2(G, Terminals):
    """
    Input: Graph G, set Terminals of vertices
    Output: Nothing (Modifies the input graph)
    
    If a non-terminal vertex v exactly two neighbors, say x, y then the reduction rule deletes v and add edge xy which is of weight min{wt(xy), wt(vx) + wt(vy)}. 
    """
    deg2_vertices_names = [v["name"] for v in G.vs() if G.degree(v) == 2]
    for v_name in deg2_vertices_names:
        v = G.vs().find(name=v_name)
        #print G.degree(v)
        if v in Terminals:
            continue
        if G.degree(v) == 1:
            G.delete_vertices(v)
            continue
		
        x, y = G.vs()[G.neighbors(v)[0]], G.vs()[G.neighbors(v)[1]]
        wt_vx = G.es()[G.get_eid(v, y)]["weight"]
        wt_vy = G.es()[G.get_eid(v, y)]["weight"]
        if G.are_connected(x, y):
            #Get edge id if they are connecte. Otherwise handle in error.
            wt_xy = G.es()[G.get_eid(x, y)]["weight"] 
            G.es()[G.get_eid(x, y)]["weight"] = min(wt_xy, wt_vx + wt_vy)
        else:
            G.add_edge(x, y, weight = wt_vx + wt_vy)
        G.delete_vertices(v)
    
########### End: Kernelization ##############
    
    
########### Star : Steiner Tree #############

g = igraph.Graph()
Terminals = set([])
killer = Killer()

#Reading input
while True:
    try:
        str = raw_input()       
    except EOFError:
        break

    if str == "":
        continue
  
    String_break = str.split()
    if String_break[0] == "Section":
        continue
  
    if String_break[0] == "E":
        #x, y, wt = map(int, String_break[1:])
        x, y = String_break[1], String_break[2]
        wt = int(String_break[3])
        
        try:
            g.vs.find(name=x)
        except ValueError:
            g.add_vertex(name=x)

        try:
            g.vs.find(name=y)
        except ValueError:
            g.add_vertex(name=y)
            
        #Add edge with weight
        g.add_edge(x, y, weight = wt)

    if String_break[0] == "T":
        x = String_break[1]
        Terminals.add(g.vs().find(name=x))
      
    if String_break[0] == "Nodes":
        N = String_break[1]

    if String_break[0] == "Edges":
        E = String_break[1]
          
    if String_break[0] == "Terminals":
        T = String_break[1]
              
    if String_break[0] == "EOF":
        break


#print g.vcount(), "\t", g.ecount(), "\t", len(Terminals)
#print "Total_edge weight", sum([e['weight'] for e in g.es()])
    
#pendent_vertices = [v for v in g.vs() if g.degree(v) <= 1]
#deg2_vertices = [v for v in g.vs() if g.degree(v) == 2]

#while [v for v in g.vs() if (g.degree(v) <= 2 and not v in Terminals)]:
#    rr_deg1(g, Terminals)
#    rr_deg2(g, Terminals)

#print g.vcount(), "\t", g.ecount(), "\t", len(Terminals)
#wt_steinter_tree, tree_edge_ids = approx_mst(g)

g1 = g.copy()
#Min_wt_st_random, Min_tree_edges_random = random_select(g1, Terminals)

Min_wt_st_random, Min_tree_edges_random = approx_mst_pruning(g1, Terminals)

#COUNTER = 0 # Measures the number of times experiment is run
while True:
    #COUNTER += 1
    g1 = g.copy()
    wt_st_random, tree_edges_random = random_select(g1, Terminals)
    if Min_wt_st_random > wt_st_random:
        Min_wt_st_random = wt_st_random
        Min_tree_edges_random = tree_edges_random

#print "VALUE", wt_st_random

#wt_st1, tree_edges = approx_mst_pruning(g1, Terminals)
#print "VALUE prunning", wt_st1

#for x, y in tree_edges_random:
#    print x, y

#print wt_st1
#for edge_id in tree_edge_ids:
#    edge = g.es()[edge_id]
    #print g.are_connected(x, y)
#    x, y = edge.source, edge.target
#    v = g.vs().find(name=x)
#    u = g.vs().find(name=y)
#    print u["name"], v["name"]
    #tot += g.es()[g.get_eid(v, u)]["weight"]
    
sys.exit()
########### End : Steiner Tree #############
