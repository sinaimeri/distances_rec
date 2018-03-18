######################################################################
# Different distances between reconciliations are developed and test #
######################################################################

from sys  import argv
import pdb
from ete3 import Tree
import time
import random
import itertools
import sys
######import the things for the triplet for multifurcations
import csv
import numpy as np
from numpy import *
import scipy
#from pyTQDist import tripletDistance,allPairsTripletDistance,pairsTripletDistance


def usage():
    print "USAGE: "
# print argv[0], "file_with_reconciliations.txt", "file_with_corresponding_distances_AllAgainstAll.txt"
# print "The file of reconciliations should be the same format of the output of Eucalypt tool"


##############e#######################################################
#            Definition of functions to manipulate trees             #
######################################################################


#calculates the diamenter of a tree. To be fasten by implementing a BFS visit.
def diameter(T):
    diameter=0
    current_dist=0
    #get the set of leaves in T
    Leaves=set()
    Leaves=T.get_leaf_names(is_leaf_fn=None)
    #iterate over all pairs of leaves
    for item in list(itertools.combinations(Leaves, 2)):
        (a,b)=item
        current_dist=T.get_distance(a,b, topology_only=True) + 1
        if (diameter< current_dist):
            diameter = current_dist
    print diameter
    return diameter




##############e#######################################################
#            Definition of Different metrics                         #
#            - path metric                                           #
#            - triplet based metric                                  #
#            - Robinson-Foulds  metric                               #
######################################################################

def path_metric(rec1, rec2, HostTree):
    d=0
    for i in range(0,len(rec1)):
        if (rec1[i][0] != rec2[i][0]):
            print "Error: The two reconcilaitions are not sorted in the same way"

        elif (rec1[i][1] != rec2[i][1]):
            node1=HostTree.search_nodes(name=rec1[i][1])[0]
            node2=HostTree.search_nodes(name=rec2[i][1])[0]

            #get_distance gets the number of nodes between two nodes. As the distance is defined on edges we add 1.
            d=d+1+ HostTree.get_distance(node1, node2, topology_only=True)
    return d

#########################
#creates a multifutcated tree from a reconciliation
def multifurcated_tree(rec, HTree, PTree):
    mulTree=HTree.copy()
    #get the set of leaves in Htree.
    Leaves=set()
    Leaves=set(HTree.get_leaf_names(is_leaf_fn=None))
    host_leaves_already_used=list()
    for i in range(0,len(rec)):
        p,h= rec[i]
        #add only the internal nodes of the parasite
        if PTree.search_nodes(name=p):
                P=PTree.search_nodes(name=p)[0]
                #We add the pendant edge also for the leaves of P.
                #if (P.is_leaf()==False):
                #get the host node
                if not (p == 'P0'):
                    if mulTree.search_nodes(name=h):
                        H=mulTree.search_nodes(name=h)[0]
                        #p is mapped into an internal node h of the host we add a pendant edge h->p
                        if (H.is_leaf()==False):
                                H.add_child(name=p)
                        #p is mapped into a leaf of the host, then we add the new vertex to the father
                        elif (H.is_leaf()==True):
                            #we add p to the parent of the leaf h
                            #H.add_sister(name=p)
                            #get the parent of H
                            par=H.up
                            #if I have already mapped something to this leaf.
                            if (h in host_leaves_already_used):
                                par.add_child(name=p)
                            #print par
                            #create a new subtree
                            else:
                                t=Tree()
                                t.add_child(name=p)
                                t.add_child(name=h)
                                #print t
                                par.add_child(t)
                                H.delete()
                                host_leaves_already_used.append(h)
                        else:
                            print "error"
#print mulTree
    return mulTree

#########################
#define the Triplet_distance for multifurcated rooted trees. The implementation is trivial
def Triplet_Distance(tree1, tree2):
    dist=0
    #get the set of leaves in tree1. They are the same as tree2.
    Leaves=set()
    Leaves=tree1.get_leaf_names(is_leaf_fn=None)
    #iterate over all triplets of t1
    for item in list(itertools.combinations(Leaves, 3)):
        (a,b,c)=item
        if ((a!=b) and (a!=c) and (b!=c)):
                        #check the triple type
                        x1=tree1.get_common_ancestor(a,b)
                        y1=tree1.get_common_ancestor(a,c)
                        z1=tree1.get_common_ancestor(c,b)
                        x2=tree2.get_common_ancestor(a,b)
                        y2=tree2.get_common_ancestor(a,c)
                        z2=tree2.get_common_ancestor(c,b)
                        #determine one of the 4 topologies for each triple
                        if ((x1==y1) and (y1==z1)): triple1=0
                        elif ((x1==y1) and (y1 != z1)): triple1=1
                        elif ((x1==z1) and (y1 != z1)): triple1=2
                        elif ((z1==y1) and (y1 != x1)): triple1=3
                        # determine the triple topology for the tree2
                        if ((x2==y2) and (y2==z2)): triple2=0
                        elif ((x2==y2) and (y2 != z2)): triple2=1
                        elif ((x2==z2) and (y2 != z2)): triple2=2
                        elif ((z2==y2) and (y2 != x2)): triple2=3
                        # determine the distance. We consider ONLY resolved triplets.
                        if (triple1!=triple2):
                            if (triple1==0) or (triple2==0):
                                dist=dist+1
                            else:
                                dist=dist+2
    if (dist==0):
        print "Triplet distance 0"
    return dist



#########################
#Defines the clusters for each internal node of the tree
def tree_clusters(t):
    
    clusters_set=set()
    for node in t.traverse("postorder"):
        if (node.is_leaf()!= True):
            clusters_set.add(frozenset(node.get_leaf_names()))
    #print clusters_set

    return clusters_set

#define the RF for multifurcated rooted trees. It coincides with the cluster distance
def RF(tree1, tree2):
    
    #have to turn lists into tuples as lists are not hashable
    
    first_set = tree_clusters(tree1)
    secnd_set = tree_clusters(tree2)
    
    #diff=first_set.symmetric_difference(secnd_set)
    diff=first_set ^ secnd_set
    #print tree1
    #print first_set
    #print tree2
    #print secnd_set
    #print "diff"
    #print diff
    #control if we get two identical trees for two different reconciliaitons. This should not happen
    if (len(diff) < 1.0000):
        print diff
        print tree1
        print tree2
    
    return len(diff)/2

##############e#######################################################
#                                 MAIN                               #
######################################################################

argc = len(argv) - 1

if ((argc != 2) or (int(argv[2]) > 5)):
    print argc
    usage()
    exit(0)

#measures the time of the program
start = time.time()

if argc == 2:
    script, reconciliations_file, distance_chosen = argv

#reconciliations= open(reconciliations_file, "w")
all_distances_file=reconciliations_file+ "_distances.txt"
#we dont' need to write all the multifurctaed trees
multifurcated_file=reconciliations_file+"_multifurcatedTrees.txt"
multifurcated_trees_file=open(multifurcated_file, "w")


#This requires the file of the reconciliaitons to be in the format of the output of EUCALYPT
with open(reconciliations_file) as reconciliations:

    line=reconciliations.readline()
    #read the host tree
    line=reconciliations.readline()
    line = line.rstrip('\n')
    records=line.split("=",1)[1]
    records+=(';')
    Host_tree=Tree(records,format=8)
    
    #read the parasite tree
    line=reconciliations.readline()
    line = line.rstrip('\n')
    records=line.split("=",1)[1]
    records+=(';')
    Parasite_tree=Tree(records,format=8)
    
    ######################################################################
    #            Parameters to keep                                      #
    ######################################################################
    
    #read the number of leaves in the host tree
    line=reconciliations.readline()
    line = line.rstrip('\n')
    records=line.split('=')
    #nr of vertices, internal vertices and leaves of the host
    nr_vertices_host=float(records[1])
    nr_leaves_host=(nr_vertices_host+1)/2
    nr_internal_vertices_host=nr_vertices_host - nr_leaves_host
    #read the number of leaves in the parasite tree
    line=reconciliations.readline()
    line = line.rstrip('\n')
    records=line.split('=')
    nr_vertices_parasite=float(records[1])
    nr_leaves_parasite=(nr_vertices_parasite+1)/2
    nr_internal_vertices_parasite=nr_vertices_parasite - nr_leaves_parasite
    #normalization of the distances
    normalize_path=nr_internal_vertices_parasite* diameter(Host_tree)
    normalize_discrete=2* nr_internal_vertices_parasite
    normalize_triple= 2*(scipy.special.binom((nr_vertices_parasite+nr_leaves_host),3) - scipy.special.binom(nr_leaves_host,3))
    normalize_RF= nr_vertices_host -1
    
    #reads the reconciliations
    #attention: the newick format doesn't include a name for the root, so as the mapping root to root is forced  we do not consider the mapping of the root of the parasite P0
    #Keeps the list of all reconciliaiotns as a list of sets of pairs (p,h) in Reconciliaiton_Pairs_List
    Reconciliation_Pairs_List=list()
    Multifurcated_Trees_List=list()
    for line in reconciliations:
        li=line.strip()
        if not (li.startswith("#") or li.startswith("[")):
            line = line.rstrip('\n')
            records = line.split(', ') # o records = line.split('\t')
            reconciliation=list()
            for r in records:
                #remove the space before the parasite
                p , h = r.split('@')
                #add the pair (p,h) to the list
                reconciliation.append((p,h))
            #sort the reconcilaition by parasite node label
            reconciliation.sort(key=lambda pair: pair[0])

            #append the reconcilaiotn to the list of reconciliaiotns
            Reconciliation_Pairs_List.append(reconciliation)
            #append the multifurcated tree to the list
            multifurcated= multifurcated_tree(reconciliation, Host_tree, Parasite_tree)
            Multifurcated_Trees_List.append(multifurcated)
            #pdb.set_trace()
            # mulTree.write(format=8,output=outputFile)
            multifurcated_trees_file.write(multifurcated.write(format=9))
            multifurcated_trees_file.write('\n')

#closing reconciliation file
reconciliations.close()
multifurcated_trees_file.close()
    
##############e#######################################################
#            Calculation of Different metrics                        #
#            - discreate metric                                      #
#            - path metric                                           #
#            - triplet based metric (multifurcated trees)            #
#            - RF metric (multifurcated trees)                       #
######################################################################

#n is the number of reconciliaitons we will calculate (n^2-n)/2
n=len(Reconciliation_Pairs_List)

if (n<2): sys.exit(0)


##writing the distances in two different files:
discreteAndpath_distances_file=reconciliations_file+ "_Disc_path_distances.txt"
minRiconcil_path_distance_file=reconciliations_file+ "_minRicons_path_distances.txt"


distances_file = open(all_distances_file, "w")
#minRiconcil_path_file=open(minRiconcil_path_distance_file,"w")

#discreteAndpath_distances_file=open(discreteAndpath_distances_file, "w")



#initialization of the distances
discrete_distance=list()
path_distance=list()
triplet_distance=list()
RF_distance=list()

#minRicons_path_list=list()
#min_path_dist=1.0000


#########################
#TODO if n >5000 then I should randomly pick 5000 reconciliations

###To uncomment and finish the writing of the distances in this case. it will choose the distances random
if (n<=0):
    #nr_pairs=(100*(100-1)/2)
    #for N in range(nr_pairs):
    #   i=random.randint(0,n-1)
    #   j=random.randint(0,n-1)
    #   #calculate the discrete metric
    #   dist=float(((len(set(Reconciliation_Pairs_List[i]).symmetric_difference(Reconciliation_Pairs_List[j])))) / float(2 * nr_vertices_parasite)/(2 * nr_vertices_parasite))
    #   discrete_distance.append(dist)
        #calculate the path metric
        #   dist=float(path_metric(Reconciliation_Pairs_List[i], Reconciliation_Pairs_List[j], Host_tree)/max_path_distance)
        #path_distance.append(dist)
  n=0

#########################
else:
    #starts keeping time
    start = time.time()
    for i in range(0,n):
        for j in range(i,n):
            if (i != j):
                #calculate the discrete metric
                dist=float(((len(set(Reconciliation_Pairs_List[i]).symmetric_difference(Reconciliation_Pairs_List[j])))) / float(normalize_discrete))
                if ((dist>1) or (dist==0)):
                    print " discrete: \n"
                    print dist
                    print Multifurcated_Trees_List[i]
                    print Multifurcated_Trees_List[j]
                    print Reconciliation_Pairs_List[i]
                    print Reconciliation_Pairs_List[j]
                discrete_distance.append(dist)
                #discreteAndpath_distances_file.write('%f\t'%dist)
                #calculate the path metric
                dist=float(path_metric(Reconciliation_Pairs_List[i], Reconciliation_Pairs_List[j], Host_tree))/float(normalize_path)
                if ((dist>1) or (dist==0)):
                    print " path: \n"
                    print dist
                    print Multifurcated_Trees_List[i]
                    print Multifurcated_Trees_List[j]
                    print Reconciliation_Pairs_List[i]
                    print Reconciliation_Pairs_List[j]
                path_distance.append(dist)
                dist= float(RF(Multifurcated_Trees_List[i],Multifurcated_Trees_List[j])) / float(normalize_RF)
                if ((dist>1) or (dist==0)):
                    print " RF: \n"
                    print dist
                    print nr_vertices_host
                    print nr_internal_vertices_host
                    print Multifurcated_Trees_List[i]
                    print Multifurcated_Trees_List[j]
                    print Reconciliation_Pairs_List[i]
                    print Reconciliation_Pairs_List[j]
                RF_distance.append(dist)
                dist= float(Triplet_Distance(Multifurcated_Trees_List[i],Multifurcated_Trees_List[j])) / float(normalize_triple)
                if ((dist>1) or (dist==0)):
                    print " Triplet: \n"
                    print dist
                    print Multifurcated_Trees_List[i]
                    print Multifurcated_Trees_List[j]
                    print Reconciliation_Pairs_List[i]
                    print Reconciliation_Pairs_List[j]
                triplet_distance.append(dist)
        print i



#######################
#   RF- distance      #
#######################
#for i in range(0,len(Multifurcated_Trees_List)):
#   for j in range(i, len(Multifurcated_Trees_List)):
#        if (i!=j):
#           dist= RF(Multifurcated_Trees_List[i],Multifurcated_Trees_List[j])
#           RF_distance.append(dist)

#######################
#  Triplet distance   #
#######################
#for i in range(0,len(Multifurcated_Trees_List)):
#   for j in range(i, len(Multifurcated_Trees_List)):
#       if (i!=j):
#           dist= Triplet_Distance(Multifurcated_Trees_List[i],Multifurcated_Trees_List[j])
#           triplet_distance.append(dist)

#########################
#write the distances in the file  distanceFile
distances_file.write('Discrete_Dist\tPath\tTriplet\tRF\n')

#print len(discrete_distance)

for i in range(0,len(discrete_distance)):
    distances_file.write('%f\t'%discrete_distance[i])
    distances_file.write('%f\t'%path_distance[i])
    distances_file.write('%f\t'%triplet_distance[i])
    distances_file.write('%f\n'%RF_distance[i])

#########################
#write the distances in a square table for the multiscaling
#one file for each distance

discrete_multiscaling=reconciliations_file+"_discrete_multiscaling.txt"
discrete_multiscaling_file=open(discrete_multiscaling, "w")

path_multiscaling=reconciliations_file+"_path_multiscaling.txt"
path_multiscaling_file=open(path_multiscaling, "w")

triplet_multiscaling=reconciliations_file+"_triplet_multiscaling.txt"
triplet_multiscaling_file=open(triplet_multiscaling, "w")

RF_multiscaling=reconciliations_file+"_RF_multiscaling.txt"
RF_multiscaling_file=open(RF_multiscaling, "w")

#Should be choosen 1 and repeted for all the distances
#put the names of the reconciliations
discrete_multiscaling_file.write('X\t')
path_multiscaling_file.write('X\t')
triplet_multiscaling_file.write('X\t')
RF_multiscaling_file.write('X\t')

for i in range(0,n):
    if (i!=n-1):
        discrete_multiscaling_file.write('rec%d\t'%i)
        path_multiscaling_file.write('rec%d\t'%i)
        triplet_multiscaling_file.write('rec%d\t'%i)
        RF_multiscaling_file.write('rec%d\t'%i)

    else:
        discrete_multiscaling_file.write('rec%d\n'%i)
        path_multiscaling_file.write('rec%d\n'%i)
        triplet_multiscaling_file.write('rec%d\n'%i)
        RF_multiscaling_file.write('rec%d\n'%i)

#recall that for n reconciliaitons there are n(n-1)/2 distances
#recall that n is the number of reconciliaiotns
#to change the calculation of indices.
for i in range(0,n):
    discrete_multiscaling_file.write('rec%d\t'%i)
    path_multiscaling_file.write('rec%d\t'%i)
    triplet_multiscaling_file.write('rec%d\t'%i)
    RF_multiscaling_file.write('rec%d\t'%i)
    for j in range(0,n):
        if (i==j):
            discrete_multiscaling_file.write('%f'%0)
            path_multiscaling_file.write('%f'%0)
            triplet_multiscaling_file.write('%f'%0)
            RF_multiscaling_file.write('%f'%0)
        elif (i<j):
            pos=(((j+1)*j/2) - ((j-i)*(j-i-1)/2) - 1 + (i)*(n-j-1))
            discrete_multiscaling_file.write('%f'%(float(discrete_distance[pos])))
            path_multiscaling_file.write('%f'%(float(path_distance[pos])))
            triplet_multiscaling_file.write('%f'%(float(triplet_distance[pos])))
            RF_multiscaling_file.write('%f'%(float(RF_distance[pos])))
        elif (i>j):
            pos=(((i+1)*i/2) - ((i-j)*(i-j-1)/2) - 1+ (j)*(n-i-1))
            discrete_multiscaling_file.write('%f'%(float(discrete_distance[pos])))
            path_multiscaling_file.write('%f'%(float(path_distance[pos])))
            triplet_multiscaling_file.write('%f'%(float(triplet_distance[pos])))
            RF_multiscaling_file.write('%f'%(float(RF_distance[pos])))
        if (j==(n-1)):
            discrete_multiscaling_file.write('\n')
            path_multiscaling_file.write('\n')
            triplet_multiscaling_file.write('\n')
            RF_multiscaling_file.write('\n')
        else:
            discrete_multiscaling_file.write('\t')
            path_multiscaling_file.write('\t')
            triplet_multiscaling_file.write('\t')
            RF_multiscaling_file.write('\t')



#closing files
distances_file.close()
discrete_multiscaling_file.close()
path_multiscaling_file.close()
triplet_multiscaling_file.close()
RF_multiscaling_file.close()


#

#print the time
end = time.time()

print end - start

