import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import sympy as sp
from math import comb
from tqdm import tqdm
from collections import defaultdict
import math
import time
import itertools
import random
import os

def largerBin_old(num): #Given a decimal integer, returns a supersetof all possible states that it can reach in one iteration
    nums = []
    binary = d2b(num)
    zeroIndices = []
    for i in range(len(binary)):
        if binary[i] == "0":
            zeroIndices.append(n-i-1) 
    for i in range(len(zeroIndices)+1):
        for subset in itertools.combinations(zeroIndices, i):
            temp = num
            for j in subset:
                temp += 2**j
            nums.append(temp)
    return sorted(nums)

def largerBin(num, graphs): #This one checks the edges for possibilities first. Gives all possible states that it can reach in one iteration
    nums = []
    binary = d2b(num)
    blueIndices = []
    whiteIndices = []
    possibleIndices = []
    for i in range(len(binary)):
        if binary[i] == "0":
            whiteIndices.append(n-i-1)
        else:
            blueIndices.append(n-i-1)
    for i in whiteIndices:
        if len(set(blueNeighbors(graphs[num], i)).intersection(set(blueIndices))) > 0:
            possibleIndices.append(i)
    for i in range(len(possibleIndices)+1):
        for subset in itertools.combinations(possibleIndices, i):
            temp = num
            for j in subset:
                temp+= 2**j
            nums.append(temp)
    return sorted(nums)

def largerBin_directed_newRule(num, graphs):
    nums = []
    binary = d2b(num)
    blueIndices = []
    whiteIndices = []
    possibleIndices = []
    for i in range(len(binary)):
        if binary[i] == "0":
            whiteIndices.append(n-i-1)
        else:
            blueIndices.append(n-i-1)
    for i in blueIndices:
        for j in graphs[num].successors(i):
            if j not in blueIndices and j not in possibleIndices:
                possibleIndices.append(j)
    for i in range(len(possibleIndices)+1):
        for subset in itertools.combinations(possibleIndices, i):
            temp = num
            for j in subset:
                temp+= 2**j
            nums.append(temp)        
    return sorted(nums)

def sameOnes(num): #Given a decimal integer, returns all other numbers with the same number of 1 bits. 
    ones = bin(num).count("1")
    ans = []
    for i in range(2**n):
        if bin(i).count("1") == ones:
            ans.append(i)
    return ans

def sameOnes_larger(num):
    ones = bin(num).count("1")
    ans = []
    for i in range(num, 2**n):
        if bin(i).count("1") == ones:
            ans.append(i)
    return ans

def extract(lst):
    return [item[0] for item in lst]

def d2b(i): #Function to return the string of the binary representation of an integer with n digits. 
    binaryString = str(bin(i))[2:]
    while(len(binaryString) < n):
        binaryString = "0" + binaryString
    return binaryString

def b2d(bin): #Function to return the decimal represenation of a binary number
    dec = 0
    for i in range(len(bin)):
        dec += int((bin[n-i-1]))*(2**i)
    return dec

def transitionPossible(i,j): #Function to return whether graph i (based on the binary rep. of i) can transition to j
    graphi = d2b(i)
    graphj = d2b(j)
    for k in range(len(graphi)):
        if graphi[k] == "1" and graphj[k] == "0":
            return False
    return True

def differingDigits(i,j): #Given 2 integers, function returns the digits/nodes that need to be forced
    num1Rev = d2b(i)[::-1]
    num2Rev = d2b(j)[::-1]
    output = []
    for i in range(len(num1Rev)):
        if num1Rev[i] != num2Rev[i]:
            output.append(i)
    return output

def sameZeros(i,j): #Given 2 integers,returns the indices of digits that are zero in both
    num1Rev = d2b(i)[::-1]
    num2Rev = d2b(j)[::-1]
    output = []
    for i in range(len(num1Rev)):
        if num1Rev[i] == num2Rev[i] and num1Rev[i] == "0":
            output.append(i)
    return output

def blueNeighbors(graph, node): #Given a node and a graph, returns a list of its neighbors that are blue
    output = []
    for i in graph.neighbors(node):
        if graph.nodes[i]['color'] == 'blue':
            output.append(i)
    return output

def blueSuppliers(graph, node):
    output = []
    for i in list(graph.predecessors(node)):
        if graph.nodes[i]['color'] == 'blue':
            output.append(i)
    return output

def numBlueNeighbors(graph): #Given a graph, returns an array of how many blue neighbors each node has
    output = []
    for i in range(n):
        output.append(len(blueNeighbors(graph, i)))
    return output

def numBlueSuppliers(graph):
    output = []
    for i in range(n):
        output.append(len(blueSuppliers(graph, i)))
    return output

def forcedProb(graph, node, numBlueNeighbors):
    temp = 1
    for i in blueNeighbors(graph, node):
        deg = graph.degree[i]
        numBlueNeighborsi = numBlueNeighbors[i]
        iforcesnodeprob = (numBlueNeighborsi+1)/deg
        idoesntforcenodeprob = 1 - iforcesnodeprob
        temp *= idoesntforcenodeprob
    return (1-temp)

def forcedProb_Directed(graph, node):
    suppliers = list(graph.predecessors(node))
    if len(suppliers) == 0:
        return 0
    isWeighted = True
    for edge in graph.in_edges(node):
        if 'weight' not in edge:
            isWeighted = False 
    if not isWeighted:
        return len(blueSuppliers(graph, node))/len(suppliers)
    else:
        weightSum = 0
        blueWeightSum = 0
        for i in suppliers:
            weight = graph[i][node].get('weight')
            weightSum += weight
            if graph[i].get('color') == 'blue':
                blueWeightSum += weight
        return weightSum/blueWeightSum

def numswithbitcount(max, ones):
    nums = []
    for i in range(max):
        if d2b(i).count("1") == ones:
            nums.append(i)
    return nums

def print_graph(graph):
    nx.draw(graph, with_labels=True, node_color='skyblue', edge_color='gray', node_size=1500, font_size=12)
    plt.show()

def save_graph_radius_ver(graph, radius, node, time, name):
    directory = "/Users/noah01px2019/Desktop/graphs"
    plt.figure()
    plt.title("Radius = " + str(radius) + " Starting node = " + str(node) + " Propogation time = " +str(time))
    file_path = os.path.join(directory, name)
    nx.draw(graph, with_labels=True)
    plt.savefig(file_path)
    plt.close()

def graph_gen(targetGraph):
    #start = time.process_time()
    n = targetGraph.number_of_nodes()
    numGraphs = 2**n
    graphs = [None]*(numGraphs)
    for i in range(numGraphs):
        binaryString = d2b(i)
        graph = targetGraph.copy()
        for j in range(n):
            if binaryString[n-j-1] == "1":
                graph.nodes[j]['color'] = "blue"
            else:
                graph.nodes[j]['color'] = "white"
        graphs[i] = graph
    #print("Time needed for graph creation = " + str(time.process_time() - start) + " seconds")
    return graphs

def tm_generation_old(graphs): #outdated, quite slow
    start = time.process_time()
    tm = np.zeros([2**n,2**n])
    #The states of the transition matrix is stored as follows - the nth node being blue refers to a 1 in the nth digit of the binary expression of the state.
    for i in range(2**n): #Computing the transition rate matrix
        for j in range(2**n):
            stateigraph = graphs[i]
            numBlueNeighborsi = numBlueNeighbors(stateigraph)
            transitionCalculated = False
            while transitionCalculated == False:
                if transitionPossible(i,j) == False:
                    tm[i][j] = 0
                    transitionCalculated = True
                else:  
                    tm[i][j] = 1
                    for k in differingDigits(i, j): #all differing digits must go from 0 to 1, otherwise transition impossible. First consdiering probabilities of vertices being forced
                        tm[i][j] *= forcedProb(stateigraph, k, numBlueNeighborsi)
                    for m in sameZeros(i,j): #all 0s that stay as 0s must not have been forced!)
                        tm[i][j] *= (1-forcedProb(stateigraph, m, numBlueNeighborsi))
                    transitionCalculated = True
    tm[2**n-1][2**n-1] = 1
    #print("Time needed for Transition Matrix Generation = " + str(time.process_time() - start) + " seconds")
    return tm

def tm_generation(graphs):
    start = time.process_time()
    tm = np.zeros([2**n,2**n])
    #The states of the transition matrix is stored as follows - the nth node being blue refers to a 1 in the nth digit of the binary expression of the state.
    for i in range(2**n): #Computing the transition rate matrix, tqdm here
        for j in largerBin(i, graphs):
            stateigraph = graphs[i]
            numBlueNeighborsi = numBlueNeighbors(stateigraph)
            tm[i][j] = 1
            for k in differingDigits(i, j): #all differing digits must go from 0 to 1, otherwise transition impossible. First consdiering probabilities of vertices being forced
                tm[i][j] *= forcedProb(stateigraph, k, numBlueNeighborsi)
            for m in sameZeros(i,j): #all 0s that stay as 0s must not have been forced!)
                tm[i][j] *= (1-forcedProb(stateigraph, m, numBlueNeighborsi))
    tm[2**n-1][2**n-1] = 1
    #print("Time needed for Transition Matrix Generation = " + str(time.process_time() - start) + " seconds")
    return tm

def tm_generation_directed(graphs):
    start = time.process_time()
    tm = np.zeros([2**n,2**n])
    for i in range(2**n):
        for j in largerBin_directed_newRule(i, graphs):
            stateigraph = graphs[i]
            tm[i][j] = 1
            for k in differingDigits(i,j):
                tm[i][j] *= forcedProb_Directed(stateigraph, k)
            for m in sameZeros(i,j):
                tm[i][j] *= (1-forcedProb_Directed(stateigraph, m))
    tm[2**n-1][2**n-1] = 1
    return tm

def tm_generation_sub(graphs, startingSet):
    start = time.process_time()
    tm = np.zeros([2**n,2**n])
    states = largerBin(startingSet, graphs)
    for i in states: #tqdm here
        for j in largerBin(i, graphs):
            stateigraph = graphs[i]
            numBlueNeighborsi = numBlueNeighbors(stateigraph)
            tm[i][j] = 1
            for k in differingDigits(i, j): #all differing digits must go from 0 to 1, otherwise transition impossible. First consdiering probabilities of vertices being forced
                tm[i][j] *= forcedProb(stateigraph, k, numBlueNeighborsi)
            for m in sameZeros(i,j): #all 0s that stay as 0s must not have been forced!)
                tm[i][j] *= (1-forcedProb(stateigraph, m, numBlueNeighborsi))
    tm[2**n-1][2**n-1] = 1
    subtm = tm[np.ix_(states, states)]
    #print("Time needed for Transition Matrix Generation with submatrix of size " + str(len(states)) + " in " + str(time.process_time() - start) + " seconds")
    return subtm
"""
#sanity check
rowSums = np.zeros(2**n)
for i in range(2**n):
    sum = 0
    for j in range(2**n):
        sum += tm[i][j]
    rowSums[i] = sum
print(rowSums)
"""
def propogation_time_solver_algebraic(tm): #Turns out this one is really slow. Use the other one!
    #Solving the system of equations to find expected number of steps until absorption
    mus = sp.symbols('mu0:%d'%2**n)
    a = sp.Matrix(tm)
    eqns = []
    for i in range(2**n):
        if i != (2**n-1):
            eqns.append(a.row(i).dot(mus) + 1 - (mus[i])) 
        else:
            eqns.append(a.row(i).dot(mus))
    ans = list(sp.linsolve(eqns[-(2**n-1):], mus[-(2**n-1):])) #Has no value for the empty starting set
    return ans

def propogation_time_solver(tm):
    start = time.process_time()
    size = len(tm)
    #The fastest one
    solutions = np.zeros(size)
    solutions[size-1] = 0
    for i in range(size-2, 0, -1):
        temp = calc_expected(tm, i, solutions)
        solutions[i] = temp
    solutions[0] = float('inf')
    #print("Solving for propogation time took " + str(time.process_time() - start) + " seconds")
    return solutions

def propogation_time_solver_inverse(tm): #This one uses the linalg method outlined in Hogben and Jesse's paper. It's Slow :(
    size = len(tm)
    for row in tm:
        row[-1] -= 1
    start = time.process_time()
    temp = np.linalg.inv(tm-np.identity(size))
    #print("Finding an inverse took " + str(time.process_time() - start) + " seconds")
    return temp[0][size-1] + 1

def calc_expected(tm, row, solutions):
    size = len(tm)
    temp = 0
    coefficient = 1-tm[row][row]
    if coefficient == 0:
        return float('inf')
    for i in range(row+1, size):
        if tm[row][i] != 0:
            temp += tm[row][i] * (solutions[i]+1)
    temp+= (1-coefficient)
    temp /= coefficient
    return temp

def can_be_calculated(tm, row, solutions):
    ans = True
    for i in range(2**n-1):
        if tm[row][i] != 0:
            if solutions[i] == 0 and i != row:
                ans = False
    return ans

#generating arrays for zero forcing sets by size of intial set

def optimal_zfs_by_size(ptimes, size):
    solutions = []
    indices = numswithbitcount(2**n,size)
    min = ptimes[indices[0]]
    for i in indices:
        if ptimes[i] < min:
            min = ptimes[i]
    for i in indices:
        if ptimes[i] == min:
            solutions.append((i, min))
    return solutions #returns list with two elements. first element contains all zfs that achieve the min time, second elment is the min time

def print_zfs(targetGraph, num, graphs, transitionTimes):
    color_map = []
    for i in range(n):
        if graphs[num].nodes[i]['color'] == 'blue':
            color_map.append('blue')
        else: 
            color_map.append('white')      
    plt.title("Expected Propogation Time = " + str(transitionTimes[num]))
    nx.draw(targetGraph, node_color=color_map, with_labels=True)
    plt.show()

def save_zfs(targetGraph, num, graphs, title, transitionTimes):
    color_map = []
    for i in range(n):
        if graphs[num].nodes[i]['color'] == 'blue':
            color_map.append('blue')
        else: 
            color_map.append('white')  
    pos = nx.spring_layout(graphs[num], k = 1)
    plt.figure()    
    plt.title("Expected Propogation Time = " + str(transitionTimes[num]))
    nx.draw(targetGraph, node_color=color_map, with_labels=True)
    plt.savefig(title)

def return_automorphic_states(graphs,index,automorphicKeys):
    key = tuple(blue_degree_sequence(graphs, index)) + tuple(sorted(numBlueNeighbors(graphs[index])))
    ans = []
    for state in automorphicKeys[key]:
        if nx.is_isomorphic(graphs[index], graphs[state], node_match=color_match):
            ans.append(state)
            global succesfulChecks
            succesfulChecks += 1
        else:
            global failedChecks
            failedChecks += 1
    return ans

def automorphism_groups(graphs):
    start = time.process_time()
    automorphicKeys = pre_automorphism_checking_dictionary(graphs)
    states = range(2**n)
    groups = []
    while len(states) > 0: 
        group = return_automorphic_states(graphs,states[0], automorphicKeys)
        groups.append(group)
        states = [item for item in states if item not in group]
    print("Time needed to generate automorphism groups " + str(time.process_time() - start) + " seconds")
    return groups

def pre_automorphism_checking_dictionary(graphs): #this one creates a dictionary of all things to check before checking for automorphisms
    groups = defaultdict(list)
    for i in range(2**n):
        key = tuple(blue_degree_sequence(graphs, i)) + tuple(sorted(numBlueNeighbors(graphs[i])))
        groups[key].append(i)
    return groups

def color_match(n1, n2):
    return n1['color'] == n2['color']

def blue_degree_sequence(graphs, index): #Returns the degree sequence of the blue nodes within a graph.
    seq = []
    graph = graphs[index]
    for v in graph.nodes:
        if graph.nodes[v]['color'] == 'blue':
            seq.append(graph.degree[v])
    seq.sort(reverse=True)
    return seq

def blue_degree_sequence_groups(graphs):
    groups = defaultdict(list)
    for i in range(2**n):
        degSeq = tuple(blue_degree_sequence(graphs,i))
        groups[degSeq].append(i)
    return groups

def degree_sequence_groups(graphs):
    groups = defaultdict(list)
    for i in range(2**n):
        degseq = tuple(blue_then_white_degree_sequence(graphs,i))
        groups[degseq].append(i)
    return groups

def blue_then_white_degree_sequence(graphs, index):
    blueSeq = []
    whiteSeq = []
    graph = graphs[index]
    for v in graph.nodes:
        if graph.nodes[v]['color'] == 'blue':
            blueSeq.append(graph.degree[v])
        else:
            whiteSeq.append(graph.degree[v])
    blueSeq.sort(reverse=True)
    whiteSeq.sort(reverse=True)
    #the zero acts as the split
    return tuple(blueSeq), tuple(whiteSeq)

def white_degree_sequence(graphs, index):
    seq = []
    graph = graphs[index]
    for v in graph.nodes:
        if graph.nodes[v]['color'] == 'white':
            seq.append(graph.degree[v])
    seq.sort(reverse=True)
    return seq

def white_degree_sequence_groups(graphs):
    groups = defaultdict(list)
    for i in range(2**n):
        degSeq = tuple(white_degree_sequence(graphs,i))
        groups[degSeq].append(i)
    return groups

def generate_automorphism_matrix(automorphism_groups, graphs):
    start = time.process_time()
    numStates = len(automorphism_groups)
    matrix = np.zeros([numStates,numStates])
    matrix[0][0] = 1
    for i in tqdm(range(numStates)):
        startingState = (automorphism_groups[i])[0]
        startingStateGraph = graphs[startingState]
        possibleStates = largerBin(startingState, graphs)
        numBlueNeighborsStartingState = numBlueNeighbors(startingStateGraph)
        for j in range(i, numStates):
            prob = 0
            targetGroup = automorphism_groups[j]
            transitionPossibleSet = sorted(list(set(possibleStates).intersection(set(targetGroup))))
            for graphState in transitionPossibleSet:
                temp = 1
                for k in differingDigits(startingState, graphState): #all differing digits (nodes) must go from 0 to 1, otherwise transition impossible. First consdiering probabilities of vertices being forced
                    probability = forcedProb(startingStateGraph, k, numBlueNeighborsStartingState)
                    temp *= probability
                for m in sameZeros(startingState,graphState): #all 0s that stay as 0s must not have been forced!)
                    temp *= (1-forcedProb(startingStateGraph, m, numBlueNeighborsStartingState))
                prob += temp
            matrix[i][j] = prob
    print("Time needed to generate automorphism transition matrix " + str(time.process_time() - start) + " seconds")
    return matrix

def calculate_transition_probability(startingGraph, startingGraphIndex, endingGraphIndex):
    needToBeForced = differingDigits(startingGraphIndex, endingGraphIndex) #vertices that need to be forced
    needToStayUnforced = sameZeros(startingGraphIndex, endingGraphIndex) #vertices that need to stay unforced
    temp = 1
    for vertex in needToBeForced:
        temp =1
    return temp

def transition_times_by_maxclique(targetGraph, transitionTimes, zeroForcingSetSize): #currently breaks for complete graphs
    n = targetGraph.number_of_nodes()
    maxClique = nx.approximation.max_clique(targetGraph)
    cliqueTransitionTimeAvg = 0
    nonCliqueTransitionTimeAvg = 0
    sets = numswithbitcount(2**n, zeroForcingSetSize)
    cliqueNodeSets = list(itertools.combinations(maxClique, zeroForcingSetSize))
    cliqueSets = []
    for i in cliqueNodeSets:
        temp = 0
        for j in i:
            temp += 2**j
        cliqueSets.append(temp)
    for i in sets:
        if i in cliqueSets:
            cliqueTransitionTimeAvg += transitionTimes[i]
        else:
            nonCliqueTransitionTimeAvg += transitionTimes[i]
    cliqueTransitionTimeAvg /= len(cliqueSets)
    nonCliqueTransitionTimeAvg /= (len(sets)-len(cliqueSets))
    return cliqueTransitionTimeAvg, nonCliqueTransitionTimeAvg
    #print("The average transition time average with the starting node being in the maxclique (size " + str(len(maxClique)) + ") is " + str(cliqueTransitionTimeAvg) + " iterations")
    #print("The average transition time average with the starting node not being in the maxclique is " + str(nonCliqueTransitionTimeAvg) + " iterations")

def transition_times_by_maxclique_minimum(targetGraph, transitionTimes, zeroForcingSetSize): #currently breaks for complete graphs
    n = targetGraph.number_of_nodes()
    maxClique = nx.approximation.max_clique(targetGraph)
    cliqueTransitionTimeAvg = 0
    nonCliqueTransitionTimeAvg = 0
    sets = numswithbitcount(2**n, zeroForcingSetSize)
    cliqueNodeSets = list(itertools.combinations(maxClique, zeroForcingSetSize))
    cliqueSets = []
    for i in cliqueNodeSets:
        temp = 0
        for j in i:
            temp += 2**j
        cliqueSets.append(temp)
    min_index = min(sets, key=lambda i: transitionTimes[i])
    min_value = transitionTimes[min_index]
    if min_index in cliqueSets:
        return True
    else:
        return False
    #print("The average transition time average with the starting node being in the maxclique (size " + str(len(maxClique)) + ") is " + str(cliqueTransitionTimeAvg) + " iterations")
    #print("The average transition time average with the starting node not being in the maxclique is " + str(nonCliqueTransitionTimeAvg) + " iterations")

def return_transition_times(targetGraph):
    graphs = graph_gen(targetGraph)
    if nx.is_directed(targetGraph):
        transitionMatrix = tm_generation_directed(graphs)
    else:
        transitionMatrix = tm_generation(graphs)
    transitionTimes = propogation_time_solver(transitionMatrix)
    return transitionTimes

def has_multiple_zero_in_degree_nodes(G):
    zero_in_degree_count = sum(1 for node in G.nodes if G.in_degree(node) == 0)
    return zero_in_degree_count > 1

def has_cycle_without_node(graph, node):
    allcycles = nx.simple_cycles(graph)
    for i in allcycles:
        if node not in i:
            return True
    return False

def length_of_all_longest_paths(graph, startNode):
    pathLengths = [0]
    for i in [node for node in graph if node != startNode]:
        if nx.has_path(graph, startNode, i):
            paths = nx.all_simple_paths(graph, startNode, i)
            max_length = max(len(path) for path in paths)
            pathLengths.append(max_length -1 )
    return pathLengths

def all_nodes_reachable(graph, startNode):
    for i in [node for node in graph if node != startNode]:
        if not nx.has_path(graph, startNode, i):
            return False
    return True

def create_tournament(size):
    G = nx.DiGraph()
    G.add_nodes_from([0,size-1])
    for i in range(size-1):
        for j in range(i+1, size):
            G.add_edge(i,j)
    return(G)

def return_nonzeros(tm):
    output = []
    for i in range(len(tm)):
        temp = [(i)]
        for j in range(len(tm)):
            if tm[i][j] != 0:
                temp.append((tm[i][j], j))
        output.append(temp)
    return output

def tournament_erzf():
    erzf_dic = {}
    return tournament_erzf_state(1, erzf_dic), erzf_dic

def transition_prob_tournament_erzf(startingState, endingState):
    prob = 1
    start = d2b(startingState)
    end = d2b(endingState)
    #for i in range(len(start)-2,-1,-1):
    return 0

def tournament_erzf_state(state, expectedTimes):
    if state in expectedTimes:
        return expectedTimes[state]
    if state == 2**n-1:
        expectedTimes[state] = 0
        return 0
    else:
        et = 0
        for possibleState in largerBin_old(state)[1:]:
            start = d2b(state)
            end = d2b(possibleState)
            prob = 1
            numOnes = 1
            for i in range(n-2,-1,-1):
                if start[i] == '0':
                    if end[i] == '0':
                        prob *= (1-numOnes/(n-i-1))
                    else:
                        prob *= numOnes/(n-i-1)
                else: 
                    numOnes += 1
            et += prob * (tournament_erzf_state(possibleState, expectedTimes) +1)
    expectedTimes[state] = et
    return et

def tournament_probability_matrix_method_erzf():
    mat = np.zeros((n,n)) #rows = timestep starting at 0, columns the nodes
    for i in range(n):
        for j in range(i+1):
            mat[i][j] = 1
    for i in range(1,n):
        for j in range(i+1, n):
            temp = 0
            for k in range(j):
                temp += mat[i-1][k]
            mat[i][j] = mat[i-1][j] + (1-mat[i-1][j])*temp/(j)

    forcedProbVector = [0]
    cumForcedProbVector = [0]
    for i in range(1,n):
        temp = 1
        for j in range(n):
            temp *= mat[i][j]
        forcedProbVector.append(temp-cumForcedProbVector[i-1])
        cumForcedProbVector.append(temp)

    erzf = 0
    for i in range(n):
        erzf += i*forcedProbVector[i]
    print(sum(forcedProbVector))
    return mat

def longest_string_matrix():
    mat = np.zeros((n,n))
    mat[0][0] = 1
    for i in range(1,n):
        for j in range(i,n-1):
            prob = 1
            for k in range(i,j+1):
                prob *= i/k
            prob *= 1-(i)/(j+1)
            mat[i][j] = prob
    for i in range(n):
        mat[i][n-1] = 1-sum(mat[i])
    return mat

def create_bidirectional_path_graph():
    path_graph = nx.path_graph(n)
    directed_graph = nx.DiGraph()
    for u, v in path_graph.edges():
        directed_graph.add_edge(u, v) 
        directed_graph.add_edge(v, u) 
    return directed_graph

def create_bidirectional_complete_graph():
    complete_graph = nx.complete_graph(n)
    directed_graph = nx.DiGraph()
    for u, v in complete_graph.edges():
        directed_graph.add_edge(u, v) 
        directed_graph.add_edge(v, u) 
    return directed_graph

def create_bidirectional_complete_bipartite_graph(x,y):
    G = nx.complete_bipartite_graph(x, y)
    directed_graph = nx.DiGraph()
    for u, v in G.edges():
        directed_graph.add_edge(u, v)
        directed_graph.add_edge(v, u)
    return directed_graph

def create_bidirectional_cycle_graph():
    cycle_graph = nx.cycle_graph(n)
    directed_graph = nx.DiGraph()
    for u, v in cycle_graph.edges():
        directed_graph.add_edge(u, v) 
        directed_graph.add_edge(v, u) 
    return directed_graph

def create_bidirectional_sun_graph(): #n nodes, n/2 in the cycle, n/2 pendant nodes, cycle nodes labelled 0 to n/2-1
    cycleLength = int(n/2)
    G = nx.cycle_graph(cycleLength)
    for i in range(cycleLength):
        pendant_node = cycleLength + i  
        G.add_node(pendant_node)
        G.add_edge(i, pendant_node)
    directed_graph = nx.DiGraph()
    for u, v in G.edges():
        directed_graph.add_edge(u, v) 
        directed_graph.add_edge(v, u) 
    return directed_graph

def create_bidirectional_grid_graph(x, y):
    grid_graph = nx.grid_2d_graph(x, y)
    directed_graph = nx.DiGraph()
    for u, v in grid_graph.edges():
        directed_graph.add_edge(u, v)
        directed_graph.add_edge(v, u)
    mapping = { (i, j): i * y + j for i in range(x) for j in range(y) }
    relabeled_graph = nx.relabel_nodes(directed_graph, mapping)
    return relabeled_graph

def create_bidirectional_hypercube_graph():
    d = int(math.log2(n))
    hypercube_graph = nx.hypercube_graph(d)
    mapping = {node: i for i, node in enumerate(hypercube_graph.nodes())}
    hypercube_graph = nx.relabel_nodes(hypercube_graph, mapping)
    directed_graph = nx.DiGraph()
    for u, v in hypercube_graph.edges():
        directed_graph.add_edge(u, v)
        directed_graph.add_edge(v, u)
    return directed_graph

def all_nodes_blue(graph):
    return all(data.get('color') == 'blue' for _, data in graph.nodes(data=True))

def rzf_simulation(targetGraph, startingSet):
    iterations = 0
    start = 0
    for num in startingSet:
        start += 2**num
    binaryString = d2b(start)
    graph = targetGraph.copy()
    for j in range(n):
        if binaryString[n-j-1] == "1":
            graph.nodes[j]['color'] = "blue"
        else:
            graph.nodes[j]['color'] = "white"
    while not all_nodes_blue(graph):
        iterations += 1
        turnedBlue = []
        for i in range(n):
            in_degree = graph.in_degree(i)
            if in_degree > 0:
                blue_predecessors = len(blueNeighbors(graph,i))
                probability = blue_predecessors / in_degree
                if random.random() <= probability:
                    turnedBlue.append(i)
        for node in turnedBlue:
            graph.nodes[node]['color'] = 'blue'
    return iterations

n = 2**16
G = create_bidirectional_hypercube_graph()
"""start = time.process_time()
tt = return_transition_times(G)
zfs_size_one_indices = [i for i in range(2**n) if bin(i).count('1') == 1]
print(tt[zfs_size_one_indices])
print("Time needed to generate transition times = " + str(time.process_time() - start) + " seconds")"""
sum = 0
zfsSet = [0]
start = time.process_time()
for i in range(1000):
    sum += rzf_simulation(G, zfsSet)
print(sum)
print("Time needed to run 1000 simulations = " + str(time.process_time() - start) + " seconds")

#for simulation, track blue nodes and only consider successors to blue nodes