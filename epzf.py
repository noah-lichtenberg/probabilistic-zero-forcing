import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import sympy as sp
from math import comb
from tqdm import tqdm
from collections import defaultdict
import time
import itertools
import random

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

def numBlueNeighbors(graph): #Given a graph, returns an array of how many blue neighbors each node has
    output = []
    for i in range(n):
        output.append(len(blueNeighbors(graph, i)))
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

def numswithbitcount(max, ones):
    nums = []
    for i in range(max):
        if d2b(i).count("1") == ones:
            nums.append(i)
    return nums

def print_graph(graph):
    nx.draw(graph)
    plt.show()

def graph_gen(targetGraph):
    start = time.process_time()
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
        for j in largerBin(i, graphs):
            temp = 0 #KEEP WRITING HERE!!
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
    for i in range(size-2, -1, -1):
        temp = calc_expected(tm, i, solutions)
        solutions[i] = temp
    #print("Solving for propogation time took " + str(time.process_time() - start) + " seconds")
    return solutions

def propogation_time_solver_4(tm): #This one uses the linalg method outlined in Hogben and Jesse's paper. It's Slow :(
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
    for i in range(row+1, size):
        temp+= tm[row][i] * (solutions[i]+1)
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

def print_zfs(num, transitionTimes):
    color_map = []
    for i in range(n):
        if graphs[num].nodes[i]['color'] == 'blue':
            color_map.append('blue')
        else: 
            color_map.append('white')      
    plt.title("Expected Propogation Time = " + str(transitionTimes[num]))
    nx.draw(targetGraph, node_color=color_map, with_labels=True)
    plt.show()

def save_zfs(num, title, transitionTimes):
    color_map = []
    for i in range(n):
        if graphs[num].nodes[i]['color'] == 'blue':
            color_map.append('blue')
        else: 
            color_map.append('white')  
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
    transitionMatrix = tm_generation(graphs)
    transitionTimes = propogation_time_solver(transitionMatrix)
    return transitionTimes


target_Graph = nx.barbell_graph(4,3)
n = target_Graph.number_of_nodes()
print(return_transition_times(target_Graph))












"""
targetGraph = nx.barbell_graph(5,4) #graph is created here.   #cycle 14 takes 1k seconds for auto groups and 1 to 2 secs for trans matrix
n = targetGraph.number_of_nodes()
zeroForcingSetSize = 2 #size of the zero forcing set
graphs = graph_gen(targetGraph) #array of graphs, contains every single possible state
blueDegSeqGroups = blue_degree_sequence_groups(graphs) #dictionary of the degree sequences of blue nodes for every possible state
transitionMatrix = tm_generation(graphs)
transitionTimes = propogation_time_solver(transitionMatrix)
transition_times_by_maxclique(targetGraph, transitionTimes, zeroForcingSetSize)
"""







#trying random graphs

attempts = 0
cliqueLarger = 0
nonCliqueLarger = 0
edgelist = []
maxcliquelist = []
numNodes = 10
info = np.zeros((numNodes, 2))
while attempts < 1000:
    nodes = numNodes
    edges = random.randint(nodes-1, nodes*(nodes-1)/2 -1 ) #i think this is biased towards larger graphs
    targetGraph = nx.gnm_random_graph(nodes, edges)
    while (not nx.is_connected(targetGraph)):
        edges = random.randint(nodes -1 , nodes*(nodes-1)/2 -1 ) #i think this is biased towards larger graphs
        targetGraph = nx.gnm_random_graph(nodes, edges)
    edgelist.append(edges)
    maxCliqueSize = len(nx.approximation.max_clique(targetGraph))
    maxcliquelist.append(maxCliqueSize)
    #print_graph(targetGraph)
    n = nodes
    graphs = graph_gen(targetGraph)
    blueDegSeqGroups = blue_degree_sequence_groups(graphs) #dictionary of the degree sequences of blue nodes for every possible state
    transitionMatrix = tm_generation(graphs)
    transitionTimes = propogation_time_solver(transitionMatrix)
    info[maxCliqueSize][0] += 1
    info[maxCliqueSize][1] += transition_times_by_maxclique_minimum(targetGraph, transitionTimes, 1)
    #cliqueLarger += transition_times_by_maxclique_minimum(targetGraph, transitionTimes, 1)
    attempts += 1

filtered_data = info[2:]
x = np.arange(2, len(filtered_data) + 2)
y1 = filtered_data[:, 1] / filtered_data[:, 0]
y2 = x / numNodes
plt.figure(figsize=(8, 6))
plt.plot(x, y1, 'bo-', label='Percent of Graphs where minzeroforcingnode is in Maxclique')
plt.plot(x, y2, 'r--', label='Theoretical percentage if maxclique doesnt affect zeroforcingtime')
plt.xlabel('Row Number (Starting from 2)')
plt.ylabel('Percentage Chance ')
plt.legend()
plt.show()
#print(nonCliqueLarger) #avg time so far for most graphs lean towards not starting in clique as faster



"""
#checking efficiency of automorphism algorithm
succesfulChecks = 0
failedChecks = 0
automorphismGroups = automorphism_groups(graphs)[1:]
automorphismMatrix = generate_automorphism_matrix(automorphismGroups, graphs)
automorphicTransitionTimes = propogation_time_solver(automorphismMatrix)
print(automorphicTransitionTimes)
print(succesfulChecks)
print(failedChecks)
"""

#Compare values in automorphism transition times and regular transition times


"""
start = time.process_time()
for i in automorphism_groups(graphs):
    print(i)
    print(transitionTimes[i[0]])
print("Time needed for algebraic epzf solving = " + str(time.process_time() - start) + " seconds")
"""

"""
optimalSol = optimal_zfs_by_size(transitionTimesFixed, zeroForcingSetSize)
print_zfs((optimalSol[0])[0])
"""

#next steps, symmetric graphs, cut edges