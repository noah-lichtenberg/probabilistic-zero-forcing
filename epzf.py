import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import sympy as sp
from math import comb
from tqdm import tqdm
import time
import itertools

n=5 #Size of the graph. Many functions depend on this.
targetGraph = nx.path_graph(n) #graph is created here.
zeroForcingSetSize = 1 #size of the zero forcing set

def largerBin(num): #Given a decimal integer, retunrs all possible states that it can reach in one iteration
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

def extract(lst):
    return [item[0] for item in lst]

def d2b(i): #Function to return the string of the binary representation of an integer with n digits. 
    binaryString = str(bin(i))[2:]
    while(len(binaryString) < n):
        binaryString = "0" + binaryString
    return binaryString

print(largerBin(13))

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
    print("Time needed for graph creation = " + str(time.process_time() - start) + " seconds")
    return graphs

def tm_generation(graphs):
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
    print("Time needed for Transition Matrix Generation = " + str(time.process_time() - start) + " seconds")
    return tm

def tm_generation_sub(startingSetDec):
    states = largerBin(startingSetDec)
    numStates = len(states)

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

def propogation_time_solver(tm): #Turns out this one is really slow. Use the other one!
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

def propogation_time_solver_2(tm):
    #Solving the system of equations without using sympy
    solutions = np.zeros(2**n) #Rewrite to solve from bottom to top, 
    needs_to_be_calculated = list(range(1,2**n-1))
    while len(needs_to_be_calculated) > 0:
        for i in needs_to_be_calculated:
            if can_be_calculated(tm, i, solutions):
                temp = calc_expected(tm, i, solutions)
                solutions[i] = temp
                needs_to_be_calculated.remove(i)
    return solutions

def propogation_time_solver_3(tm):
    #The fastest one
    solutions = np.zeros(2**n)
    solutions[2**n-1] = 0
    for i in range(2**n-2, 0, -1):
        temp = calc_expected(tm, i, solutions)
        solutions[i] = temp
    return solutions

def calc_expected(tm, row, solutions):
    temp = 0
    coefficient = 1-tm[row][row]
    for i in range(row+1, 2**n):
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

def print_zfs(num):
    color_map = []
    for i in range(n):
        if graphs[num].nodes[i]['color'] == 'blue':
            color_map.append('blue')
        else: 
            color_map.append('white')      
    plt.title("Expected Propogation Time = " + str(transitionTimesFixed[num]))
    nx.draw(targetGraph, node_color=color_map, with_labels=True)
    plt.show()

def save_zfs(num, title):
    color_map = []
    for i in range(n):
        if graphs[num].nodes[i]['color'] == 'blue':
            color_map.append('blue')
        else: 
            color_map.append('white')  
    plt.figure()    
    plt.title("Expected Propogation Time = " + str(transitionTimesFixed[num]))
    nx.draw(targetGraph, node_color=color_map, with_labels=True)
    plt.savefig(title)
                
graphs = graph_gen(targetGraph)
transitionMatrix = tm_generation(graphs)

"""
emptyZeroForcingSetVal = tuple("0")
start = time.process_time()
transitionTimesFixed = emptyZeroForcingSetVal + (propogation_time_solver(transitionMatrix)[0])
print("Time needed for old algorithm = " + str(time.process_time() - start) + " seconds")

start = time.process_time()
transitionTimesFixed = propogation_time_solver_2(transitionMatrix)
print("Time needed for calculating expected forcing time = " + str(time.process_time() - start) + " seconds")
print(transitionTimesFixed)
"""

start = time.process_time()
transitionTimesFixed = propogation_time_solver_3(transitionMatrix)
print("Time needed for NEW calculating expected forcing time = " + str(time.process_time() - start) + " seconds")
print(transitionTimesFixed)


"""
optimalSol = optimal_zfs_by_size(transitionTimesFixed, zeroForcingSetSize)
print_zfs((optimalSol[0])[0])
"""