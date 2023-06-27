import math
import random
import time
from collections import defaultdict
import numpy as np
import copy
import itertools
from numpy.linalg import norm
from scipy.special import comb
from itertools import combinations

class Hypergraph:
    def __init__(self, hyperedges, weightedEdges=False):
        self.addEdges(hyperedges, weightedEdges)
        self.deleteDegenerateHyperedges()
        self.findHyperedgeSizes()
        self.generateNeighbors()
        self.nodeLabels = list(self.nodes.keys())

    def __iter__(self):
        """Iterate over the nodes. Use: 'for n in G'.
        Returns
        -------
        niter : iterator
            An iterator over all nodes in the graph.
        Examples
        --------
        >>> G = nx.path_graph(4)  # or DiGraph, MultiGraph, MultiDiGraph, etc
        >>> [n for n in G]
        [0, 1, 2, 3]
        >>> list(G)
        [0, 1, 2, 3]
        """
        return iter(self.nodes)

    def __contains__(self, n):
        """Returns True if n is a node, False otherwise. Use: 'n in G'.
        Examples
        --------
        >>> G = nx.path_graph(4)  # or DiGraph, MultiGraph, MultiDiGraph, etc
        >>> 1 in G
        True
        """
        try:
            return n in self.nodes
        except TypeError:
            return False

    def __len__(self):
        """Returns the number of nodes in the graph. Use: 'len(G)'.
        Returns
        -------
        nnodes : int
            The number of nodes in the graph.
        See Also
        --------
        number_of_nodes, order  which are identical
        Examples
        --------
        >>> G = nx.path_graph(4)  # or DiGraph, MultiGraph, MultiDiGraph, etc
        >>> len(G)
        4
        """
        return len(self.nodes)

    def addEdges(self, hyperedges, weightedEdges):
        # unweighted format for hyperedges: {"id0":{"members":(1,2,3)}, "id1":{"members":(1,2)},...}
        # weighted format for hyperedges: {"id0":{"members":(1,2,3),"weight":1.1}, "id1":{"members":(1,2),"weight":0.5},...}
        self.weightedEdges = weightedEdges
        self.nodes = dict()
        nodes = set()
        # if list of tuples
        if isinstance(hyperedges, list):
            self.hyperedges = dict()
            uid = 0
            for hyperedge in hyperedges:
                if self.weightedEdges:
                    self.hyperedges[uid] = {"members":hyperedge[:-1],"weight":hyperedge[-1]}
                else:
                    self.hyperedges[uid] = {"members":hyperedge}
                    nodes.update(hyperedge)
                uid += 1

        elif isinstance(hyperedges, dict):
            self.hyperedges = hyperedges.copy()
            for edgeData in self.hyperedges.values():
                nodes.update(edgeData["members"])

        for nodeLabel in list(nodes):
            self.nodes[nodeLabel] = dict()
        # need a better way to check whether the format is correct

    def addNodeAttributes(self, nodeAttributes):
        # find unique nodes in the hyperedges
        for label, attribute in nodeAttributes.items():
            try:
                self.nodes[label] = attribute
            except:
                print("invalid label")

    def deleteDegenerateHyperedges(self):
        cleanedHyperedges = dict()
        for uid, hyperedge in self.hyperedges.items():
            if len(hyperedge["members"]) >= 2:
                cleanedHyperedges[uid] = hyperedge
        self.hyperedges = cleanedHyperedges

    def number_of_nodes(self):
        return len(self.nodes)

    def has_node(self, n):
        try:
            return n in self.nodes
        except TypeError:
            return False

    def findHyperedgeSizes(self):
        hyperedgeSizes = set()
        for edgeData in list(self.hyperedges.values()):
            hyperedgeSizes.add(len(edgeData["members"]))
        self.hyperedgeSizes = list(hyperedgeSizes)

    def getHyperedgeSizes(self):
        return self.hyperedgeSizes

    def generateNeighbors(self):
        self.neighbors = defaultdict(dict)
        if self.weightedEdges:
            self.generateWeightedNeighbors()
        else:
            self.generateUnweightedNeighbors()

    def generateUnweightedNeighbors(self):
        for uid, edgeData in self.hyperedges.items():
            try:
                members = edgeData["members"]
            except:
                print("Incorrect input format for hyperedge list")
                break
            for index in range(len(members)):
                self.neighbors[members[index]][uid] = {"neighbors":tuple([members[i] for i in range(len(members)) if i != index])}

    def generateWeightedNeighbors(self):
        for uid, edgeData in self.hyperedges.items():
            try:
                members = edgeData["members"]
                weight = edgeData["weight"]
            except:
                print("Incorrect input format for weighted hyperedge list")
                break
            for index in range(len(members)):
                self.neighbors[members[index]][uid] = {"neighbors":tuple([members[i] for i in range(len(members)) if i != index]), "weight":weight}

    def getDegreeSequenceBySize(self, hyperedgeSize):
        degreeSequence = defaultdict(lambda: 0)
        for uid, hyperedge in self.hyperedges.items():
            if len(hyperedge["members"]) == hyperedgeSize:
                for node in hyperedge["members"]:
                    degreeSequence[node] += 1
        return degreeSequence

    def getHyperedgesBySize(self, hyperedgeSize):
        hyperedges = list()
        for uid, data in self.hyperedges.items():
            if len(data["members"]) == hyperedgeSize:
                hyperedges.append(data["members"])
        return hyperedges


class HypergraphGenerator:
    def __init__(self, data, type="structure"):
        if type == "structure":
            self.generateHyperdegreeSequence(data)
            self.generateHyperedges()
        elif type == "hyperedge-list":
            self.hyperedgeListToDictionary(data)
            self.getHyperdegreeSequence()
        elif type == "hyperedge-dictionary":
            self.hyperedges = copy.deepcopy(data)
            self.getHyperdegreeSequence()
        else:
            print("invalid option")

    def getHyperedges(self):
        return self.hyperedges

    def getHyperedgeList(self):
        return [e["members"] for e in list(self.hyperedges.values())]

    def hyperedgeListToDictionary(self, hyperedges):
        self.hyperedges = dict()
        uid = 0
        for hyperedge in hyperedges:
            self.hyperedges[uid] = {"members":hyperedge}
            uid += 1

    def getHyperdegreeSequence(self):
        try:
            return self.hyperdegreeSequence
        except:
            self.hyperdegreeSequence = defaultdict(lambda: defaultdict(lambda: 0))
            for uid, data in self.hyperedges.items():
                members = data["members"]
                hyperedgeSize = len(members)
                for node in members:
                    self.hyperdegreeSequence[node][hyperedgeSize] += 1

    def generateHyperdegreeSequence(self, structure):
        self.hyperdegreeSequence = defaultdict(lambda: defaultdict(lambda: 0))
        self.hyperedgeSizes = list()
        correlatedSequence = list()
        for info in structure:
            try:
                hyperedgeSize = info["hyperedge-size"]
                self.hyperedgeSizes.append(hyperedgeSize)
            except:
                print("Error in specified distribution parameters")

            if info["degree-distribution"] == "power-law":
                try:
                    numNodes = info["size"]
                    minDegree = info["min-degree"]
                    maxDegree = info["max-degree"]
                    exponent = info["exponent"]
                except:
                    print("Error in specified distribution parameters")
                    break
                sequence = self.generatePowerLawDegreeSequence(numNodes, minDegree, maxDegree, exponent)
            elif info["degree-distribution"] == "uniform":
                try:
                    numNodes = info["size"]
                    minDegree = info["min-degree"]
                    maxDegree = info["max-degree"]
                except:
                    print("Error in specified distribution parameters")
                    break
                sequence = self.generateUniformDegreeSequence(numNodes, minDegree, maxDegree)
            elif info["degree-distribution"] == "poisson":
                try:
                    numNodes = info["size"]
                    meanDegree = info["mean-degree"]
                except:
                    print("Error in specified distribution parameters")
                    break
                sequence = self.generatePoissonDegreeSequence(numNodes, meanDegree)
            else:
                print("Invalid selection")
                break
            try:
                isCorrelated = info["is-correlated"]
            except:
                print("Specify whether this hyperedge size is correlated or not")
                break
            if isCorrelated:
                if correlatedSequence == []:
                    correlatedSequence = sequence
                self.updateHyperdegreeSequence(correlatedSequence, hyperedgeSize)
                #self.hyperdegreeSequence[hyperedgeSize] = correlatedSequence
            else:
                self.updateHyperdegreeSequence(sequence, hyperedgeSize)
                #self.hyperdegreeSequence[hyperedgeSize] = sequence

    def updateHyperdegreeSequence(self, sequence, hyperedgeSize):
        try:
            for node in range(len(sequence)):
                self.hyperdegreeSequence[node][hyperedgeSize] = sequence[node]
        except:
            self.hyperdegreeSequence = defaultdict(lambda: defaultdict(lambda: 0))
            for node in range(len(sequence)):
                self.hyperdegreeSequence[node][hyperedgeSize] = sequence[node]

    def generateHyperedges(self):
        self.hyperedges = dict()
        for hyperedgeSize in self.hyperedgeSizes:
            self.hyperedges.update(self.generateHyperedgesBySize(hyperedgeSize))

    def generateHyperedgesBySize(self, hyperedgeSize):
        import string
        k = dict()
        for node, hyperdegree in self.hyperdegreeSequence.items():
            k[node] = hyperdegree[hyperedgeSize]

        # Making sure we have the right number of stubs
        if (sum(k.values()) % hyperedgeSize) != 0:
            remainder = sum(k.values()) % hyperedgeSize
            for i in range(int(round(hyperedgeSize - remainder))):
                j = random.randrange(len(k))
                k[j] = k[j] + 1

        stubs = list()
        hyperedges = dict()
        # Creating the list to index through
        for index in range(len(k)):
            stubs.extend([index]*int(k[index]))

        while len(stubs) != 0:
            uid = ''.join(random.choice(string.ascii_lowercase) for i in range(8))
            u = random.sample(range(len(stubs)), hyperedgeSize)
            hyperedge = list()
            for index in u:
                hyperedge.append(stubs[index])

            hyperedges[uid] = {"members":tuple(hyperedge)}

            for index in sorted(u, reverse=True):
                del stubs[index]
        return hyperedges

    def generatePowerLawDegreeSequence(self, numNodes, minDegree, maxDegree, exponent):
        degreeSequence = list()
        for i in range(numNodes):
            u = random.uniform(0, 1)
            degreeSequence.append(round(self.invCDFPowerLaw(u, minDegree, maxDegree, exponent)))
        return degreeSequence # originally this was sorted but I'm worried about between-size correlations

    def invCDFPowerLaw(self, u, minDegree, maxDegree, exponent):
        return (minDegree**(1-exponent) + u*(maxDegree**(1-exponent) - minDegree**(1-exponent)))**(1/(1-exponent))

    def generateUniformDegreeSequence(self, numNodes, minDegree, maxDegree):
        degreeSequence = list()
        for i in range(numNodes):
            u = random.randrange(round(minDegree), round(maxDegree))
            degreeSequence.append(round(u))
        return degreeSequence

    def generatePoissonDegreeSequence(self, numNodes, meanDegree):
        return np.random.poisson(lam=meanDegree, size=numNodes).tolist()

    def getHyperedgesBySize(self, hyperedgeSize, ids=True):
        if ids:
            hyperedgesBySize = dict()
            for uid, hyperedgeData in self.hyperedges.items():
                if len(hyperedgeData["members"]) == hyperedgeSize:
                    hyperedgesBySize[uid] = dict(hyperedgeData)
            return hyperedgesBySize
        else:
            hyperedgesBySize = list()
            for hyperedgeData in list(self.hyperedges.values()):
                if len(hyperedgeData["members"]) == hyperedgeSize:
                    hyperedgesBySize.append(hyperedgeData["members"])
            return hyperedgesBySize

    def updateHyperedgesAfterShuffle(self, hyperedgeListBySize):
        for key, hyperedge in hyperedgeListBySize.items():
            self.hyperedges[key] = dict(hyperedge)

    def updateMeanDegreeProduct(self, meanDegreeProduct, degrees1, degrees2, newDegrees1, newDegrees2, m, numEdges):
        numCombos = comb(m, 2)
        degreeProduct1 = sum([k[0]*k[1]/numCombos for k in combinations(degrees1, 2)])
        degreeProduct2 = sum([k[0]*k[1]/numCombos for k in combinations(degrees2, 2)])
        newDegreeProduct1 = sum([k[0]*k[1]/numCombos for k in combinations(newDegrees1, 2)])
        newDegreeProduct2 = sum([k[0]*k[1]/numCombos for k in combinations(newDegrees2, 2)])
        return meanDegreeProduct + (newDegreeProduct1 + newDegreeProduct2 - degreeProduct1 - degreeProduct2)/numEdges

    def getMeanDegreeProduct(self, edgeList, k, m):
        numEdges = len(edgeList)
        numCombos = comb(m, 2)
        meanDegreeProduct = sum([np.sum(np.prod(k[list(indices)]) for indices in combinations(edge, 2))/(numEdges*numCombos) for edge in edgeList])
        return meanDegreeProduct

    def getMeanPowerOfDegree(self, k, power=1):
        return np.mean(np.power(k, power))

    def getHyperdegreeSequenceBySize(self, m):
        d =  {node: hyperdegree[m] for node, hyperdegree in self.hyperdegreeSequence.items()}
        degreeSequence = np.zeros(max(d.keys()) + 1)
        for node, degree in d.items():
            degreeSequence[node] = degree
        return degreeSequence

    def shuffleHyperedges(self, m, maxIterations=10000):

        shuffledEdges = self.getHyperedgesBySize(m, ids=True)

        iteration = 0

        HEdges = [list(np.sort(shuffledEdges[idx]["members"])) for idx in shuffledEdges.keys()]

        while (iteration < maxIterations):
            uid1, uid2 = random.sample(shuffledEdges.keys(), 2)
            member1, member2 = random.choices(range(m), k=2)
            edge1 = shuffledEdges[uid1]["members"]
            edge2 = shuffledEdges[uid2]["members"]
            newEdge1 = list(edge1)
            newEdge2 = list(edge2)
            newEdge1[member1] = edge2[member2]
            newEdge2[member2] = edge1[member1]

            newEdge1.sort()
            newEdge2.sort()
            #HEdges=[list(np.sort(shuffledEdges[idx]["members"])) for idx in shuffledEdges.keys()]

            # ensure that we don't create fully coincident hyperedges, that size(e)=size(e')=size(f)=size(f'), i.e. we don't make a swap with a node already present in the hyperedge
            if len(set(edge1)) == len(set(newEdge1)) and len(set(edge2)) == len(set(newEdge2)) and (newEdge1 not in HEdges) and (newEdge2 not in HEdges):
                shuffledEdges[uid1]["members"] = tuple(newEdge1)
                shuffledEdges[uid2]["members"] = tuple(newEdge2)
                iteration += 1

                HEdges[list(shuffledEdges.keys()).index(uid1)]=newEdge1
                HEdges[list(shuffledEdges.keys()).index(uid2)]=newEdge2

        #print("Number of double-edge swaps: " + str(iteration), flush=True)
        self.updateHyperedgesAfterShuffle(shuffledEdges)
        return;

    def boltzmannProbability(self, difference, temperature):
        if difference >= 0:
            return 1
        if temperature == 0:
            return 0
        else:
            return math.exp(difference/temperature)

    def getAssortativity(self, m):
        edgeList = self.getHyperedgesBySize(m, ids=False)
        k = self.getHyperdegreeSequenceBySize(m)
        k1 = self.getMeanPowerOfDegree(k, power=1)
        k2 = self.getMeanPowerOfDegree(k, power=2)
        kk1 = self.getMeanDegreeProduct(edgeList, k, m)        

        return (m - 1)*(k1/k2)**2*kk1/2 - 1