import heapq

import sys


class PathGraph:
    def __init__(self):
        self.vertices = {}

    def add_vertex(self, name):
        self.vertices[name] = dict()

    def add_edge(self, from_vertex, to_vertex, distance):
        if from_vertex not in self.vertices:
            self.vertices[from_vertex] = {}
        if to_vertex not in self.vertices:
            self.vertices[to_vertex] = {}

        self.vertices[from_vertex][to_vertex] = distance
        self.vertices[to_vertex][from_vertex] = distance

    def shortest_path(self, start, finish):
        distances = {}  # Distance from start to node
        previous = {}  # Previous node in optimal path from source
        nodes = []  # Priority queue of all nodes in Graph
        for vertex in self.vertices:
            if vertex == start:  # Set root node as distance of 0
                distances[vertex] = 0
                heapq.heappush(nodes, [0, vertex])
            else:
                distances[vertex] = sys.maxsize
                heapq.heappush(nodes, [sys.maxsize, vertex])
            previous[vertex] = None

        while nodes:
            smallest = heapq.heappop(nodes)[1]  # Vertex in nodes with smallest distance in distances
            if smallest == finish:  # If the closest node is our target we're done so print the path
                path = []
                while previous[smallest]:  # Traverse through nodes til we reach the root which is 0
                    path.append(smallest)
                    smallest = previous[smallest]
                return path
            if distances[smallest] == sys.maxsize:  # All remaining vertices are inaccessible from source
                print('Breaking at', smallest)
                print(nodes)
                break

            for neighbor in self.vertices[smallest]:  # Look at all the nodes that this vertex is attached to
                #if self.vertices[smallest][neighbor] <= 1: print('ebesinin ami')
                alt = distances[smallest] + self.vertices[smallest][neighbor]  # Alternative path distance
                if alt < distances[neighbor]:  # If there is a new shortest path update our priority queue (relax)
                    distances[neighbor] = alt
                    previous[neighbor] = smallest
                    heapq.heappush(nodes, [alt, neighbor])
                    #for n in nodes:
                    #    if n[1] == neighbor:
                    #        n[0] = alt
                    #        break
                    #heapq.heapify(nodes)
        return distances

    def write_to_file(self, filename):
        f = open(filename, 'w')
        unique_vertices = set(self.vertices.keys())
        f.write(' ')
        for i in unique_vertices:
            for j in unique_vertices:
                if i == j:
                    f.write('MER_%s ' % i)
                elif j in self.vertices[i]:
                    f.write('%d ' % self.vertices[i][j])
                else:
                    f.write('# ')
            f.write('\n')

    def open_file(self, filename):
        f = open(filename, 'r')
        matrix = []
        for line in f:
            matrix.append(line.strip().split(' '))
            #print(matrix[-1])
        for i ,row in enumerate(matrix):
            for j, col in enumerate(row):
                if i!= j and col != '#':
                    self.add_edge(matrix[i][i], matrix[j][j], int(col))
                  


#graph = PathGraph()
#graph.open_file('deneme.out')
#print(graph.shortest_path('MER_8', 'MER_1398'))
#print(graph.vertices.keys())
