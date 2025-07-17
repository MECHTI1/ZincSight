#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 26 11:54:16 2023

@author: mechti
"""

# Python program for Kruskal's algorithm to find
# Minimum Spanning Tree of a given connected,
# undirected and weighted graph

# Class to represent a graph


class Graph:

	def __init__(self, vertices):
		self.V = vertices # No. of vertices
		self.graph = []
		# to store graph

	# function to add an edge to graph
	def addEdge(self, u, v, w):
		self.graph.append([u, v, w])

	# A utility function to find set of an element i
	# (truly uses path compression technique)
	def find(self, parent, i):
		if parent[i] != i:
		# Reassignment of node's parent to root node as
		# path compression requires
			parent[i] = self.find(parent, parent[i])
		return parent[i]

	# A function that does union of two sets of x and y
	# (uses union by rank)
	def union(self, parent, rank, x, y):

		# Attach smaller rank tree under root of
		# high rank tree (Union by Rank)
		if rank[x] < rank[y]:
			parent[x] = y
		elif rank[x] > rank[y]:
			parent[y] = x

		# If ranks are same, then make one as root
		# and increment its rank by one
		else:
			parent[y] = x
			rank[x] += 1

	# The main function to construct MST using Kruskal's
		# algorithm
	def KruskalMST(self):

		result = [] # This will store the resultant MST

		# An index variable, used for sorted edges
		i = 0

		# An index variable, used for result[]
		e = 0

		# Step 1: Sort all the edges in
		# non-decreasing order of their
		# weight. If we are not allowed to change the
		# given graph, we can create a copy of graph
		self.graph = sorted(self.graph,
							key=lambda item: item[2])

		parent = []
		rank = []

		# Create V subsets with single elements
		for node in range(self.V):
			parent.append(node)
			rank.append(0)

		# Number of edges to be taken is less than to V-1
		while e < self.V - 1:

			# Step 2: Pick the smallest edge and increment
			# the index for next iteration
			u, v, w = self.graph[i]
			i = i + 1
			x = self.find(parent, u)
			y = self.find(parent, v)

			# If including this edge doesn't
			# cause cycle, then include it in result
			# and increment the index of result
			# for next edge
			if x != y:
				e = e + 1
				result.append([u, v, w])
				self.union(parent, rank, x, y)
			# Else discard the edge

		minimumCost = 0
		print("Edges in the constructed MST")
		for u, v, weight in result:
			minimumCost += weight
			print("%d -- %d == %d" % (u, v, weight))
		print("Minimum Spanning Tree", minimumCost)

		return(result)


# Driver's code
def main(list_nodes, list_edges):

	print("running kruskal algorithm...")
	print ("list_nodes",list_nodes)
	print ("list_edges",list_edges)
	g = Graph(len(list_nodes))
	for edge in list_edges:
		u, v, w = edge
		g.addEdge(list_nodes.index(u), list_nodes.index(v), w)

	# Function call
	Kruskal_result=g.KruskalMST()

	lst_edges_by_resi_names_not_num = [[list_nodes[sublist[0]], list_nodes[sublist[1]], sublist[2]] for sublist in Kruskal_result]

	print(lst_edges_by_resi_names_not_num)
	return (lst_edges_by_resi_names_not_num)


if __name__ == '__main__':
	list_nodes= ['A_273', 'A_277', 'A_465']
	list_edges =[['A_273', 'A_277', 6], ['A_273', 'A_465', 9], ['A_277', 'A_465', 8]]
	main(list_nodes, list_edges)


#Gilad- type_chain_resi1, type_chain_resi2, alpha_distance


# This code is contributed by Neelam Yadav
# Improved by James Graça-Jones


