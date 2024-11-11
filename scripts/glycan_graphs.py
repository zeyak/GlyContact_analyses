import pickle

with open('//data/glycan_graphs.pkl', 'rb') as file:
    data = pickle.load(file)  # Must be 'rb' mode

# Access a specific graph by its key (structure name)
structure_name = 'GlcNAc(b1-6)GalNAc'
#structure_name = 'Neu5Ac(a2-3)Gal(b1-4)GlcNAc(b1-2)[Neu5Ac(a2-3)Gal(b1-4)GlcNAc(b1-4)]Man(a1-3)[Neu5Ac(a2-3)Gal(b1-4)GlcNAc(b1-2)[Neu5Ac(a2-3)Gal(b1-4)GlcNAc(b1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc'
graph = data[structure_name]

# To view all nodes in this graph
nodes = list(graph.nodes)
print("Nodes:", nodes)

# To view all edges in this graph
edges = list(graph.edges)
print("Edges:", edges)

# Viewing nodes with attributes
nodes_with_attributes = list(graph.nodes(data=True))
print("Nodes with attributes:", nodes_with_attributes)


import networkx as nx
import matplotlib.pyplot as plt

# Draw the specific graph with labels
plt.figure(figsize=(8, 6))
nx.draw(graph, with_labels=True, node_color='lightblue', font_size=10)
plt.show()


