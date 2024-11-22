import pandas as pd
import pickle

# Load glycan binding data
glycan_binding = pd.read_csv('data/glycan_binding.csv')
binding_10 = glycan_binding.head(10).dropna(axis=1)

with open('data/glycan_graphs.pkl', 'rb') as file:
    data = pickle.load(file)  # Must be 'rb' mode

structure_name = 'GlcNAc(b1-6)GalNAc'
graph = data[structure_name]

# Viewing nodes with attributes
nodes_with_attributes = list(graph.nodes(data=True))
# Create a DataFrame with nodes and attributes
nodes_with_attributes_df = pd.DataFrame.from_records(
    [(node, attrs) for node, attrs in graph.nodes(data=True)],
    columns=['Node', 'Attributes']
)

#make nodes_with_attributes_df columns attributes into separate columns
nodes_with_attributes_df = pd.concat([nodes_with_attributes_df.drop(['Attributes'], axis=1),
                                      nodes_with_attributes_df['Attributes'].apply(pd.Series)], axis=1)


nodes_with_attributes_df.to_csv('/Users/xakdze/PycharmProjects/GlycoShape/out/nodes_with_attributes.csv',sep="\t", index=False)

binding_10.to_csv('/Users/xakdze/PycharmProjects/GlycoShape/out/binding_10.csv', sep="\t",index=False)