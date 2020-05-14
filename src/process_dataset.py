import pandas as pd
import numpy as np
import json
import networkx
import obonet

from constants import *

if __name__ == "__main__":
    df = pd.read_csv(snakemake.input[0])
    df = df.loc[df[COLUMNS.DATASET_ID.value] == snakemake.wildcards[COLUMNS.DATASET_ID.value]]
    
    graph = obonet.read_obo(CL_OBO_URL)
    
    # make sure there are no cycles
    assert(networkx.is_directed_acyclic_graph(graph))

    id_to_name = {id_: data.get('name') for id_, data in graph.nodes(data=True)}
    name_to_id = {data['name']: id_ for id_, data in graph.nodes(data=True) if 'name' in data}


    # Get ancestors of B cell
    ancestor_term_set = networkx.descendants(graph, 'CL:0000236') # counterintuitive

    # make sure this is a cell
    assert('CL:0000000' in ancestor_term_set)

    
    curr_node_id = 'CL:0000236' # B cell
    ancestors = [curr_node_id]
    direct_ancestors = list(graph.out_edges(curr_node_id, keys=True))
    if len(direct_ancestors) > 1:
        print(f"WARNING: {curr_node_id} has {len(direct_ancestors)} parents.")
    _, parent_id, relationship = direct_ancestors[0]

    while curr_node_id != 'CL:0000000':
        direct_ancestors = list(graph.out_edges(curr_node_id, keys=True))
        if len(direct_ancestors) > 1:
            print(f"WARNING: {curr_node_id} has {len(direct_ancestors)} parents.")
        _, parent_id, relationship = direct_ancestors[0]
        ancestors.append(parent_id)
        curr_node_id = parent_id


    print(ancestors)
    print([ id_to_name[a] for a in ancestors ])