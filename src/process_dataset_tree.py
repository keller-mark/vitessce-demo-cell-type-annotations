import pandas as pd
import numpy as np
import json
from constants import *
from utils import *

if __name__ == "__main__":
    df = pd.read_csv(snakemake.input[0])
    df = df.loc[df[COLUMNS.DATASET_ID.value] == snakemake.wildcards[COLUMNS.DATASET_ID.value]]

    tree = init_tree()
    
    # Load the cell ontology DAG
    graph, id_to_name, name_to_id = load_co()

    
    for cell_type, cell_type_df in df.groupby(COLUMNS.ANNOTATION.value):

        try:
            node_id = name_to_id[cell_type]
        except KeyError:
            print(f"ERROR: annotation '{cell_type}' does not match any node in the cell ontology.")

        # Get ancestors of the cell type
        ancestor_term_set = networkx.descendants(graph, node_id) # counterintuitive

        # Make sure the cell type has an ancestor with the 'cell' root ID
        assert(CL_ROOT_ID in ancestor_term_set)
        
        # Initialize the current node ID to the cell type of interest
        curr_node_id = node_id
        
        # Construct a list of ancestors, inclusive. [node_id, parent_id, grandparent_id, ...]
        ancestors = [curr_node_id]

        # Get the parents of the current node.
        curr_node_parents = list(graph.out_edges(curr_node_id, keys=True))
        
        # Warn if the current node has multiple parents.
        if len(curr_node_parents) > 1:
            print(f"WARN: {curr_node_id} has {len(curr_node_parents)} parents.")

        # Select the first (hopefully only) parent node.
        _, curr_parent_id, relationship = curr_node_parents[0]
        assert(relationship == "is_a")

        while curr_node_id != CL_ROOT_ID:
            # Get the parents of the current node.
            curr_node_parents = list(graph.out_edges(curr_node_id, keys=True))
            
            # Warn if the current node has multiple parents.
            if len(curr_node_parents) > 1:
                print(f"WARN: {curr_node_id} has {len(curr_node_parents)} parents.")

            # Select the first (hopefully only) parent node.
            _, curr_parent_id, relationship = curr_node_parents[0]
            assert(relationship == "is_a")
            ancestors.append(curr_parent_id)

            # Set the current node to its parent to prepare for the next iteration.
            curr_node_id = curr_parent_id


        print(ancestors)
        print([ id_to_name[a] for a in ancestors ])

        """
        tree["tree"][0]["children"].append({
            "name": cell_type,
            "set": cell_type_df[COLUMNS.CELL_ID.value].unique().tolist()
        })
        """