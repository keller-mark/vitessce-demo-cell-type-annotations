#!/usr/bin/env python3

import networkx
from pprint import pprint

from constants import COLUMNS, CL_ROOT_ID
from utils import init_tree, load_cl_obo_graph


def generate_flat_cell_sets(df):
    tree = init_tree()

    leiden_clusters_children = []
    for cluster_name, cluster_df in df.groupby("leiden"):
        leiden_clusters_children.append({
            "name": cluster_name,
            "set": cluster_df[COLUMNS.CELL_ID.value].unique().tolist(),
            "itemtype": "static"
        })

    tree["tree"].append({
        "name": "Leiden Clustering",
        "children": leiden_clusters_children
    })

    cell_type_annotation_children = []
    for cell_type, cell_type_df in df.groupby(COLUMNS.ANNOTATION.value):
        set_cell_ids = cell_type_df[COLUMNS.CELL_ID.value].values.tolist()
        set_cell_scores = cell_type_df[COLUMNS.PREDICTION_SCORE.value].values.tolist()
        set_value = [ list(x) for x in zip(set_cell_ids, set_cell_scores) ]

        cell_type_annotation_children.append({
            "name": cell_type,
            "set": set_value,
            "itemtype": "probabilistic"
        })

    tree["tree"].append({
        "name": "Cell Type Annotations (flat)",
        "children": cell_type_annotation_children
    })
    return tree


def generate_hierarchical_cell_sets(df, cl_obo_file):
    tree = init_tree()

    # Load the cell ontology DAG
    graph, id_to_name, name_to_id = load_cl_obo_graph(cl_obo_file)


    def check_multi_parents(node_id, node_parent_ids):
        # Warn if the current node has multiple parents.
        num_parents = len(node_parent_ids)
        if num_parents > 1:
            parent_names = "[ " + ", ".join([ f"{p_id} ({id_to_name[p_id]})" for p_id, _, _ in node_parent_ids ]) + " ]"
            """
            print((
                f"WARN: {node_id} ({id_to_name[node_id]}) has "
                f"{num_parents} parents: {parent_names}."
            ))
            """
    
    def get_parents(node_id):

        if node_id == CL_ROOT_ID:
            return [
                [node_id]
            ]

        # Get ancestors of the cell type
        # (counterintuitive that the function is called descendants)
        ancestor_term_set = networkx.descendants(graph, node_id)

        # Make sure the cell type has an ancestor
        # with the 'cell' root ID
        assert(CL_ROOT_ID in ancestor_term_set)

        # Get the parents of the current node.
        node_parents = list(graph.out_edges(node_id, keys=True))

        up_dag_paths = []
        for node_parent in node_parents:
            _, curr_parent_id, relationship = node_parent
            if relationship == "is_a":
                parent_paths = get_parents(curr_parent_id)
                for parent_path in parent_paths:
                    up_dag_paths.append([node_id] + parent_path)
        return up_dag_paths

    ancestors_and_sets = []

    for cell_type, cell_type_df in df.groupby(COLUMNS.ANNOTATION.value):

        try:
            node_id = name_to_id[cell_type]
        except KeyError:
            print((
                f"ERROR: annotation '{cell_type}' does "
                "not match any node in the cell ontology."
            ))

        # Get ancestors of the cell type
        # (counterintuitive that the function is called descendants)
        ancestor_term_set = networkx.descendants(graph, node_id)

        # Make sure the cell type has an ancestor
        # with the 'cell' root ID
        assert(CL_ROOT_ID in ancestor_term_set)

        # Initialize the current node ID to the cell type of interest
        curr_node_id = node_id

        paths_up = get_parents(node_id)
        named_paths_up = [ [id_to_name[n_id] for n_id in path_nodes] for path_nodes in paths_up ]
        print()
        print(f"{id_to_name[node_id]} has {len(paths_up)} paths up to {CL_ROOT_ID} ({id_to_name[CL_ROOT_ID]}):")
        for named_path_nodes in named_paths_up:
            print(named_path_nodes)
        print()

        # Construct a list of ancestors, inclusive.
        # [node_id, parent_id, grandparent_id, ...]
        ancestors = [curr_node_id]

        # Get the parents of the current node.
        curr_node_parents = list(graph.out_edges(curr_node_id, keys=True))

        check_multi_parents(curr_node_id, curr_node_parents)

        # Select the first (hopefully only) parent node.
        _, curr_parent_id, relationship = curr_node_parents[0]
        assert(relationship == "is_a")

        while curr_node_id != CL_ROOT_ID:
            # Get the parents of the current node.
            curr_node_parents = list(graph.out_edges(curr_node_id, keys=True))

            # Warn if the current node has multiple parents.
            check_multi_parents(curr_node_id, curr_node_parents)

            # Select the first (hopefully only) parent node.
            _, curr_parent_id, relationship = curr_node_parents[0]
            assert(relationship == "is_a")
            ancestors.append(curr_parent_id)

            # Set the current node to its parent to
            # prepare for the next iteration.
            curr_node_id = curr_parent_id

        named_ancestors = [id_to_name[a] for a in ancestors]
        named_ancestors_reversed = list(reversed(named_ancestors))
        #print(named_ancestors_reversed)

        set_cell_ids = cell_type_df[COLUMNS.CELL_ID.value].values.tolist()
        set_cell_scores = cell_type_df[COLUMNS.PREDICTION_SCORE.value].values.tolist()
        set_value = [ list(x) for x in zip(set_cell_ids, set_cell_scores) ]

        ancestors_and_sets.append((
            named_ancestors_reversed,
            set_value
        ))

    # Pop off all ancestors that are the same for all cell types.
    # e.g. 'cell', 'native cell', ...
    ancestor_list_lens = [len(x[0]) for x in ancestors_and_sets]
    min_ancestor_list_len = min(ancestor_list_lens)
    assert(min_ancestor_list_len >= 1)
    for level in range(min_ancestor_list_len - 1):
        unique_level_cell_types = set()
        for ancestors, cell_set in ancestors_and_sets:
            unique_level_cell_types.add(ancestors[0])

        if len(unique_level_cell_types) == 1:
            for ancestors, cell_set in ancestors_and_sets:
                ancestors.pop(0)
        else:
            #print(unique_level_cell_types)
            break

    # Construct a hierarchy of cell types.
    def find_or_create_parent(d, keys, child):
        key = keys[0]

        if key in d and isinstance(d[key], dict):
            result = d[key]
        else:
            result = d[key] = dict()

        if len(keys) == 1:
            result["any"] = child
            return result
        else:
            new_keys = keys.copy()
            new_keys.pop(0)
            return find_or_create_parent(result, new_keys, child)

    h = dict()
    for ancestors, cell_set in ancestors_and_sets:
        find_or_create_parent(h, ancestors, cell_set)

    def to_tree(name, value):
        if isinstance(value, dict):
            return {
                "name": name,
                "children": [
                    to_tree(child_name, child_value)
                    for child_name, child_value in value.items()
                ]
            }
        else:
            return {
                "name": name,
                "set": value,
                "itemtype": "probabilistic"
            }

    tree["tree"] = [
        to_tree("Cell Type Annotations (hierarchical)", h)
    ]
    return tree
