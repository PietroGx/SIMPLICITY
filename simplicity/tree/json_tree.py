#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar  9 13:53:52 2025

@author: pietro
"""
import json
from anytree import Node

def get_node_attributes(node):
    """Extract all JSON-serializable attributes from an anytree Node, printing warnings for user-defined non-serializable ones."""
    serializable_attrs = {}
    for k, v in vars(node).items():
        if k.startswith("_NodeMixin"):  # Ignore internal anytree attributes
            continue
        if k in ('parent', 'children'):  # Explicitly exclude parent/children references
            continue
        try:
            json.dumps(v)  # Check if the value is JSON serializable
            serializable_attrs[k] = v
        except (TypeError, OverflowError):
            print(f"Warning: Attribute '{k}' of node '{node.name}' was not stored (non-serializable: {type(v).__name__})")
            # Attribute is skipped
    return serializable_attrs

def serialize(node):
    node_data = get_node_attributes(node)  # Get dynamic attributes
    node_data["name"] = node.name
    node_data["children"] = [serialize(child) for child in node.children]
    return node_data

def write_json_tree_file(root, json_tree_filepath):
    with open(json_tree_filepath, "w") as f:
        json.dump(serialize(root), f, indent=4)

def deserialize(data, parent=None):
    node_attributes = {k: v for k, v in data.items() if k not in ("name", "children")}
    node = Node(data["name"], parent=parent, **node_attributes)  # Recreate node
    for child in data.get("children", []):
        deserialize(child, node)  # Recursively add children
    return node

def import_tree(json_tree_filepath):
    """Import a tree from JSON and reconstruct the hierarchy with attributes."""
    with open(json_tree_filepath, "r") as f:
        data = json.load(f)
    return deserialize(data)


