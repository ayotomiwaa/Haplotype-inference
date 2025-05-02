import pandas as pd
import numpy as np
import networkx as nx
from itertools import product
from scipy.stats import norm
import os
import heapq
from functools import lru_cache
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.lines import Line2D

def generate_possible_cn_pairs(max_copy_number=5):
    """Generate possible copy number pairs with a+b <= max_copy_number"""
    return [(a, b) for a, b in product(range(max_copy_number + 1), repeat=2) 
            if 0 < a + b <= max_copy_number]

def calculate_expected_baf(a, b):
    """Calculate expected BAF given copy numbers, with zero division protection"""
    if a + b == 0:
        return 0.5
    return min(a, b) / (a + b)

def baf_likelihood(observed_baf, expected_baf, sigma):
    """Calculate likelihood of observed BAF given expected BAF using Gaussian model"""
    return -np.log(norm.pdf(observed_baf, expected_baf, sigma) + 1e-10)

def calculate_scale_factor(total_cn, rdr):
    """Calculate scale factor with zero division protection"""
    if rdr <= 0:
        return float('inf')
    return total_cn / rdr

def compute_edge_weight(a1, b1, a2, b2, rdr1, rdr2):
    """Compute edge weight based on scale factor consistency"""
    try:
        # Cache total copy numbers
        total_cn1 = a1 + b1
        total_cn2 = a2 + b2
        
        # Safely calculate scale factors
        if rdr1 <= 0 or rdr2 <= 0:
            return 10.0  # Default high weight for unreliable data
            
        scale1 = total_cn1 / rdr1
        scale2 = total_cn2 / rdr2
        
        return abs(scale1 - scale2)
    except (ZeroDivisionError, TypeError, ValueError):
        return 10.0  # Default high weight for any calculation errors

def build_multipartite_graph(rdr_values, baf_values, possible_cn_pairs, baf_sigma):
    """
    Build a multipartite graph for copy number estimation with enhanced efficiency
    
    Parameters:
    - rdr_values: RDR values for each bin in a cell
    - baf_values: BAF values for each bin in a cell
    - possible_cn_pairs: List of possible (a,b) copy number pairs
    
    Returns:
    - G: NetworkX DiGraph representing the multipartite graph
    """
    num_bins = len(rdr_values)
    G = nx.DiGraph()
    
    # Precompute total copy numbers for each pair
    total_cn_map = {(a, b): a + b for a, b in possible_cn_pairs}
    
    # Add nodes for each bin and possible copy number pair
    for bin_idx in range(num_bins):
        for a, b in possible_cn_pairs:
            node_id = (bin_idx, a, b)
            
            # Get observed BAF for this bin
            observed_baf = baf_values[bin_idx]
            
            # Calculate expected BAF based on copy numbers
            expected_baf = calculate_expected_baf(a, b)
            
            # Use Gaussian model for BAF likelihood
            node_weight = baf_likelihood(observed_baf, expected_baf, sigma=baf_sigma)
            
            G.add_node(node_id, weight=node_weight)
    
    # Add edges between adjacent bins
    for bin_idx in range(num_bins - 1):
        for a1, b1 in possible_cn_pairs:
            node1 = (bin_idx, a1, b1)
            
            for a2, b2 in possible_cn_pairs:
                node2 = (bin_idx + 1, a2, b2)
                
                # Edge weight based on scale factor consistency
                edge_weight = compute_edge_weight(
                    a1, b1, a2, b2, 
                    rdr_values[bin_idx], rdr_values[bin_idx + 1]
                )
                
                G.add_edge(node1, node2, weight=edge_weight)
    
    return G

def find_shortest_path(G, num_bins, possible_cn_pairs):
    """
    Find shortest path through the multipartite graph using optimized Dijkstra
    
    Parameters:
    - G: NetworkX DiGraph
    - num_bins: Number of bins
    - possible_cn_pairs: List of possible (a,b) copy number pairs
    
    Returns:
    - path: List of (bin_idx, a, b) tuples representing the shortest path
    """
    # Create a super source node
    super_source = 'source'
    G.add_node(super_source, weight=0)
    
    # Connect super source to all nodes in first bin
    for a, b in possible_cn_pairs:
        first_node = (0, a, b)
        if G.has_node(first_node):
            G.add_edge(super_source, first_node, weight=0)
    
    # Create a super sink node
    super_sink = 'sink'
    G.add_node(super_sink, weight=0)
    
    # Connect all nodes in last bin to super sink
    last_bin = num_bins - 1
    for a, b in possible_cn_pairs:
        last_node = (last_bin, a, b)
        if G.has_node(last_node):
            G.add_edge(last_node, super_sink, weight=0)
    
    try:
        # Combine edge and node weights for Dijkstra
        for u, v in G.edges():
            if u != super_source and v != super_sink:
                G[u][v]['weight'] += G.nodes[v]['weight']
        
        # Find shortest path using single Dijkstra call
        path = nx.shortest_path(G, super_source, super_sink, weight='weight')
        
        # Remove super source and sink
        path = path[1:-1]
        
    except (nx.NetworkXNoPath, nx.NodeNotFound) as e:
        print(f"Error finding path: {e}")
        path = None
    
    # Clean up: remove super source and sink
    G.remove_node(super_source)
    G.remove_node(super_sink)
    
    return path

def load_and_prepare_data(rdr_file, baf_file):
    """
    Load and prepare RDR and BAF data with error handling
    
    Parameters:
    - rdr_file: Path to RDR TSV file
    - baf_file: Path to BAF TSV file with pre-computed BAF values
    
    Returns:
    - cell_data: Dictionary with cell data organized by cell ID
    """
    try:
        # Load RDR data
        rdr_df = pd.read_csv(rdr_file, sep='\t')
        
        # Handle different column formats flexibly
        if len(rdr_df.columns) == 7:
            if 'CELL' not in rdr_df.columns:
                rdr_df.columns = ['CHROMOSOME', 'START', 'END', 'CELL','NORMAL', 'COUNT', 'RDR']
        else:
            raise ValueError(f"Unexpected RDR file format with {len(rdr_df.columns)} columns")
        # print(f"Loaded RDR data with {rdr_df.head()} rows")
        # Load BAF data - format: CHROMOSOME START END CELL BAF
        try:
            # Try with header first
            baf_df = pd.read_csv(baf_file, sep='\t', header=None)
            if 'BAF' not in baf_df.columns and len(baf_df.columns) == 5:
                # No header, need to set column names
                baf_df.columns = ['CHROMOSOME', 'START', 'END', 'CELL', 'BAF']
        except:
            # If that fails, try without header
            baf_df = pd.read_csv(baf_file, sep='\t')
            if len(baf_df.columns) == 5:
                baf_df.columns = ['CHROMOSOME', 'START', 'END', 'CELL', 'BAF']
            else:
                raise ValueError(f"Unexpected BAF file format with {len(baf_df.columns)} columns")
        
        # Create a unique bin identifier for each chromosome/start/end combination
        rdr_df['BIN_ID'] = rdr_df['CHROMOSOME'] + '_' + rdr_df['START'].astype(str) + '_' + rdr_df['END'].astype(str)
        baf_df['BIN_ID'] = baf_df['CHROMOSOME'] + '_' + baf_df['START'].astype(str) + '_' + baf_df['END'].astype(str)
        
        # Get unique cells - prioritize BAF cells, as they're the ones we need
        cells = baf_df['CELL'].unique()
        print(cells)
        
        # Create a dictionary to store data by cell
        cell_data = {}
        
        for cell in cells:
            try:
                # Filter RDR data for this cell
                if cell in rdr_df['CELL'].values:
                    cell_rdr = rdr_df[rdr_df['CELL'] == cell].sort_values(['CHROMOSOME', 'START'])
                    
                    # Filter BAF data for this cell
                    cell_baf = baf_df[baf_df['CELL'] == cell].sort_values(['CHROMOSOME', 'START'])
            
                    
                    # Ensure we have both RDR and BAF for each bin
                    common_bins = set(cell_rdr['BIN_ID']).intersection(set(cell_baf['BIN_ID']))
                    
                    # Skip cells with no common bins
                    if not common_bins:
                        print(f"No common bins for cell {cell}, skipping.")
                        continue
                    
                    # Filter to keep only bins with both RDR and BAF data
                    cell_rdr = cell_rdr[cell_rdr['BIN_ID'].isin(common_bins)]
                    cell_baf = cell_baf[cell_baf['BIN_ID'].isin(common_bins)]
                    
                    # Sort both by bin ID to ensure alignment
                    cell_rdr = cell_rdr.sort_values('BIN_ID')
                    cell_baf = cell_baf.sort_values('BIN_ID')
                    
                    # Store cell data
                    cell_data[cell] = {
                        'rdr_values': cell_rdr['RDR'].values,
                        'baf_values': cell_baf['BAF'].values,
                        'chromosomes': cell_rdr['CHROMOSOME'].values,
                        'starts': cell_rdr['START'].values,
                        'ends': cell_rdr['END'].values,
                        'bin_ids': cell_rdr['BIN_ID'].values
                    }
            except Exception as e:
                print(f"Error processing cell {cell}: {e}")
                continue
        
        return cell_data
        
    except Exception as e:
        print(f"Error loading data: {e}")
        return {}

def process_cell(cell, data, possible_cn_pairs, baf_sigma):
    """Process a single cell to find optimal CN path"""
    try:
        # Build graph for this cell
        G = build_multipartite_graph(data['rdr_values'], data['baf_values'], possible_cn_pairs, baf_sigma)
        
        # Find shortest path
        path = find_shortest_path(G, len(data['rdr_values']), possible_cn_pairs)
        
        return path
    except Exception as e:
        print(f"Error processing cell {cell}: {e}")
        return None

def infer_haplotype_specific_copy_numbers(rdr_file, baf_file, output_file, max_copy_number=5, baf_sigma=0.05):
    """
    Main function to infer haplotype-specific copy numbers
    
    Parameters:
    - rdr_file: Path to RDR TSV file
    - baf_file: Path to BAF TSV file with pre-computed BAF values
    - output_file: Path to save output TSV file
    - max_copy_number: Maximum total copy number to consider
    """
    try:
        # Generate possible copy number pairs - do this once
        possible_cn_pairs = generate_possible_cn_pairs(max_copy_number)
        
        # Load and prepare data
        cell_data = load_and_prepare_data(rdr_file, baf_file)
        
        if not cell_data:
            print("No valid cell data found. Exiting.")
            return None
        
        # Process each cell
        results = []
        
        for cell, data in cell_data.items():
            print(f"Processing cell {cell}...")
            
            # Skip cells with no data
            if len(data['rdr_values']) == 0:
                print(f"  No data for cell {cell}, skipping.")
                continue
            
            # Process cell
            path = process_cell(cell, data, possible_cn_pairs, baf_sigma)
            
            if path is None:
                print(f"  No valid path found for cell {cell}, skipping.")
                continue
            
            # Extract results
            for node in path:
                try:
                    bin_idx, a, b = node
                    
                    results.append({
                        'CHROMOSOME': data['chromosomes'][bin_idx],
                        'START': data['starts'][bin_idx],
                        'END': data['ends'][bin_idx],
                        'CELL': cell,
                        'HAPLOTYPE_A': a,
                        'HAPLOTYPE_B': b,
                        'TOTAL_CN': a + b,
                        'BAF': data['baf_values'][bin_idx],
                        'RDR': data['rdr_values'][bin_idx]
                    })
                except Exception as e:
                    print(f"Error extracting results for node {node}: {e}")
                    continue
        
        # Save results
        if results:
            results_df = pd.DataFrame(results)
            results_df.to_csv(output_file, sep='\t', index=False)
            print(f"Results saved to {output_file}")
            return results_df
        else:
            print("No results generated. Check data and parameters.")
            return None
            
    except Exception as e:
        print(f"Error in haplotype inference: {e}")
        return None

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="Infer haplotype-specific copy numbers using path-generation method")
    parser.add_argument("--rdr", required=True, help="Path to RDR TSV file")
    parser.add_argument("--baf", required=True, help="Path to BAF TSV file with pre-computed BAF values")
    parser.add_argument("--output", required=True, help="Path to output TSV file")
    parser.add_argument("--baf-sigma", type=float, default=0.05, help="Standard deviation for BAF likelihood calculation (default: 0.05)")
    parser.add_argument("--max-cn", type=int, default=5, help="Maximum total copy number to consider")
    
    args = parser.parse_args()
    
    infer_haplotype_specific_copy_numbers(args.rdr, args.baf, args.output, args.max_cn, args.baf_sigma)