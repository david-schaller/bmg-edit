# -*- coding: utf-8 -*-

"""
Extended Best Hits Method.

Implementation of the Extended Best Hits method for best match inference.
"""

import os, subprocess, time

from asymmetree.analysis.BestMatches import extended_best_hits
from asymmetree.file_io.ScenarioFileIO import (parse_bmg_edges,
                                               matrix_to_phylip,
                                               species_to_genes)


__author__ = 'David Schaller'


# --------------------------------------------------------------------------
#                           PYTHON IMPLEMENTATION
#
# --------------------------------------------------------------------------
    
def ebh(leaves, D, epsilon=1e-8):
    """Compute BMG and RBMG from a distances matrix D.
    
    Keyword arguments:
        epsilon -- epsilon for relative BM threshold: (x,y) in BMG if
                   D(x,y) <= (1+epsilon) * min d(x,y'),
                   default=1e-8 (for limited float precision).
    """
    
    return extended_best_hits(leaves, D, epsilon=epsilon, supply_rbmg=True)


# --------------------------------------------------------------------------
#                            EXTERNAL C++ PROGRAM
#
# --------------------------------------------------------------------------
 
def ebh_qinfer(scenario,
               matrix_filename, species_filename,
               epsilon=0.5,
               benchmark_file=None,
               binary_path=None):
    """Compute BMG and RBMG from a distances matrix D using 'qinfer'.
    
    Keyword arguments:
        epsilon -- epsilon for relative BM threshold: (x,y) in BMG if
                   D(x,y) <= (1+epsilon) * min d(x,y'),
                   default=10E-8 (for limited float precision).
        benchmark_file -- activate benchmarking in 'qinfer' and
                          specify the filename
        binary_path -- path to 'qinfer' binary (if not available
                       within path)
    """
    
    if not binary_path:
        qinfer_command = "qinfer"
    elif os.path.exists(binary_path):
        qinfer_command = binary_path
    else:
        raise FileNotFoundError("path to qinfer binary file '{}' does not exist".format(binary_path))
    
    output = -1
    command = [qinfer_command, matrix_filename, species_filename,
               "--disable-quartet", "--epsilon=" + str(epsilon)]

    if benchmark_file is not None:
        command.append( "--benchmark=" + benchmark_file )
    
    # call 'qinfer' and measure execution time
    start = time.time()
    
    try:
        output = subprocess.run(command, stdout=subprocess.PIPE)
    except:
        raise FileNotFoundError("calling qinfer failed")
    
    exec_time = time.time() - start
    
    if output == -1:
        raise Exception("no output from qinfer")
    
    bmg = parse_bmg_edges(output.stdout.decode(), scenario)
    
    return bmg, bmg.to_undirected(reciprocal=True), exec_time


def ebh_from_scenario(scenario, epsilon=0.5):
    """Compute BMG and RBMG from a scenario using 'qinfer'.
    
    Keyword arguments:
        epsilon -- epsilon for relative BM threshold: (x,y) in BMG if
                   D(x,y) <= (1+epsilon) * min d(x,y'),
                   default=10E-8 (for limited float precision).
    """

    matrix_filename = "temp.phylip"
    species_filename = "temp_species.txt"
    
    matrix = scenario.get_distance_matrix()
    matrix_to_phylip(matrix_filename, scenario.genes, matrix)
    species_to_genes(species_filename, scenario)
    
    bmg, rbmg, exec_time = ebh_qinfer(scenario,
                                      matrix_filename, species_filename,
                                      epsilon=epsilon)
    
    os.remove(matrix_filename)
    os.remove(species_filename)
    
    return bmg, rbmg, exec_time