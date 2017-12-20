# Co-occurrence network analysis
This folder contains R and python scripts for correlation-based network analysis, starting from an abundance table of microbial entities (e.g., OTUs, genus and genes). The abundance table is a tab-delimited text file in which each row represents a microbial entity and each column represents a sampling point.

Note: a minimum number of 20 to 25 samples is suggested for achieving sufficient specificity of co-occurrence networks (reference: doi: 10.3389/fmicb.2014.00219)


Functions

1.Pairwise_correlations.R

This scrpit defines the function "co_occurrence_network" and proceeds in three major steps:
a) calculates all pairwise Spearman's correlations between abundance of microbial entities (e.g., OTUs, gene, transcripts).
b) filtrate the correlations by user defined cutoffs for coefficient (alpha) and FDR-adjusted P-value (p.cutoff)
c) generate a gml-formatted network file from filtered correlations

2.Network_Analysis.R
This scrpit reads the abundance table of microbial entities, generates gml-formatted network files, and calculates topological properties of the observed co-occurrence network:
a) Filtrate OTUs by occurrence frequency (i.e.,number of samples an OTU is Present)
b) Generatte gml files of co-occurrence network, which can be visulized in Gephi (https://gephi.org/)
c) Calculating network topological properties, such as number of edges (e), number of nodes (v), clustering coefficient (cc), short path length (spl), modularity (md),network diameter (nd), graph.densit (GD), etc. These properties can be further used to explore the 'small-world' properties in the observed co-ocurrence network, as shown in the two references (see the end of this file).

3. Random_vs_observed_cooccurrence.py
This script uses a map file and a gml file as input to calculated the random and observed incidences of co-occurrence patterns between microbial entities. The random incidence is calculated by considering frequency (number) of each type of entities and assuming random connections between any two nodes. The observed incidence is calculted by observed number of edges (connections) dividied by the total number of edges in a co-occurrence network. Details on the calculation methods for random vs. observed incidences of co-occurrence can be found in the methods of Ju et al., The ISME Journal (2015) 9, 683–695 (2015)

The map file is a tab-delimited file with node ID (col 1, e.g., OTU name/ID) and type of node (col 2, e.g., genus/family names of each node).
The gml file is a tab-delimited file is the gml file generated from 2.Network_Analysis.R

4.Random_network_simulator.R
This script generates a large number of random networks based on Erdős–Rényi model and calculate typical topological properties of random networks to enable their comparsion with a real network.

Users can set the size of random network (i.e., number of network nodes and edges) and the number of random networks to create.

References: 
1. Ju F, Xia Y, Guo F, Wang ZP, Zhang T. 2014. Taxonomic relatedness shapes bacterial assembly in activated sludge of globally distributed wastewater treatment plants. Environmental Microbiology. 16(8):2421-2432
2. Ju F, Zhang T. 2015. Bacterial assembly and temporal dynamics in activated sludge of a full-scale municipal wastewater treatment plant. The ISME Journal. 9: 683-695 
3. Hu AY, Ju F, Hou LY, Li JW, Yang XY, Wang HJ, Mulla SI, Sun Q, Bürgmann H, Yu CP. 2017. Strong impact of anthropogenic contamination on the co-occurrence patterns of a riverine microbial community. Environmental Microbiology (2017) 19(12), 4993–500


