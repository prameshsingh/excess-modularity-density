# Excess Modularity Density

Computatioanl maximization of excess modularity density ($Q_x$) for weighted networks.


If you use this code, please cite the following preprint: 
Resolution limit revisited: community detection using generalized modularity density.
T. Chen, P. Singh, K. E. Bassler
Network community detection using modularity density measures,
*J. Stat. Mech.*, 053406 (2018)

To use the code, follow the steps below:

1. Compile the code using something like:


		`gcc -o a.out poc_Qx_matrix_byComp_wei.c Qx_noSP_wei.c -lm -O3`

	This will generate an executable `a.out`.

3. Run
	run it with 3 arguments, respectively as:
	* argument 1: number of nodes
	* argument 2: weight threshold as fraction of the maximum weight in the network.
	* argument 3: nput network as weighted adjacency matrix.
	* Example:

		`./a.out 100 0.1 example.txt`

4. Result

	* file: `Qx_t=0.1_partition_example`

	containing cluster label for each node in the order they appear in the adajcecny matrix.