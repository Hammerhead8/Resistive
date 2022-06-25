/* resistive.h
 * Type definitions for Resistive. */
#ifndef RESISTIVE_H
#define RESISTIVE_H

/* Type for Resistive circuits.
 * Each instance has the number of nodes, a matrix
 * of the conductance between nodes in the circuit, and a DC
 * voltage source.
 *
 * G matrix has the form:
 *     | G_11 G_12 ... G_1N |
 * G = | G_21 G_22 ... G_2N |
 *     | ...  ...  ... ...  |
 *     | G_N1 G_N2 ... G_NN |
 *
 * where G_ij is the conductance beween nodes i and j, where G == 1 / R.
 * The matrix is constructed the same was as an admittance matrix used
 * for power flow analysis. G_ii = g_i1 + g_i2 + ... + g_iN and
 * G_ij = -g_ij. g_ij = 0 if there is no connection between nodes
 * i and j.
 *
 * Currently only one DC input source is supported.
 */

#include <vector>

class
Circuit
{
	public:
		/* Class constructor */
		Circuit (const unsigned int nodes, const double *cond, const double V, const unsigned int in);

		/* Print the node voltages */
		void printNodeVoltages ();

		/* Used to calculate the node voltages */
		int calcNodeVoltages ();

	private:
		unsigned int N; /* Number of nodes in the circuit, including ground */
		std::vector<std::vector<double>> G; /* Conductance matrix */
		double Vin; /* DC voltage source */
		unsigned int inNode; /* Node where the input source is connected */

		std::vector<double> vNode; /* Vector of node voltages */
};

#endif
