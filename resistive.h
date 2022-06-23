/* resistive.h
 * Type definitions for Resistive. */
#ifndef RESISTIVE_H
#define RESISTIVE_H

/* Type for Resistive circuits.
 * Each instance has the number of nodes, a matrix
 * of the resistors in the circuit, and a DC
 * voltage source.
 *
 * R matrix has the form:
 *     | R_11 R_12 ... R_1N R_10 |
 * R = | R_21 R_22 ... R_2N R_20 |
 *     | ...			 |
 *     | R_N1 R_N2 ... R_NN R_N0 |
 *
 * where R_ij is the resistance between nodes i and j. If i < j then
 * R_ij is positive. If i > j then R_ij is negative.
 *
 * Currently only one DC input source is supported.
 */

#include <vector>

class
Circuit
{
	public:
		/* Class constructor */
		Circuit (unsigned int nodes, double *resis, double V);

		/* Used to calculate the node voltages */
		friend void calcNodeVoltages (Circuit *c);

	private:
		unsigned int N; /* Number of nodes in the circuit, including ground */
		std::vector<std::vector<double>> R; /* Resistance matrix */
		double Vin; /* DC voltage source */
		unsigned int inNode; /* Node where the input source is connected */

		std::vector<double> vNode; /* Vector of node voltages */
};

#endif
