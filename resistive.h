/* An attempt to reimplement Resistive with the ability to designate which nodes
 * to which voltage sources and resistors are connected. */
#ifndef RESISTIVE_H
#define RESISTIVE_H

#include <vector>

class
Circuit
{
	public:
		/* Class constructor */
		Circuit ();

		/* Add a voltage source.
		 * Currently one terminal has to be connected to ground */
		int addSource (unsigned int n1, unsigned int n2, double value);

		/* Add a resistor between nodes n1 and n2 */
		int addResistor (unsigned int n1, unsigned int n2, double value);

		/* Calculate the node voltages */
		int calcNodeVoltages ();

		/* Print the voltages at all the nodes */
		void printNodeVoltages ();

	private:
		unsigned int N; /* Number of nodes in the circuit */
		unsigned int Ns; /* Number of voltage sources */

		std::vector<std::vector<double>> G; /* Conductance matrix */

		std::vector<double> Vin; /* DC voltage sources */

		std::vector<unsigned int> inNode; /* Nodes to which the DC sources are connected */

		std::vector<double> vNode; /* Vector of node voltages */
};

#endif
