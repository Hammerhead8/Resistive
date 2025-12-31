/* An attempt to reimplement Resistive with the ability to designate which nodes
 * to which voltage sources and resistors are connected so the conductance
 * matrix and source vector don't need to be calculated by the user.
 * Also use exceptions for error handling. */
#ifndef RESISTIVE_H
#define RESISTIVE_H

#include <vector>
#include <complex>

class
Circuit
{
	public:
		/* Class constructor */
		Circuit (double freq);
		
		Circuit (int lowerExp, int upperExp, int numPoints);

		/* Add a voltage source.
		 * Currently one terminal has to be connected to ground */
		void addSource (unsigned int n1, unsigned int n2, double value);

		/* Add a resistor between nodes n1 and n2 */
		void addResistor (unsigned int n1, unsigned int n2, double value);
		
		/* Add an inductor between nodes n1 and n2 */
		void addInductor (unsigned int n1, unsigned int n2, double value);
		
		/* Add a capacitor between nodes n1 and n2 */
		void addCapacitor (unsigned int n1, unsigned int n2, double value);

		/* Calculate the node voltages */
		int calcNodeVoltages ();

		/* Print the voltages at all the nodes */
		void printNodeVoltages (int node = 0);

	private:
		/* Calculate node voltages for DC circuits */
		int calcDCNodes ();
		
		/* Calculate node voltages for AC circuits */
		int calcACNodes ();
		
		unsigned int N; /* Number of nodes in the circuit */
		unsigned int Ns; /* Number of voltage sources */
		std::vector<double> w; /* Frequency of the circuit in rad/s */

//		std::vector<std::vector<double>> G; /* Conductance matrix */
		std::vector<std::vector<std::vector<std::complex<double>>>> G; /* Conductance matrix */

		std::vector<double> Vin; /* DC voltage sources */

		std::vector<unsigned int> inNode; /* Nodes to which the DC sources are connected */

//		std::vector<double> vNode; /* Vector of node voltages */
		std::vector<std::vector<std::complex<double>>> vNode; /* Vector of node voltages */

//		std::vector<double> iBranch; /* Vector of loop currents */
};

#endif
