/* Example circuit with three nodes. */
#include "resistive.h"

int
main ()
{
	unsigned int nodes = 3; /* Number of nodes */
	double g[9] = {1, -1, 0, -1, 3, -1, 0, -1, 2}; /* Conductance matrix */
	double vIn = 5; /* Input source voltage */
	unsigned int in = 1; /* Source is connected to node 1 */

	/* Create a DC circuit */
	Circuit testCircuit (0);

	testCircuit.addSource (1, 0, 5);
	testCircuit.addResistor (1, 2, 1);
	testCircuit.addResistor (2, 0, 1);
	testCircuit.addResistor (2, 3, 1);
	testCircuit.addResistor (3, 0, 1);

	/* Calculate the node voltages */
	testCircuit.calcNodeVoltages ();

	/* Print the node voltages */
	testCircuit.printNodeVoltages ();

	return 0;
}
