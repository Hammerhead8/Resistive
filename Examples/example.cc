#include "../resistive.h"

int
main ()
{
	Circuit c (0);

	c.addSource (1, 0, 5);

	c.addResistor (1, 2, 1);
	c.addResistor (2, 0, 1);
	c.addResistor (2, 3, 1);
	c.addResistor (3, 0, 1);

	c.calcNodeVoltages ();
	c.printNodeVoltages ();

	return 0;
}
