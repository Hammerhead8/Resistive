# Resistive
Resistive is a C++ tool for calculating the node voltages in DC and AC circuits. The syntax is similar to SPICE, where the user provides the nodes where the element is connected and its value. This means that the user does not need to calculate the conductance matrix themselves. Instead, the conductance matrix is updated automatically as elements are added to the circuit. Resistive uses nodal analysis to calculate the node voltages.

Currently, Resistive supports DC and AC sources, but only if all sources are of the same frequency and all sources need to have one node connected to ground.

# Features
* Support for multiple voltages sources
* Support for both DC and AC circuits
* Support for AC sweep analysis with node voltages printed to a text file for analysis
* Able to print voltages for either the whole circuit or a particular node

# Limitations
Currently Resistive has the following limitations. Note that some of these are no plans to add support for some of these since they are outside the intended scope of the tool.
* No checks for invalid connections such as components not connected to the rest of the circuit
* All voltages sources must have ground as one of their nodes
* No support for Monte Carlo simulation using random values for a particular component
* No support for transient analysis
* No built-in support for plotting result of AC sweep

# Dependencies
The only external dependancy is LAPACKE. If you're using Linux then it should be available through your software
repositories. Otherwise, it can be downloaded and built from https://netlib.org/lapack.

# Examples
The usage of Resistive can be shown using a basic example. Consider the following circuit:
![example circuit](/Examples/example.svg)

Node 1 is between Vin and R1, node 2 is between R1 and R2, and node 3 is between R3 and R4. The node voltages can be easily calculated by hand as V1 = 5V, V2 = 2V, and V3 = 1V. To solve this using Resistive, we use the following program:

```
#include "resistive.h"

int
main ()
{
    /* Create an empty DC circuit */
    Circuit c (0);

    /* Add the voltage source */
    c.addSource (1, 0, 12);

    /* Now add the resistors */
    c.addResistor (1, 2, 1);
    c.addResistor (2, 0, 1);
    c.addResistor (2, 3, 1);
    c.addResistor (3, 0, 1);

    /* Calculate the node voltages */
    c.calcNodeVoltages ();

    /* Print the node voltages */
    c.printNodeVoltages ();

    return 0;
}
```

This program creates the following output:
```
V1 = 5
V2 = 2
V3 = 1
```

Which matches the expected solution. Note that the solutions are given in phasor form, where the phase shift is given in degrees.

Another example using an RLC circuit is given below. Using the same circuit as before, replace R2 with a 1 mH inductor and R3 with a 1 uF capacitor. Letting the input frequency be 60 Hz (377 rad/s), we can write the following program:

```
#include "resistive.h"

int
main ()
{
    /* Create an empty circuit with input frequency of 377 rad/s */
    Circuit c (377);

    /* Add the voltage source */
    c.addSource (1, 0, 5);

    /* Now add the resistors, inductor, and capacitor */
    c.addResistor (1, 2, 1);
    c.addInductor (2, 0, 1e-3);
    c.addResistor (2, 3, 1);
    c.addCapacitor (3, 0, 1e-6);

    /* Calculate the node voltages */
    c.calcNodeVoltages ();

    /* Print the node voltages */
    c.printNodeVoltages ();

    return 0;
}
```

This program creates the following output:
```
V1 = 5<0
V2 = 1.50497<52.9876
V3 = .000567373<217.034
```

# Advanced usage
Resistive can throw exceptions when adding elements to the circuit. If an invalid connection is made (connecting both terminals to the same node, negative resistance, etc.) then a runtime exception is thrown and its message can be printed. Using the DC circuit above, this can be added as follows:

```
#include "resistive.h"

int
main ()
{
    /* Create an empty DC circuit */
    Circuit c (0);

    /* Catch any exceptions which may occur when adding elements */
    try {
        /* Add the voltage source */
        c.addSource (1, 0, 12);

        /* Now add the resistors */
        c.addResistor (1, 2, 1);
        c.addResistor (2, 0, 1);
        c.addResistor (2, 3, 1);
        c.addResistor (3, 0, 1);
    }

    catch (std::exception& e) {
        std::cerr << e.what () << std::endl;
    }

    /* Calculate the node voltages */
    c.calcNodeVoltages ();

    /* Print the node voltages */
    c.printNodeVoltages ();

    return 0;
}
```

The same is done when adding inductor and capacitors.
