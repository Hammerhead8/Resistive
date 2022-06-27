# Resistive
Resistive is a C++ library for calculating the node voltages and branch currents in purely resistive cicuits.
The user provides the number of nodes in the circuit, a conductance matrix, and the value and location of a DC
input voltage source. Resistive using nodal analysis to calculate the node voltages. Using these, it can calculate
the branch currents.

Currently, Resistive only supports circuits with one DC voltage source.

# Dependencies
The only external dependancy is LAPACKE. If you're using Linux then it should be available through your software
repositories. Otherwise, it can be downloaded and built from https://netlib.org/lapack.

# Examples
The usage of Resistive can be shown using a basic example. consider the following circuit:

![example circuit](/Examples/example.svg)

Node 1 is between Vin and R1, node 2 is between R1 and R3, and node 3 is between R3 and R4. The node voltages can
be easily calculated by hand as V1 = 5, V2 = 2, and V3 = 1. To solve for these voltages using Resistive, we need
a [conductance matrix](https://en.wikipedia.org/wiki/Nodal_admittance_matrix#Construction). This is constructed the same way as an admittance matrix when performing power flow analysis.
Along the diagonal is the sum of the conductances of all resistors connected to that node, and across the rows are
the conductances between nodes i and j. For this example circuit, the conductance matrix would be
```
    |  1/R1        -1/R1                 0     |
G = | -1/R1  1/R1 + 1/R2 + 1/R3       -1/R3    |
    | 0            -1/R3           1/R3 + 1/R4 |
```
Substituting the values for the resistors gives
```
    |  1  -1   0 |
G = | -1   3  -1 |
    |  0  -1   2 |
```

Using G as well as Vin = 5 and the fact that Vin is connected at node 1, we can call Resistive using
```
#include <resistive.h>

int
main ()
{
  unsigned int nodes = 3;
  double G[9] = {1, -1, 0, -1, 3, -1, 0, -1, 2};
  double vIn = 5;
  unsigned int in = 1;

  /* Create a circuit using the above values */
  Circuit c (nodes, g, vIn, in);

  /* Calculate the node voltages of the circuit */
  c.calcNodeVoltages ();

  /* Print the node voltages */
  c.printNodeVoltages ();
}
```

This program creates the following output:
```
V1 = 5
V2 = 2
V3 = 1
```

Which matches what was calculated by hand.
