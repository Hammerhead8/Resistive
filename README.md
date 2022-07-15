# Resistive
Resistive is a C++ library for calculating the node voltages and branch currents in purely resistive cicuits.
The user provides the number of nodes in the circuit, a conductance matrix, and the value and location of a DC
input voltage source. Resistive using nodal analysis to calculate the node voltages. Using these, it can calculate
the branch currents.

It is also possible to calculate the loop currents using mesh analysis. The user provides a resistance matrix and a vector of the loop voltages, which will be illustrated in the example below.

Currently, Resistive supports multiple DC voltage sources, but only if the sources have a terminal connected to ground.

# Dependencies
The only external dependancy is LAPACKE. If you're using Linux then it should be available through your software
repositories. Otherwise, it can be downloaded and built from https://netlib.org/lapack.

# Examples
The usage of Resistive can be shown using a basic example. consider the following circuit:

![example circuit](/Examples/example.svg)

Node 1 is between Vin and R1, node 2 is between R1 and R3, and node 3 is between R3 and R4. The node voltages can
be easily calculated by hand as V1 = 5V, V2 = 2V, and V3 = 1V. The loop currents can also be easily calculated as I1 = 3A and I2 = 1A using the convention that positive current flows clockwise around the loop. To solve for these voltages using Resistive, we need
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

To calculate the loop currents as well, we can construct a resistance matrix as
```
    | R1 + R2       -R2      |
R = |   -R2     R2 + R3 + R4 |
```

and the loop voltage vector as
```
    | Vin |
V = |  0  |
```

Substituting the values of R1, R2, and Vin gives
```
    | 2 -1 |
R = | -1 3 |

    | 5 |
V = | 0 |
```

Using G as well as Vin = 5 and the fact that Vin is connected at node 1, we can call Resistive using
```
#include "resistive.h"

int
main ()
{
  unsigned int nodes = 3;
  unsigned int loops = 2;
  double G[9] = {1, -1, 0, -1, 3, -1, 0, -1, 2};
  double R[4] = {2, -1, -1, 3};
  unsigned int sources = 1;
  double vIn[1] = {5};
  double vLoop[2] = {5, 0};
  unsigned int in = 1;

  /* Create a circuit using the above values */
  Circuit c (nodes, g, sources, vIn, in);

  /* Create the resistance matrix for calculating the loop currents */
  c.createResistMatrix (loops, R, vLoop);

  /* Calculate the node voltages of the circuit */
  c.calcNodeVoltages ();

  /* Calculate the loop currents */
  c.calcLoopCurrents ();

  /* Print the node voltages */
  c.printNodeVoltages ();

  /* Print the loop currents */
  c.printLoopCurrents ();
}
```

This program creates the following output:
```
V1 = 5
V2 = 2
V3 = 1
I1 = 3
I2 = 1
```

Which matches what was calculated by hand.
