/* Function definitions for the reimplementation of Resistive.
 *
 * Node voltages are found using nodal analysis with a conductance matrix */

#include "resistive2.h"
#include <iostream>
#include <string>
#include <exception>
#include <cmath>
#include <lapacke.h>

/* Return the larger of two integers */
#define max(x, y) (x > y ? x : y)

/* Class constructor.
 * An empty circuit is created and all variables are initialized. */
Circuit::Circuit ()
{
	/* Set the values of N and Ns to zero */
	this->N = 0;
	this->Ns = 0;

	/* Create the vectors Vin and inNode */
	this->Vin.resize (0);
	this->inNode.resize (0);

	/* Create the conductance matrix */
	this->G.resize (1);
	this->G[0].resize (1);
}

/* Add a DC voltage source */
int
Circuit::addSource (unsigned int n1, unsigned int n2, double value)
{
	/* If neither of the nodes is zero then we cannot proceed */
	if (n1 != 0 && n2 != 0) {
//		std::cerr << "Voltage source must have one node connected to ground." << std::endl;
		throw std::runtime_error ("Voltage source must have one node connected to ground.");
	}

	/* Make n1 the larger node and n2 the ground */
	if (n1 == 0) {
		n1 = n2;
		n2 = 0;
		value *= -1;
	}

	/* Increment the value of Ns, which is the number of sources in the circuit */
	this->Ns++;

	/* Now add the source to the Vin vector and update the inNode vector */
	this->Vin.resize (this->Ns);

	/* Resize inNode */
	this->inNode.resize (this->Ns);

	this->Vin[this->Ns-1] = value;
	this->inNode[this->Ns-1] = n1;

	return 0;
}

/* Add a resistor to the circuit and update the conductance matrix */
int
Circuit::addResistor (unsigned int n1, unsigned int n2, double value)
{
	unsigned int i;
	unsigned int temp;

	/* Check if n1 == n2. If so then throw an exception. */
	if (n1 == n2) {
//		std::cerr << "Positive and negative terminal of a resistor must be connected to different nodes." << std::endl;
		throw std::runtime_error ("Positive and negative terminal of a resistor must be connected to different nodes.");
	}

	/* If the resistance is negative then we cannot proceed */
	if (value <= 0) {
//		std::cerr << "Zero or negative resistance between nodes " << n1 << " and " << n2 << " is not valid." << std::endl;
		throw std::runtime_error (std::string ("Zero or negative resistance between nodes " + std::to_string (n1) + " and " + std::to_string (n2)));
	}

	/* Now update the conductance matrix.
	 * First check if the larger node is part of the circuit.
	 * If not then resize the conductance matrix so it is.
	 * If it is then the first node is also in the current circuit. */
	if ((n1 > this->N) || (n2 > this->N)) {
		this->N = max (n1, n2);
	}
	
	if (this->N > this->G.size ()) {
		this->G.resize (this->N);
		
		for (i = 0; i < this->G.size (); ++i) {
			this->G[i].resize (this->N);
		}
	}

	/* Add the new resistor to the conductance matrix.
	 * If one of the nodes is 0, meaning ground, then we
	 * only need to add it to its corresponding diagonal element. */
	if (n2 == 0) {
		this->G[n1-1][n1-1] += (1 / value);
	}

	else if (n1 == 0) {
		this->G[n2-1][n2-1] += (1 / value);
	}

	/* Otherwise we need to add it to multiple diagonals and two
	 * off-diagonal elements. */
	else {
		/* First add it to the diagonal elements */
		this->G[n1-1][n1-1] += (1 / value);
		this->G[n2-1][n2-1] += (1 / value);

		/* Now subtract it from the off-diagonal elements */
		this->G[n1-1][n2-1] -= (1 / value);
		this->G[n2-1][n1-1] -= (1 / value);
	}
	 
	return 0;
}


/* Calculate the node voltages */
int
Circuit::calcNodeVoltages ()
{
	unsigned int i, j; /* Loop counters */
	unsigned int in;
	unsigned int gSize;
	unsigned int vSize;
	unsigned int nodes = this->N; /* Number of nodes */
	unsigned int sources = this->Ns; /* Number of sources */
	double *gMatrix; /* Copy of the conductance matrix */
	double *vVector; /* Copy of the source vector */
	lapack_int ipvt[nodes]; /* Needed by LAPACKE_dgesv to find the node voltages */
	lapack_int err; /* Returned value from LAPACKE_dgesv */

	/* Create copies of the vectors to preserve their values */
	gSize  = nodes * nodes;
	vSize = nodes;

	gMatrix = new double [gSize];
	if (gMatrix == 0x00) {
		return -1;
	}

	vVector = new double [vSize];
	if (vVector == 0x00) {
		delete[] gMatrix;
		return -1;
	}

	/* Now fill gMatrix and vVector */
	for (i = 0; i < nodes; ++i) {
		vVector[i] = 0;

		for (j = 0; j < nodes; ++j) {
			gMatrix[i * nodes + j] = this->G[i][j];
		}
	}

	/* Set the non-zero values in vVector from vIn */
	for (i = 0; i < sources; ++i) {
		in = this->inNode[i] - 1;
		vVector[in] = this->Vin[i];

		/* Set the diagonal element to 1 and the others to zero
		 * in the conductance matrix for the row corresponding to
		 * the source */
		for (j = 0; j < nodes; ++j) {
			if (j == in) {
				gMatrix[in * nodes + j] = 1;
			}

			else {
				gMatrix[in * nodes + j] = 0;
			}
		}
	}

	/* Allocate vNode */
	this->vNode.resize (nodes);

	/* Now that gMatrix and vVector are filled we can calculate the nodes voltages
	 * by solving the system of equations. The result is stored in vVector. */
	err = LAPACKE_dgesv (LAPACK_ROW_MAJOR, nodes, 1, gMatrix, nodes, ipvt, vVector, 1);

	/* Check for errors */
	if (err == 0) {
		/* If no errors occured then copy the results from vVector to vNode */
		for (i = 0; i < nodes; ++i) {
			this->vNode[i] = vVector[i];
		}
	}

	else if (err < 0) {
		std::cerr << "Argument " << -1 * err << " is invalid." << std::endl;
		for (i = 0; i < nodes; ++i) {
			this->vNode[i] = NAN;
		}
	}

	else {
		std::cerr << "Conductance matrix is singular." << std::endl;
		for (i = 0; i < nodes; ++i) {
			this->vNode[i] = NAN;
		}
	}

	/* The last step is to delete gMatrix and vVector */
	delete[] gMatrix;
	delete[] vVector;

	/* Return the status of LAPACKE_dgesv */
	return err;
}

/* Print the voltages at all nodes */
void
Circuit::printNodeVoltages ()
{
	unsigned int i;

	for (i = 0; i < this->N; ++i) {
		std::cout << "V" << i + 1 << " = " << this->vNode[i] << std::endl;
	}
}

