/* resistive.cc
 * Function definitions for Resistive.
 *
 * Analysis is done by nodal analysis
 * using a conductance matrix */

#include "resistive.h"
#include <iostream>
#include <cmath> /* For nan */
#include <lapacke.h> /* For dgesv */

/* Class constructor.
 * Set the values of N and Vin as well as create the resistance matrix */
Circuit::Circuit (const unsigned int nodes, const double *cond, const double V, const unsigned int in)
{
	unsigned int i, j;

	/* Set the values of N, Vin, and inNode. Also set the size of vNode. */
	this->N = nodes;
	this->Vin = V;
	this->inNode = in;
	this->vNode.resize (nodes);
	this->G.resize (nodes);

	/* Create the resistance matrix */
	for (i = 0; i < nodes; ++i) {
		this->G[i].resize (nodes);
	}

	/* Now set the values in the matrix */
	for (i = 0; i < nodes; ++i) {
		for (j = 0; j < nodes; ++j) {
			this->G[i][j] = cond[i * nodes + j];
		}
	}
}

void
Circuit::printNodeVoltages ()
{
	unsigned int i;

	for (i = 0; i < this->N; ++i) {
		std::cout << "V" << i+1 << " = " << this->vNode[i] << std::endl;
	}
}

/* Calculate the node voltages using modified nodal analysis.
 * Returns 0 if the node voltages are successfully found, -1 if
 * either memory allocation failed or there is an illegal argument
 * to LAPACKE_dgesv, or 1 if the conductance matrix is signular. */
int
Circuit::calcNodeVoltages ()
{
	unsigned int i, j; /* Loop counters */
	const unsigned int nodes = this->N;
	double *gMatrix; /* Conductance matrix for the system of equations */
	double *vVector; /* Vector for node voltages */
	unsigned int gSize; /* Number of rows and columns in rMatrix */
	unsigned int vSize; /* Number of elements in vVector */
	lapack_int ipvt[this->N]; /* Used by LAPACKE_dgesv */
	lapack_int err; /* Returned value from LAPACKE_dgesv */

	/* Create the system of equations from the node equations */
	gSize = nodes * nodes;
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

	/* Fill gMatrix and vVector */
	for (i = 0; i < nodes; ++i) {
		/* If the voltage source is applied to node i, then
		 * we need to set G_ij = 1 and the other elements in
		 * the row to 0. */
		if (this->inNode == i + 1) {
			vVector[i] = this->Vin;

			for (j = 0; j < nodes; ++j) {
				if (j == i) {
					gMatrix[i * nodes + j] = 1;
				}
				else {
					gMatrix[i * nodes + j] = 0;
				}
			}

		}

		/* Otherwise the G matrix is just copied. */
		else {
			vVector[i] = 0;

			for (j = 0; j < nodes; ++j) {
				gMatrix[i * nodes + j] = -1 * this->G[i][j];
			}
		}
	}

	/* Now that gMatrix and vVector are filled we can solve the
	 * system of equations for the node voltages. The result is stored
	 * in vVector. */
	err = LAPACKE_dgesv (LAPACK_ROW_MAJOR, nodes, 1, gMatrix, nodes, ipvt, vVector, 1);

	/* Check for errors */
	if (err == 0) {
		this->vNode;
		for (i = 0; i < this->N; ++i) {
			this->vNode[i] = vVector[i];
		}
	}

	else if (err < 0) {
		fprintf (stderr, "Argument %d is invalid.\n", -1 * err);
		for (i = 0; i < this->N; ++i) {
			this->vNode[i] = NAN;
		}
	}

	else {
		fprintf (stderr, "G matrix is singular.\n");
		for (i = 0; i < this->N; ++i) {
			this->vNode[i] = NAN;
		}
	}

	/* The last step is to delete gMatrix and vVector */
	delete[] gMatrix;
	delete[] vVector;

	/* Return 0 if success, -1 if illegal argument, and 1 for singular matrix */
	return err;
}
