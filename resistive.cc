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
Circuit::Circuit (const unsigned int nodes, const double *cond, unsigned int sources, const double *V, const unsigned int *in)
{
	unsigned int i, j;

	/* Set the values of N, Vin, and inNode. Also set the size of vNode. */
	this->N = nodes;
	this->Ns = sources;
	this->Vin.resize (sources);
	this->inNode.resize (sources);
	this->vNode.resize (nodes);
	this->G.resize (nodes);

	/* Set the values of the DC sources */
	for (i = 0; i < sources; ++i) {
		Vin[i] = V[i];
		inNode[i] = in[i];
	}

	/* Create the conductance matrix and
	 * set the values in the matrix */
	for (i = 0; i < nodes; ++i) {
		this->G[i].resize (nodes);

		for (j = 0; j < nodes; ++j) {
			this->G[i][j] = cond[i * nodes + j];
			
			if (i != j) {
				if (this->G[i][j] > 0) {
					std::cerr << "G[" << i+1 << "][" << j+1 << "] is non-negative. Correcting." << std::endl;
					this->G[i][j] *= -1;
				}
			}
		}
	}
}

/* Create a resistance matrix. This function needs to be called before
 * using the function calcLoopCurrents.
 * loops is the number of loops in the circuit.
 * R is the resistance matrix
 * vLoop is a vector of the total voltage sources around the loop. In other words, it is the
 * right-hand side of the system of mesh equations. */
void
Circuit::createResistMatrix (const unsigned int loops, const double *R, const double *loopVolt)
{
	unsigned int i, j;

	this->Nl = loops;

	this->R.resize (loops);
	this->vLoop.resize (loops);
	this->iBranch.resize (loops);

	for (i = 0; i < loops; ++i) {
		this->R[i].resize (loops);

		/* Set the value for the voltage vector */
		this->vLoop[i] = loopVolt[i];

		/* Set the values in the resistance matrix */
		for (j = 0; j < loops; ++j) {
			this->R[i][j] = R[i * loops + j];
			
			if (i != j) {
				if (this->R[i][j] > 0) {
					std::cerr << "R[" << i+1 << "][" << j+1 << "] is non-negative. Correcting." << std::endl;
					this->R[i][j] *= -1;
				}
			}
		}
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
	unsigned int in;
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
		vVector[i] = 0;

		for (j = 0; j < nodes; ++j) {
			gMatrix[i * nodes + j] = -1 * this->G[i][j];
		}
	}

	/* Set the non-zero values in vVector from vIn */
	for (i = 0; i < this->Ns; ++i) {
		in = this->inNode[i] - 1;
		vVector[in] = this->Vin[i];

		for (j = 0; j < nodes; j++) {
			if (j == in) {
				gMatrix[in * nodes + j] = 1;
			}

			else {
				gMatrix[in * nodes + j] = 0;
			}
		}
	}

	/* Now that gMatrix and vVector are filled we can solve the
	 * system of equations for the node voltages. The result is stored
	 * in vVector. */
	err = LAPACKE_dgesv (LAPACK_ROW_MAJOR, nodes, 1, gMatrix, nodes, ipvt, vVector, 1);

	/* Check for errors */
	if (err == 0) {
		//this->vNode;
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

/* Return the voltage at node N */
double
Circuit::voltageAtNode (const unsigned int node)
{
	if (node <= N) {
		return this->vNode[node - 1];
	}

	else {
		std::cerr << "Invalid node" << std::endl;
		return 0;
	}
}

/* Print the node voltages */
void
Circuit::printNodeVoltages ()
{
	unsigned int i;

	for (i = 0; i < this->N; ++i) {
		std::cout << "V" << i+1 << " = " << this->vNode[i] << std::endl;
	}
}

/* Calculate the loop currents using mesh analysis. */
int
Circuit::calcLoopCurrents ()
{
	unsigned int i, j;
	double *resist;
	double *vVector;
	int ipvt[4];
	int err;

	/* Create a copy of the Resistance matrix since it will be overwritten
	 * in the function LAPACKE_dgesv below */
	resist = new double [this->Nl * this->Nl];
	if (resist == 0x00) {
		return -1;
	}
	
	/* Create a copy of the voltage vector since it will also be overwritten
	 * by LAPACKE_dgesv */
	vVector = new double [this->Nl];
	if (vVector == 0x00) {
		delete resist;
		return -1;
	}

	/* Copy the values */
	for (i = 0; i < this->Nl; ++i) {
		vVector[i] = this->vLoop[i];

		for (j = 0; j < this->Nl; ++j) {
			resist[i * this->Nl + j] = this->R[i][j];
		}
	}

	/* Solve the system of equations to find the loop currents */
	err = LAPACKE_dgesv (LAPACK_ROW_MAJOR, this->Nl, 1, resist, this->Nl, ipvt, vVector, 1);

	/* If no errors occured */
	if (err == 0) {
		for (i = 0; i < this->Nl; ++i) {
			this->iBranch[i] = vVector[i];
		}

		delete vVector;
		delete resist;

		return 0;
	}

	else if (err < 0) {
		std::cerr << "Value of argument " << -1 * err << " is illegal" << std::endl;

		delete vVector;
		delete resist;

		return err;
	}

	else {
		std::cerr << "Resistance matrix is singular" << std::endl;

		delete vVector;
		delete resist;

		return err;
	}
}

/* Calculate the current between nodes n1 and n2 using the calculated node voltages.
 * calcNodeVoltages must be called before calling this function. */
double
Circuit::calcBranchCurrent (const unsigned int n1, const unsigned int n2)
{
	double V1, V2;
	double G12;
	double I12;
	unsigned int i;
	
	/* Check if n1 == n2, in which case the current is zero */
	if (n1 == n2) {
		return 0;
	}
	
	/* If n1 == 0 or n2 == 0, then its node voltage is zero */
	if (n1 == 0) {
		V1 = 0;
		V2 = this->vNode[n2 - 1];
		
		/* To get the conductance between n1 and
		 * ground, start with the value along the
		 * diagonal in the G matrix  and add the other
		 * elements (we add since the off-diagonal elements
		 * are negative). */
		G12 = this->G[n2 - 1][n2 - 1];
		for (i = 0; i < this->N; ++i) {
			if (i == n2 - 1) {
				continue;
			}
			
			G12 += this->G[n2 - 1][i];
		}
	}
	
	else if (n2 == 0) {
		V1 = this->vNode[n1 - 1];
		V2 = 0;
		
		/* To get the conductance between n1 and
		 * ground, start with the value along the
		 * diagonal in the G matrix  and add the other
		 * elements (we add since the off-diagonal elements
		 * are negative). */
		G12 = this->G[n1 - 1][n1 - 1];
		for (i = 0; i < this->N; ++i) {
			if (i == n1 - 1) {
				continue;
			}

			G12 += this->G[n1 - 1][i];
		}
	}

	/* Otherwise get the voltages for nodes n1 and n2 from vNode */
	else {
		V1 = this->vNode[n1 - 1];
		V2 = this->vNode[n2 - 1];

		/* Get the resistance between nodes n1 and n2 from the G matrix */
		G12 = -1 * this->G[n1 - 1][n2 - 1];
	}

	/* Calculate I12 = (V1 - V2) * G12 */
	I12 = (V1 - V2) * G12;

	return I12;
}

/* Print the loop currents */
void
Circuit::printLoopCurrents ()
{
	unsigned int i;

	for (i = 0; i < this->Nl; ++i) {
		std::cout << "I" << i + 1 << " = " << this->iBranch[i] << std::endl;
	}

	std::cout << std::endl;
}

void
Circuit::printGMatrix ()
{
	unsigned int i, j;

	for (i = 0; i < this->N; ++i) {
		for (j = 0; j < this->N; ++j) {
			std::cout << this->G[i][j] << "\t";
		}

		std::cout << std::endl;
	}
}
