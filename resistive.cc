/* Function definitions for the reimplementation of Resistive.
 *
 * Node voltages are found using nodal analysis with a conductance matrix */

#include "resistive.h"
#include <iostream>
#include <exception>
#include <cmath>
#include <lapacke.h>

/* Return the larger of two integers */
#define max(x, y) (x > y ? x : y)

/* Class constructor.
 * An empty circuit is created and all variables are initialized. */
Circuit::Circuit (double freq)
{
	/* Set the values of N and Ns to zero and w to the circuit frequency */
	this->N = 0;
	this->Ns = 0;
	this->w = freq;

	/* Create the vectors Vin and inNode */
	this->Vin.resize (0);
	this->inNode.resize (0);

	/* Create the conductance matrix */
	this->G.resize (1);
	this->G[0].resize (1);
}

/* Add a DC voltage source */
void
Circuit::addSource (unsigned int n1, unsigned int n2, double value)
{
	/* If neither of the nodes is zero then we cannot proceed */
	if (n1 != 0 && n2 != 0) {
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
}

/* Add a resistor to the circuit and update the conductance matrix */
void
Circuit::addResistor (unsigned int n1, unsigned int n2, double value)
{
	unsigned int i;
	unsigned int temp;

	/* Check if n1 == n2. If so then throw an exception. */
	if (n1 == n2) {
		throw std::runtime_error ("Positive and negative terminal of a resistor must be connected to different nodes.");
	}

	/* If the resistance is negative then we cannot proceed so throw an exception */
	if (value <= 0) {
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
}

/* Add an inductor between nodes n1 and n2 */
void
Circuit::addInductor (unsigned int n1, unsigned int n2, double value)
{
	unsigned int i;
	
	/* Check if the circuit frequency is 0 (DC).
	 * If so we can't add the inductor */
	if (this->w == 0) {
		throw std::runtime_error ("Cannot add an inductor to a DC circuit.");
	}
	
	/* Check if n1 == n2. If so then throw an exception. */
	else if (n1 == n2) {
		throw std::runtime_error ("Positive and negative terminal of an inductor must be connected to different nodes.");
	}

	/* If the resistance is negative then we cannot proceed so throw an exception */
	else if (value <= 0) {
		throw std::runtime_error (std::string ("Zero or negative inductance between nodes " + std::to_string (n1) + " and " + std::to_string (n2)));
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
	
	/* Calculate the impedance of the inductor */
	std::complex Z (0.0, this->w * value);

	/* Add the new resistor to the conductance matrix.
	 * If one of the nodes is 0, meaning ground, then we
	 * only need to add it to its corresponding diagonal element. */
	if (n2 == 0) {
		this->G[n1-1][n1-1] += (1.0 / Z);
	}

	else if (n1 == 0) {
		this->G[n2-1][n2-1] += (1.0 / Z);
	}

	/* Otherwise we need to add it to multiple diagonals and two
	 * off-diagonal elements. */
	else {
		/* First add it to the diagonal elements */
		this->G[n1-1][n1-1] += (1.0 / Z);
		this->G[n2-1][n2-1] += (1.0 / Z);

		/* Now subtract it from the off-diagonal elements */
		this->G[n1-1][n2-1] -= (1.0 / Z);
		this->G[n2-1][n1-1] -= (1.0 / Z);
	}
}

/* Add a capacitor between nodes n1 and n2 */
void
Circuit::addCapacitor (unsigned int n1, unsigned int n2, double value)
{
	unsigned int i;
	
	/* Check if the circuit frequency is 0 (DC).
	 * If it is then we can't add the capacitor so throw an exception. */
	if (this->w == 0) {
		throw std::runtime_error ("Cannot add a capacitor to a DC circuit.");
	}
	
	/* Check if n1 == n2. If so then throw an exception. */
	else if (n1 == n2) {
		throw std::runtime_error ("Positive and negative terminal of a capacitor must be connected to different nodes.");
	}

	/* If the resistance is negative then we cannot proceed so throw an exception */
	else if (value <= 0) {
		throw std::runtime_error (std::string ("Zero or negative capacitance between nodes " + std::to_string (n1) + " and " + std::to_string (n2)));
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
	
	/* Calculate the impedance of the capacitor */
	std::complex Z (0.0, -1 / (this->w * value));

	/* Add the new resistor to the conductance matrix.
	 * If one of the nodes is 0, meaning ground, then we
	 * only need to add it to its corresponding diagonal element. */
	if (n2 == 0) {
		this->G[n1-1][n1-1] += (1.0 / Z);
	}

	else if (n1 == 0) {
		this->G[n2-1][n2-1] += (1.0 / Z);
	}

	/* Otherwise we need to add it to multiple diagonals and two
	 * off-diagonal elements. */
	else {
		/* First add it to the diagonal elements */
		this->G[n1-1][n1-1] += (1.0 / Z);
		this->G[n2-1][n2-1] += (1.0 / Z);

		/* Now subtract it from the off-diagonal elements */
		this->G[n1-1][n2-1] -= (1.0 / Z);
		this->G[n2-1][n1-1] -= (1.0 / Z);
	}
}	

/* Calculate the node voltages */
int
Circuit::calcNodeVoltages ()
{
	int retVal;
	
	if (this->w == 0) {
		retVal = this->calcDCNodes ();
	}
	
	else {
		retVal = this->calcACNodes ();
	}
	
	return retVal;
}

/* Calculate node voltages for DC circuits */
int
Circuit::calcDCNodes ()
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
			gMatrix[i * nodes + j] = this->G[i][j].real ();
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

/* Calculate the node voltages for AC circuits */
int
Circuit::calcACNodes ()
{
	unsigned int i, j, k; /* Loop counters */
	unsigned int in;
	const unsigned int nodes = this->N; /* Number of nodes */
	const unsigned int sources = this->Ns; /* Number of sources */
	const unsigned int gSize = nodes * nodes;
	const unsigned int vSize = nodes;
	std::complex<double> gMatrix[gSize]; /* Copy of the conductance matrix */
	std::complex<double> vVector[vSize]; /* Copy of the source vector */
	lapack_int ipvt[nodes]; /* Needed by LAPACKE_dgesv to find the node voltages */
	lapack_int err; /* Returned value from LAPACKE_dgesv */	

	/* Fill gMatrix and vVector */
	for (i = 0; i < nodes; ++i) {
		vVector[i] = std::complex (0, 0);
		
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
	 * by solving the system of equations. The result is stored in vVector.
	 *
	 * gMatrix and vVector are complex, but C++ can't use C's complex numbers.
	 * However, they are interoperable, so we can recast std::complex to
	 * lapack_complex_double, which is an alias for the C complex double type. */
	err = LAPACKE_zgesv (LAPACK_ROW_MAJOR, nodes, 1, reinterpret_cast<lapack_complex_double*>(gMatrix), nodes, ipvt, reinterpret_cast<lapack_complex_double*>(vVector), 1);

	/* Check for errors */
	if (err == 0) {
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

	/* Return the status of LAPACKE_zgesv */
	return err;
}	

/* Print the voltages at all nodes */
void
Circuit::printNodeVoltages ()
{
	double mag;
	double angle;
	const double pi = 3.14159265;
	unsigned int i;

	for (i = 0; i < this->N; ++i) {
		/* If the real part of the voltage is negative we print it normally */
		if (this->vNode[i].real () < 0) {
			/* Check if the imaginary part is zero. If it is then we only need to print the real part. */
			if (this->vNode[i].imag () == 0) {
				std::cout << "V" << i + 1 << " = " << this->vNode[i].real () << "<0" << std::endl;
			}

			/* If the imaginary part is negative then we can just print both parts */
			else if (this->vNode[i].imag () < 0) {
				/* Calculate the magnitude of the voltage */
				mag = sqrt (pow (this->vNode[i].real (), 2) + pow (this->vNode[i].imag (), 2));

				/* Calculate the phase shift of the voltage */
				angle = pi + atan (this->vNode[i].imag () / this->vNode[i].real ());

				/* Convert the angle from radians to degrees */
				angle *= (180 / pi);

				std::cout << "V" << i + 1 << " = " << mag << "<" << angle << std::endl;
			}

			/* If the imaginary part is positive then we need to print a plus sign before printing the imaginary part.
			 * This is because we don't use format specifiers with std::cout. */
			else {
				/* Calculate the magnitude of the voltage */
				mag = sqrt (pow (this->vNode[i].real (), 2) + pow (this->vNode[i].imag (), 2));

				/* Calculate the phase shift of the voltage */
				angle = pi - atan (this->vNode[i].imag () / this->vNode[i].real ());

				/* Convert the angle from radians to degrees */
				angle *= (180 / pi);

				std::cout << "V" << i + 1 << " = " << mag << "<" << angle << std::endl;
			}
		}

		/* Otherwise the real part is positive so we need to add another space after the equal sign */
		else {
			/* Check if the imaginary part is zero. If it is then we only need to print the real part. */
			if (this->vNode[i].imag () == 0) {
				std::cout << "V" << i + 1 << " = " << this->vNode[i].real () << "<0" << std::endl;
			}

			/* If the imaginary part is negative then we can just print both parts */
			else if (this->vNode[i].imag () < 0) {
				/* Calculate the magnitude of the voltage */
				mag = sqrt (pow (this->vNode[i].real (), 2) + pow (this->vNode[i].imag (), 2));

				/* Calculate the phase shift of the voltage */
				angle = atan (this->vNode[i].imag () / this->vNode[i].real ());

				/* Convert the angle from radians to degrees */
				angle *= (180 / pi);

				std::cout << "V" << i + 1 << " = " << mag << "<" << angle << std::endl;
			}

			/* If the imaginary part is positive then we need to print a plus sign before printing the imaginary part.
			 * This is because we don't use format specifiers with std::cout. */
			else {
				/* Calculate the magnitude of the voltage */
				mag = sqrt (pow (this->vNode[i].real (), 2) + pow (this->vNode[i].imag (), 2));

				/* Calculate the phase shift of the voltage */
				angle = atan (this->vNode[i].imag () / this->vNode[i].real ());

				/* Convert the angle from radians to degrees */
				angle *= (180 / pi);

				std::cout << "V" << i + 1 << " = " << mag << "<" << angle << std::endl;
			}
		}
	}
}
