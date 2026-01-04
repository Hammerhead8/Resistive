/* Function definitions for the reimplementation of Resistive.
 *
 * Node voltages are found using nodal analysis with a conductance matrix */

#include "resistive.h"
#include <iostream>
#include <fstream>
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
	this->w.resize (1);
	this->w[0] = freq;

	/* Create the vectors Vin, inNode, and vNode */
	this->Vin.resize (1);
	this->inNode.resize (1);
	this->vNode.resize (1);

	/* Create the conductance matrix */
	this->G.resize (1);
	this->G[0].resize (1);
	this->G[0][0].resize (1);
}

/* Overloaded class constructor for performing AC analysis.
 * A conductance matrix is created for each frequency being analyzed.
 *
 * lowerExp is the exponent of the starting frequency (w0 = 10^lowerExp)
 * upperExp is the exponent of the final frequency (w1 = 10^upperFreq)
 * numPoints is the number of points per decade in the sweep */
Circuit::Circuit (int lowerExp, int upperExp, int numPoints)
{
	int numDecades;
	int orderOfMagnitude;
	int i, j, k;
	int freqSize;
	double stepSize;
	
	this->N = 0;
	this->Ns = 0;
	
	/* Calculate the required size of the frequency vector */
	numDecades = upperExp - lowerExp;
	this->w.resize (numDecades * numPoints + 1);
	this->vNode.resize (numDecades * numPoints + 1);
	
	/* Fill the vector with the indivual frequencies */
	k = 0;
	orderOfMagnitude = lowerExp;
	
	for (i = 0; i < numDecades; ++i) {
		stepSize = (pow (10, orderOfMagnitude + 1) - pow (10, orderOfMagnitude)) / numPoints;
		
		for (j = 0; j < numPoints; ++j) {
			this->w[k] = pow (10, orderOfMagnitude) + j * stepSize;
			k++;
		}
		
		orderOfMagnitude++;
	}
	
	this->w[k] = pow (10, orderOfMagnitude);
	
	this->Vin.resize (1);
	this->inNode.resize (1);
	
	freqSize = (upperExp - lowerExp) * numPoints + 1;
	
	/* Create the conductance matrix */
	this->G.resize (freqSize);
	for (i = 0; i < freqSize; ++i) {
		this->G[i].resize (1);
	}
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
	unsigned int i, j;
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

	if (this->N > this->G[0].size ()) {
		for (i = 0; i < this->w.size (); ++i) {
			this->G[i].resize (this->N);

			for (j = 0; j < this->G[i].size (); ++j) {
				this->G[i][j].resize (this->N);
			}
		}
	}

	/* Add the new resistor to the conductance matrix.
	 * If one of the nodes is 0, meaning ground, then we
	 * only need to add it to its corresponding diagonal element. */
	for (i = 0; i < this->w.size (); ++i) {
		if (n2 == 0) {
			this->G[i][n1-1][n1-1] += (1 / value);
		}

		else if (n1 == 0) {
			this->G[i][n2-1][n2-1] += (1 / value);
		}

		/* Otherwise we need to add it to multiple diagonals and two
		 * off-diagonal elements. */
		else {
			/* First add it to the diagonal elements */
			this->G[i][n1-1][n1-1] += (1 / value);
			this->G[i][n2-1][n2-1] += (1 / value);

			/* Now subtract it from the off-diagonal elements */
			this->G[i][n1-1][n2-1] -= (1 / value);
			this->G[i][n2-1][n1-1] -= (1 / value);
		}
	}
}

/* Add an inductor between nodes n1 and n2 */
void
Circuit::addInductor (unsigned int n1, unsigned int n2, double value)
{
	unsigned int i, j;
	
	/* Check if the circuit frequency is 0 (DC).
	 * If so we can't add the inductor */
	if (this->w[0] == 0) {
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
	for (i = 0; i < this->w.size (); ++i) {
		if ((n1 > this->N) || (n2 > this->N)) {
			this->N = max (n1, n2);
		}

		if (this->N > this->G.size ()) {
			this->G[i].resize (this->N);

			for (j = 0; j < this->G[i].size (); ++j) {
				this->G[i][j].resize (this->N);
			}
		}
		
		/* Calculate the impedance of the inductor */
		std::complex Z (0.0, this->w[i] * value);

		/* Add the new resistor to the conductance matrix.
		 * If one of the nodes is 0, meaning ground, then we
		 * only need to add it to its corresponding diagonal element. */
		if (n2 == 0) {
			this->G[i][n1-1][n1-1] += (1.0 / Z);
		}

		else if (n1 == 0) {
			this->G[i][n2-1][n2-1] += (1.0 / Z);
		}

		/* Otherwise we need to add it to multiple diagonals and two
		 * off-diagonal elements. */
		else {
			/* First add it to the diagonal elements */
			this->G[i][n1-1][n1-1] += (1.0 / Z);
			this->G[i][n2-1][n2-1] += (1.0 / Z);

			/* Now subtract it from the off-diagonal elements */
			this->G[i][n1-1][n2-1] -= (1.0 / Z);
			this->G[i][n2-1][n1-1] -= (1.0 / Z);
		}
	}
}

/* Add a capacitor between nodes n1 and n2 */
void
Circuit::addCapacitor (unsigned int n1, unsigned int n2, double value)
{
	unsigned int i, j;
	
	/* Check if the circuit frequency is 0 (DC).
	 * If it is then we can't add the capacitor so throw an exception. */
	if (this->w[0] == 0) {
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
	for (i = 0; i < this->w.size (); ++i) {
		if ((n1 > this->N) || (n2 > this->N)) {
			this->N = max (n1, n2);
		}

		if (this->N > this->G[i].size ()) {
			this->G.resize (this->N);

			for (j = 0; j < this->G.size (); ++j) {
				this->G[i][j].resize (this->N);
			}
		}
		
		/* Calculate the impedance of the capacitor */
		std::complex Z (0.0, -1 / (this->w[i] * value));

		/* Add the new resistor to the conductance matrix.
		 * If one of the nodes is 0, meaning ground, then we
		 * only need to add it to its corresponding diagonal element. */
		if (n2 == 0) {
			this->G[i][n1-1][n1-1] += (1.0 / Z);
		}

		else if (n1 == 0) {
			this->G[i][n2-1][n2-1] += (1.0 / Z);
		}

		/* Otherwise we need to add it to multiple diagonals and two
		 * off-diagonal elements. */
		else {
			/* First add it to the diagonal elements */
			this->G[i][n1-1][n1-1] += (1.0 / Z);
			this->G[i][n2-1][n2-1] += (1.0 / Z);

			/* Now subtract it from the off-diagonal elements */
			this->G[i][n1-1][n2-1] -= (1.0 / Z);
			this->G[i][n2-1][n1-1] -= (1.0 / Z);
		}
	}
}	

/* Calculate the node voltages.
 *
 * TODO:  Add a default argument giving the node whose voltage should be output.
 * If a value less than 1 or no value is given then return all values.
 * If an invalid node was given then print a message and return all voltages. */
int
Circuit::calcNodeVoltages ()
{
	int retVal;

	/* Before calculating check if there is at least
	 * voltage source connected to the circuit.
	 * If not then throw an exception. */
	if (this->Ns == 0) {
		throw (std::runtime_error (std::string ("There must be at least one voltage source connected to the circuit.")));
	}
	
	if (this->w[0] == 0) {
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
			gMatrix[i * nodes + j] = this->G[0][i][j].real ();
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
	this->vNode[0].resize (nodes);

	/* Now that gMatrix and vVector are filled we can calculate the nodes voltages
	 * by solving the system of equations. The result is stored in vVector. */
	err = LAPACKE_dgesv (LAPACK_ROW_MAJOR, nodes, 1, gMatrix, nodes, ipvt, vVector, 1);

	/* Check for errors */
	if (err == 0) {
		/* If no errors occured then copy the results from vVector to vNode */
		for (i = 0; i < nodes; ++i) {
			this->vNode[0][i] = vVector[i];
		}
	}

	else if (err < 0) {
		std::cerr << "Argument " << -1 * err << " is invalid." << std::endl;
		for (i = 0; i < nodes; ++i) {
			this->vNode[0][i] = NAN;
		}
	}

	else {
		std::cerr << "Conductance matrix is singular." << std::endl;
		for (i = 0; i < nodes; ++i) {
			this->vNode[0][i] = NAN;
		}
	}

	/* The last step is to delete gMatrix and vVector */
	delete[] gMatrix;
	delete[] vVector;

	/* Return the status of LAPACKE_dgesv */
	return err;
}

/* Calculate the node voltages for AC circuits
 *
 * TODO:  Add a default argument giving the node whose voltage should be output.
 * If a value less than 1 or no value is given then return all values.
 * If an invalid node was given then print a message and return all voltages. */
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
	
	/* Calculate the node voltages of the circuit
	 * for each frequency in the AC analysis */
	for (i = 0; i < this->w.size (); ++i) {
		/* Allocate vNode */
		this->vNode[i].resize (nodes);
		
		/* Fill gMatrix and vVector */
		for (j = 0; j < nodes; ++j) {
			vVector[j] = std::complex (0.0, 0.0);

			for (k = 0; k < nodes; ++k) {
				gMatrix[j * nodes + k] = this->G[i][j][k];
			}
		}

		/* Set the non-zero values in vVector from vIn for each conductance matrix */
		for (j = 0; j < sources; ++j) {
			in = this->inNode[j] - 1;
			vVector[in] = this->Vin[j];

			/* Set the diagonal element to 1 and the others to zero
			 * in the conductance matrix for the row corresponding to
			 * the source */
			for (k = 0; k < nodes; ++k) {
				if (k == in) {
					gMatrix[in * nodes + k] = 1;
				}

				else {
					gMatrix[in * nodes + k] = 0;
				}
			}
		}
		
		/* Now that gMatrix and vVector are filled we can calculate the nodes voltages
		 * by solving the system of equations. The result is stored in vVector.
		 *
		 * gMatrix and vVector are complex, but C++ can't use C's complex numbers.
		 * However, they are interoperable, so we can recast std::complex to
		 * lapack_complex_double, which is an alias for the C complex double type. */
		err = LAPACKE_zgesv (LAPACK_ROW_MAJOR, nodes, 1, reinterpret_cast<lapack_complex_double*>(gMatrix), nodes, ipvt, reinterpret_cast<lapack_complex_double*>(vVector), 1);

		/* Check for errors */
		if (err == 0) {
			for (j = 0; j < nodes; ++j) {
				this->vNode[i][j] = vVector[j];
			}
		}

		else if (err < 0) {
			std::cerr << "Argument " << -1 * err << " is invalid." << std::endl;
			for (j = 0; j < nodes; ++j) {
				this->vNode[i][j] = NAN;
			}
			
			return err;
		}

		else {
			std::cerr << "Conductance matrix is singular." << std::endl;
			for (j = 0; j < nodes; ++j) {
				this->vNode[i][j] = NAN;
			}
			
			return err;
		}
	}

	/* Return the status of LAPACKE_zgesv */
	return err;
}	

/* Print the voltages at all nodes */
void
Circuit::printNodeVoltages (int node)
{
	double mag;
	double angle;
	unsigned int i, j;
	
	/* If we are only using one input frequency */
	if (this->w.size () == 1) {
		/* If we are using a DC circuit, then we don't
		 * need to worry about the phase angle */
		if (this->w[0] == 0) {
			/* If the user wants only the voltage at a specific
			 * node then only print the voltage at that node */
			if ((node > 0) && (node <= this->N)) {
				std::cout << "V" << node << " = " << this->vNode[0][node - 1].real () << std::endl;
			}
			
			/* Otherwise print all of the node voltages */
			else {
				for (i = 0; i < this->N; ++i) {
					std::cout << "V" << i + 1 << " = " << this->vNode[0][i].real () << std::endl;
				}
			}
		}
		
		/* Otherwise we are using an AC circuit so we
		 * do need to worry about the phase angle */
		else {
			/* If the user wants only the voltage at a specific
			 * node then only print the voltage at that node */
			if ((node > 0) && (node <= this->N)) {
				/* Calculate the magnitude of the voltage */
				mag = std::abs (this->vNode[0][node - 1].real ());
				
				/* Calculate the phase shift of the voltage */
				angle = std::arg (this->vNode[0][node - 1]);

				/* Convert the angle from radians to degrees */
				angle *= (180 / M_PI);
				
				std::cout << "V" << node << " = " << mag << "<" << angle << std::endl;
			}
			
			/* Otherwise print all of the node voltages */
			else {
				for (i = 0; i < this->N; ++i) {
					/* Calculate the magnitude of the voltage */
					mag = std::abs (this->vNode[0][i]);
					
					/* Calculate the phase shift of the voltage */
					angle = std::arg (this->vNode[0][i]);

					/* Convert the angle from radians to degrees */
					angle *= (180 / M_PI);
					
					std::cout << "V" << i + 1 << " = " << mag << "<" << angle << std::endl;
				}
			}
		}
		
		/* Exit the function */
		return;
	}
	
	/* TODO:  Write voltage as a function of frequency to a file.
	 * They can all be in the same file for now but there
	 * needs to be some kind of seperator between nodes
	 * in the file that makes it clear which node the values
	 * correspond to.
	 *
	 * An alternative would be to create a different file for each node. */
	 
	 /* Since we have multiple frequencies the most useful format in which
	  * we can output the node voltages is in a file, with each node seperated
	  * from the others */
	 std::ofstream fp;
	 
	 fp.open ("file.txt");
	 
	 /* Write the points to the file */
	 for (i = 0; i < this->N; ++i) {
	 	fp << "Node " << i + 1 << std::endl;
	 	
	 	for (j = 0; j < this->w.size (); ++j) {
	 		/* Calculate the magnitude of the voltage */
			mag = std::abs (this->vNode[j][i]);
	 		
	 		/* Now convert the magnitude to dB */
	 		mag = 20 * std::log10 (mag);
	 		
	 		fp << this->w[j] << "\t" << mag << std::endl;
	 	}
	 	
	 	fp << "\n\n";
	 }
	 
	 fp.close ();
	 return;
}
