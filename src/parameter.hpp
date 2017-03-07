#ifndef __PARAMETER_HPP
#define __PARAMETER_HPP

class Parameter {
public:
	/// Constructor. Initialises the parameters with default values, that
	/// represent an demonstration simulation.
	Parameter();

	/// Loads the file with the given filename and tries to extract the
	/// definitions of parameters from it.
	/// 
	/// The file should have one parameter in each line and should have no
	/// lines containing anything but a parameter definition.
	/// Each line should be formated like this:
	/// paramname = value
	/// e.g.:
	/// rho0 = 1000.0
	///
	/// The parameter names are, for now, case sensitive and need to have the
	/// exact same name as what is checked in the load method.
	/// @todo Do this more user friendly and note somewhere which parameters
	/// are available. Also allow comments or empty lines and whatnot
	///
	/// @param filename char* The name of the file to load
	void Load(const char* filename);

	int N;
	double h;
	int R;
	double tend;
	double k;
	double g;
	double dt;
	double rho0;
	double mass;
	double mu;
	double dampening;
	double FSPressure;
	double FSGravity;
	double FSViscosity;
};
#endif // __PARAMETER_HPP