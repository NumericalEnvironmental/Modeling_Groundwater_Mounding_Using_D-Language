// model for groundwater perturbation from recharge basin(s) and pumping well(s)
// based on superposition of solutions by Hantush (1967) and Theis (1935)
//
// approximations include:
//
// (1) inherent approximations in analytical solutions (e.g., saturated aquifer thickness does not change appreciably compared to ambient)
// (2) numerical integration (trapezoidal rule, distributed along geometric progression of points)
// (3) exponential integral used in Theis expression is handled by expansion of the infinite series approximation out to 20 terms


import std.stdio; 			// I/O and file system
import std.math; 			// math function library
import std.mathspecial; 	// error function
import std.conv; 			// automatic conversions between types
import std.typecons; 		// misc. meta-programming support, including type constructors, tuples, etc.
import std.string; 			// string function support


// *** classes and data structures ***


struct Model {
	// model parameters and switches
	double t_end, x0, xf, y0, yf, dt, d_max;
	int x_int, y_int;
	bool gridder, particle_trace, monitor;
	this(double t_end, bool gridder, double x0, double xf, int x_int, double y0, double yf, int y_int, bool particle_trace, double dt, double d_max, bool monitor){ 	// constructor
        this.t_end = t_end; 							// end-time for simulation
		this.gridder = gridder;							// switch to calculate final state on grid (boolean)
		this.x0 = x0; 									// grid minima, maxima, and discretization, both axes
		this.xf = xf;
		this.x_int = x_int;
		this.y0 = y0;
		this.yf = yf;
		this.y_int = y_int;
        this.particle_trace = particle_trace; 			// switch to employ particle tracking (boolean)
		this.dt = dt; 									// time-step size for particle tracking and time series
		this.d_max = d_max; 							// maximum travel distance for particle over a given time step
		this.monitor = monitor; 						// switch to calculate time series at select locations (boolean)
		}
	}


struct Aquifer {
	// basic aquifer properties
	double K, Sy, h0, dhx, dhy, phi;
	this(double K, double Sy, double h0, double dhx, double dhy, double phi){ 				// constructor
        this.K = K; 									// hydraulic conductivity
		this.Sy = Sy; 									// specific yield
		this.h0 = h0; 									// reference saturated thickness
		this.dhx = dhx; 								// ambient hydraulic gradient component parallel to x-axis
		this.dhy = dhy; 								// ambient hydraulic gradient component parallel to y-axis
		this.phi = phi; 								// porosity (only used for particle tracking)
		}	
	}

	
struct History {
	// container for operating history (time, quantity); order is important, so can't used associative array
	double[] time;
	double[] flux;
	this(double[] time, double[] flux){ 				// constructor
        this.time = time; 								// hydraulic conductivity
        this.flux = flux; 								// aquifer thickness
		}		
	}	
	
	
class Monitor {

	// properties and methods w.r.t. time series at selected points

	string[] name; 					// names of monitoring points/wells
	double[] x; 					// coordinates
	double[] y;
	
	this(){ 		// constructor
		// read monitor point locations
		string line_input;
		string[] parsed_input;	
		int line_index;
		auto input_file = File("monitor_locs.txt","r");
		while (!input_file.eof()) {
			line_input = input_file.readln();
			if (line_index > 0) { 	 					// skip header
				parsed_input = split(line_input);
				this.name ~= parsed_input[0];
				this.x ~= to!double(parsed_input[1]);
				this.y ~= to!double(parsed_input[2]);
				}
			++line_index;
			}
		input_file.close();	
		}		
	
	public void WaterLevel(Model model, Aquifer aquifer, Mound[] basin, Well[] well){

		// calculate potentiometric surface at each monitor point for each time step and write to file

		string line_string;
		double h;
		double t = 0.;

		// write output file header
		auto well_file = File("water_levels.txt", "w");
		line_string = "Time";
		for (int i = 0; i < name.length; ++i){line_string ~= "\t" ~ name[i];}
		well_file.writeln(line_string);		

		// compute superposed potentials and append to output file
		while (t <= model.t_end){
			t += model.dt;
			line_string = to!string(t);
			for (int i = 0; i < x.length; ++i){
				h = h_tot(x[i], y[i], t, aquifer, basin, well);
				line_string ~= "\t" ~ to!string(h);
				}
			well_file.writeln(line_string);
			}
	
		//close output file
		well_file.close(); 
		}
		
	}	
	
	
class Tracker {

	// properties and methods w.r.t. particle tracking
	
	double[] x; 			// particle positions; output to file at each step, so time series arrays are not created
	double[] y;

	this(){ 				// constructor
		// read initial particle positions
		string line_input;
		string[] parsed_input;	
		int line_index;
		auto input_file = File("particles_init.txt","r");
		while (!input_file.eof()) {
			line_input = input_file.readln();
			if (line_index > 0) { 	 					// skip header
				parsed_input = split(line_input);
				this.x ~= to!double(parsed_input[0]);
				this.y ~= to!double(parsed_input[1]);
				}
			++line_index;
			}
		input_file.close();	
		}	
	
	private Tuple!(double, double) Velocity(double x, double y, double t, Aquifer aquifer, Mound[] basin, Well[] well){
		// compute local pore velocity components by numerical differentiation of superimposed potential functions
		double vx, vy;
		double delta = 0.001;
		vx = -aquifer.K * (h_tot(x+0.5*delta, y, t, aquifer, basin, well)
			- h_tot(x-0.5*delta, y, t, aquifer, basin, well)) / (delta * aquifer.phi);
		vy = -aquifer.K * (h_tot(x, y+0.5*delta, t, aquifer, basin, well)
			- h_tot(x, y-0.5*delta, t, aquifer, basin, well)) / (delta * aquifer.phi);	
		return tuple(vx, vy);	
		}
	
	private Tuple!(double, double) FixPosition(double x0, double y0, double t, double vx0, double vy0, double dt, Aquifer aquifer, Mound[] basin, Well[] well){
		// iteratively determine particle position by accurate assessment of average velocity
		double xn, yn, vxn, vyn;
		double vx, vy;
		double vm0 = 0.;
		double vm = 0.;
		double conv_threshold = 0.001;
		int iter = 0;
		int max_iter = 2;
		t += 0.5*dt; 														// calculate average velocity --> middle of time step
		vx = vx0;
		vy = vy0;
		vm0 = sqrt(vx^^2 + vy^^2); 											// magnitude of current velocity estimate		
		while ((1.0 - abs(vm/vm0) > conv_threshold) && (iter < max_iter)){
			xn = x0 + vx*dt; 												// implied particle position estimate (vx, vy)
			yn = y0 + vy*dt;
			auto vn = Velocity(xn, yn, t, aquifer, basin, well); 			// implied local velocity at new position
			vxn = vn[0];
			vyn = vn[1];		
			vx = (vx0 + vxn)/2; 											// update velocity estimate: average of initial and final positions
			vy = (vy0 + vyn)/2;		
			vm = sqrt(vx^^2 + vy^^2);
			++iter;
			}
		xn = x0 + vx*dt; 													// final position estimate
		yn = y0 + vy*dt;			
		return tuple(xn, yn);
		}
	
	public void Track(Model model, Aquifer aquifer, Mound[] basin, Well[] well){
	
		// track movement of each particle; append to output file
	
		string line_string;
		double t, x0, y0, vx0, vy0, offset, dt_eff, r;
		int iter;
		int max_iter = 5;
		bool converge, live; 														// "live" --> particle is viable (i.e., not near well)
	
		// write output file header
		auto track_file = File("tracks_out.txt", "w");
		line_string = "x" ~ "\t" ~ "y" ~ "\t" ~ "time"~ "\t" ~ "particle"; 			// write tracking file column header
		track_file.writeln(line_string);		

		// for each particle ...
		for (int i = 0; i < x.length; ++i){
			writeln(to!string(i));
			t = 0.;		
			line_string = to!string(x[i]) ~ "\t" ~ to!string(y[i]) ~ "\t" ~ to!string(t)~ "\t" ~ to!string(i);
			track_file.writeln(line_string);										// post initial particle location to output file	
			live = true;
			while ((t < model.t_end) && (live == true)){
				auto v0 = Velocity(x[i], y[i], t, aquifer, basin, well); 			// use starting velocity for iteration
				vx0 = v0[0];
				vy0 = v0[1];			
				// iterate to find best estimate of particle position
				dt_eff = model.dt;
				x0 = x[i];
				y0 = y[i];	
				iter = 0;
				converge = false;
				while (converge==false){
					auto p = FixPosition(x[i], y[i], t, vx0, vy0, dt_eff, aquifer, basin, well);		
					x[i] = p[0];
					y[i] = p[1];
					assert(isnan(x[i])==false, "Particle tracking failure; check convergence.");
					offset = Dist(x[i], y[i], x0, y0);
					if ((offset <= model.d_max) || (iter >= max_iter)){converge=true;}
					else {
						dt_eff *= 0.5; 						// if no convergence, reduce time step
						x[i] = x0; 							// reset particle starting position
						y[i] = y0;						
						++iter;
						}
					}
					
				// append to output file
				t += dt_eff;
				line_string = to!string(x[i]) ~ "\t" ~ to!string(y[i]) ~ "\t" ~ to!string(t)~ "\t" ~ to!string(i);
				track_file.writeln(line_string);				
			
				// deactivate particle if it is within d_max from an extraction well
				foreach(pump; well){
					r = Dist(x[i], y[i], pump.x, pump.y);
					if (r <= model.d_max) {
						live = false; 					// finished tracking this particle, if it's close to a well
						break;
						} 			
					}
			
				} 	// end of time stepping, per particle

			} 	// end of for-loop, for particles
			
		//close output file
		track_file.close(); 			
		
		} 		// end of method definition
		
	} 		// end of class definition
	

class Grid {
	
	// x-y grid for calculating head distribution
	
	double[] x;
	double[] y;
	double[] h;
	double dx, dy;
	
	this(Model model){ 		// constructor
		// initialize grid points
		dx = (model.xf - model.x0) / model.x_int;
		dy = (model.yf - model.y0) / model.y_int;		
		for (int j = 0; j < model.y_int; ++j){
			for (int i = 0; i < model.x_int; ++i){
				this.x ~= model.x0 + (i + 0.5) * dx;
				this.y ~= model.y0 + (j + 0.5) * dy;	
				this.h ~= 0.;				
				}
			}
		}

	public void Surface(Model model, Aquifer aquifer, Mound[] basin, Well[] well){
		// compute groundwater potentiometric surface (as defined by grid of points) from superposition
		for (int i = 0; i < x.length; ++i){
			h[i] = h_tot(x[i], y[i], model.t_end, aquifer, basin, well);
			}
		}	
		
	public void WriteOutput(){
		// hydraulic head output at the end of the simulation
		string line_string;
		auto output_file = File("head_out.txt", "w");
		line_string = "x" ~ "\t" ~ "y" ~ "\t" ~ "head"; 			// write header
		output_file.writeln(line_string);
		for (int i = 0; i < x.length; ++i){
			line_string = to!string(x[i]) ~ "\t" ~ to!string(y[i]) ~ "\t" ~ to!string(h[i]);
			output_file.writeln(line_string);
			}			
		output_file.close();		
		}
		
	}
	
	
class Well {

	// properties and method w.r.t. pumping well (modeled with approximation to Theis solution)
	
	double x, y;
	Aquifer aquifer;	
	History well_history;
	
	this(double x, double y, History well_history, Aquifer aquifer){ 		// constructor
		this.x = x; 									// well location (x)
		this.y = y; 									// well location (y)
		this.well_history = well_history; 				// pumping history
		this.aquifer = aquifer; 						// basic aquifer properties
		}	

	public double dh_time_int(double r, double t, double b){
		// net change in head, per integration of all pumping periods
		int i = 0; 					// array index counter
		double dh = 0;
		double dt;
		while (i < well_history.time.length){
			if (t > well_history.time[i]) {
				dt = t - well_history.time[i];
				dh += s(r, dt, b, well_history.flux[i]);
				if (i > 0){dh -= s(r, t - dt, b, well_history.flux[i-1]);}
				}
			else {break;}			// remainder of well history is beyond t, so quit ...
			++i;
			}
		return dh;
		}
		
	private double s(double r, double t, double b, double Q){
		// Theis model
		double u;
		u = (r^^2 * aquifer.Sy) / (4 * aquifer.K * b * t);
		return Q / (4 * PI * aquifer.K * b) * W(u);
		}
	
	private double W(double u){
		// exponential integral function approximation
		double w;
		int n = 20; 				// number of terms to carry for series approximation
		double sum = 0.;
		if (u < 6){
		// summation approximation
			for (int i = 1; i <= n; ++i){sum += Odd(i) * u^^i / (i * gamma(to!double(i)+1));}
			w = -0.5772 - log(u) + sum;
			}
		else {w = 0.;}
		return w;
		}
	
	}

	
class Mound {

	// properties and method w.r.t. rectangular aquifer recharge zone

	double xm, ym, l, a, theta;
	Aquifer aquifer;
	History basin_history;
	
	this(double xm, double ym, double l, double a, double angle_deg, History basin_history, Aquifer aquifer){ 			// constructor
		this.xm = xm; 									// recharge basin/mound centroid location
		this.ym = ym;
		this.l = l; 									// half-length of recharge basin
		this.a = a; 									// half-width of recharge basin	
		this.theta = angle_deg * PI / 180;				// rotation angle from due east; convert degrees to radians
		this.basin_history = basin_history; 			// recharge history
		this.aquifer = aquifer; 						// basic aquifer properties
		}
	
	public double dh_time_int(double x, double y, double t, double b){
	
		// net change in head, per integration of all recharge periods
		
		int i = 0; 					// array index counter
		double dh = 0;
		double dt, x_mod, y_mod;
		
		// rotate coordinates to match basin orientation
		auto rotation = Rotate(x - xm, y - ym, -theta);
		x_mod = rotation[0]; 			
		y_mod = rotation[1];		
		
		// compute dh from h^2 - h0^2
		while (i < basin_history.time.length){
			if (t > basin_history.time[i]) {
				dt = t - basin_history.time[i];
				dh += sqrt(dh2(x_mod, y_mod, dt, b, basin_history.flux[i]) + aquifer.h0^^2) - aquifer.h0;
				if (i > 0){dh -= sqrt(dh2(x_mod, y_mod, dt, b, basin_history.flux[i-1]) + aquifer.h0^^2) - aquifer.h0;}
				}
			else {break;}			// remainder of well history is beyond t, so quit ...
			++i;
			}
		return dh;
		}
	
	private Tuple!(double, double) Rotate(double x, double y, double theta){
		// coordinate system rotation; theta is expressed in radians
		double x_prime, y_prime;
		x_prime = x*cos(theta) - y*sin(theta);
		y_prime = x*sin(theta) + y*cos(theta);	
		return tuple(x_prime, y_prime);
		}

	private double dh2(double x, double y, double t, double b, double w) {
	
		// compute mounding as a result of recharge as a function of x, y, and t
		
		double[2] alpha;
		double[2] beta;	
		
		auto S_terms = AlphaBeta(x, y, t, b);
		alpha = S_terms[0]; 			
		beta = S_terms[1];				
	
		// function returns h^2 - h0^2
		return w/(2.*aquifer.K) * (aquifer.K * b / aquifer.Sy) * t
			* (S(alpha[0], beta[0]) + S(alpha[0], beta[1])
			+ S(alpha[1], beta[0]) + S(alpha[1], beta[1]));
		}
	
	private Tuple!(double[2], double[2]) AlphaBeta(double x, double y, double t, double b){
		// return coefficients used in integral in Hantush (1967) solution
		double[2] alpha;
		double[2] beta;
		alpha[0] = (l + x)/sqrt(4. * aquifer.K * b / aquifer.Sy * t);
		alpha[1] = (l - x)/sqrt(4. * aquifer.K * b / aquifer.Sy * t);		
		beta[0] = (a + y)/sqrt(4. * aquifer.K * b / aquifer.Sy * t);
		beta[1] = (a - y)/sqrt(4. * aquifer.K * b / aquifer.Sy * t);			
		return tuple(alpha, beta);
		}

	private double f(double tau, double alpha, double beta){
		// integral term included in Hantush (1967) groundwater mounding analytical expression	
		return erf(alpha/sqrt(tau)) * erf(beta/sqrt(tau));
		}	

	private double S(double alpha, double beta){
		// trapezoidal rule integration of Hantush function term
		// note: use logarithmically-weighted integration points (to approach singularity at 0)	
		int n_int = 100; 							// number of integration intervals
		double pt_0 = 0.0001; 						// limits of integration are always 0 to 1
		double pt_1 = 1.;
		double dx_log, a, b, z0, z1, fa, fb;
		double sum = 0.;
		z0 = log10(pt_0);
		z1 = log10(pt_1);
		dx_log = (z1 - z0)/n_int;
		a = 10^^z0;
		fa = f(a, alpha, beta);
		for (int i = 1; i <= n_int; ++i){
			b = 10^^(z0 + i*dx_log); 
			fb = f(b, alpha, beta);
			sum += (b-a) * (fa + fb)/2;
			a = b;
			fa = fb;
			}
		return sum;
		}
	
	}


// *** utility functions ***


int Odd(int i){
	// return +1 if i is odd, and -1 if it is even
	int odd;
	if ((i & 1) == 0){odd = -1;}
	else {odd = 1;}
	return odd;
	}


double Dist(double x1, double y1, double x2, double y2){
	// Euclidian distance
	return sqrt( (x1-x2)^^2 + (y1-y2)^^2 );
	}

	
Aquifer ReadAquifer(){
	// read aquifer properties file and return structure
	double K, Sy, h0, dhx, dhy, phi;
	string line_input;
	string[] parsed_input;
	auto input_file = File("aquifer.txt","r");
	while (!input_file.eof()) {
		line_input = input_file.readln();	
		parsed_input = split(line_input);
		switch (parsed_input[0]) {
			case "Sy": 									// specific yield
				Sy = to!double(parsed_input[1]);
				break;		
			case "h0": 									// reference aquifer saturated thickness
				h0 = to!double(parsed_input[1]);
				break;	
			case "dhx": 								// ambient hydraulic gradient component parallel to x-axis
				dhx = to!double(parsed_input[1]);
				break;	
			case "dhy": 								// ambient hydraulic gradient component parallel to y-axis
				dhy = to!double(parsed_input[1]);
				break;	
			case "phi": 								// porosity
				phi = to!double(parsed_input[1]);
				break;					
			default: 									// hydraulic conductivity
				K = to!double(parsed_input[1]);				
				break;
			}		
		}
	input_file.close();		
	auto aquifer = Aquifer(K, Sy, h0, dhx, dhy, phi);
	return aquifer;
	}


Model ReadModelParams(){
	// read model parameters from file and populate data structure
	double t_end, x0, xf, y0, yf, dt, d_max;
	int x_int, y_int;
	bool gridder, particle_trace, monitor;
	string line_input;
	string[] parsed_input;
	auto input_file = File("model.txt","r");
	while (!input_file.eof()) {
		line_input = input_file.readln();	
		parsed_input = split(line_input);
		switch (parsed_input[0]) {
			case "t_end": 									// end-time of simulation
				t_end = to!double(parsed_input[1]);
				break;
			case "gridder": 								// whether or not to model groundwater elevation across grid (boolean)
				gridder = to!bool(parsed_input[1]);
				break;						
			case "x0": 										// grid origin, x-axis
				x0 = to!double(parsed_input[1]);
				break;		
			case "xf": 										// grid end, x-axis
				xf = to!double(parsed_input[1]);
				break;	
			case "x_int": 									// grid discretization density, x-axis
				x_int = to!int(parsed_input[1]);
				break;	
			case "y0": 										// grid origin, y-axis
				y0 = to!double(parsed_input[1]);
				break;		
			case "yf": 										// grid end, y-axis
				yf = to!double(parsed_input[1]);
				break;	
			case "y_int": 									// grid discretization density, y-axis
				y_int = to!int(parsed_input[1]);
				break;	
			case "particle_trace": 							// whether or not to employ particle tracking (boolean)
				particle_trace = to!bool(parsed_input[1]);
				break;		
			case "dt": 										// time-step size for particle tracking and time series
				dt = to!double(parsed_input[1]);
				break;	
			case "monitor": 								// whether or not to model time series (boolean)
				monitor = to!bool(parsed_input[1]);
				break;						
			default: 										// maximum particle displacement, per time step
				d_max = to!double(parsed_input[1]);				
				break;
			}
		}
	input_file.close();	
	Model model = Model(t_end, gridder, x0, xf, x_int, y0, yf, y_int, particle_trace, dt, d_max, monitor);
	return model;
	}


History ReadTimeSeries(string file_name){
	// read in a history file (time + flux columns) and return history object
	string line_input;
	string[] parsed_input;	
	int line_index;
	double[] time;
	double[] flux;
	auto input_file = File(file_name,"r");
	while (!input_file.eof()) {
		line_input = input_file.readln();
		if (line_index > 0) { 	 					// skip header
			parsed_input = split(line_input);
			time ~= to!double(parsed_input[0]);
			flux ~= to!double(parsed_input[1]);
			}
		++line_index;
		}
	input_file.close();	
	History time_series = History(time, flux);
	return time_series;
	}

	
Mound[] ReadMounds(Aquifer aquifer){
	// read recharge basin specs + recharge history; return array of mound objects
	string line_input;
	string[] parsed_input;	
	int line_index;
	string name;
	double xm, ym, l, a, angle_deg;
	History basin_history;
	Mound[] basin;
	auto input_file = File("basins.txt","r");
	while (!input_file.eof()) {
		line_input = input_file.readln();
		if (line_index > 0) { 	 														// skip header
			parsed_input = split(line_input);
			name = parsed_input[0]; 													// basin name (matched to input file holding recharge history)
			xm = to!double(parsed_input[1]); 											// center of basin
			ym = to!double(parsed_input[2]);			
			l = to!double(parsed_input[3])/2; 											// half-length (along x-axis) 	
			a = to!double(parsed_input[4])/2; 											// half-width (along y-axis) 				
			angle_deg = to!double(parsed_input[5]); 									// rotation angle of basin with respect to due east
			basin_history = ReadTimeSeries(name ~ "_basin.txt"); 						// read recharge history (from corresponding file)
			Mound basin_member = new Mound(xm, ym, l, a, angle_deg, basin_history, aquifer);
			basin ~= basin_member;
			}
		++line_index;
		}
	input_file.close();	
	return basin;
	}

	
Well[] ReadWells(Aquifer aquifer){
	// read well names and locations + pumping history; return array of well objects
	string line_input;
	string[] parsed_input;	
	int line_index;
	string name;
	double x, y;
	History well_history;
	Well[] well;
	auto input_file = File("wells.txt","r");
	while (!input_file.eof()) {
		line_input = input_file.readln();
		if (line_index > 0) { 	 														// skip header
			parsed_input = split(line_input);
			name = parsed_input[0]; 													// well name (matched to input file holding pumping history)
			x = to!double(parsed_input[1]); 											// well location
			y = to!double(parsed_input[2]);			
			well_history = ReadTimeSeries(name ~ "_well.txt"); 							// read pumping history (from corresponding file)
			Well well_member = new Well(x, y, well_history, aquifer);
			well ~= well_member;
			}
		++line_index;
		}
	input_file.close();	
	return well;
	}	
	
	
double h_tot(double x, double y, double t, Aquifer aquifer, Mound[] basin, Well[] well){
	// solve for total head at location (x, y) at time t_end
	double r;
	double dh = 0.;	
	double dh_m;
	double dh_w;	
	int num_m_iter = 2;
	int num_w_iter = 2;	
	foreach(recharge; basin){ 				// recharge basins
		dh_m = 0.;
		for (int i = 1; i <= num_m_iter; ++i){dh_m = recharge.dh_time_int(x, y, t, aquifer.h0 + 0.5*dh_m);}
		dh += dh_m;
		}
	foreach(pump; well){ 					// pumping wells
		dh_w = 0.;	
		r = Dist(x, y, pump.x, pump.y);
		for (int i = 1; i <= num_w_iter; ++i){dh_w = pump.dh_time_int(r, t, aquifer.h0 + 0.5*dh_w);}
		dh += dh_w;		
		}
	// add regional gradient (assumes uniform regional saturated thickness, slope of water table = slope of aquifer base)
	return aquifer.h0 + dh + x*aquifer.dhx + y*aquifer.dhy;
	}

	
// *** main function ***


void main(){

	Aquifer aquifer;
	Model model;
	Mound[] basin;
	Well[] well;
	
	// read various input files
	
	model = ReadModelParams(); 								// model parameters, returned from reader function (struct, not a class)
	writeln("Read model parameters.");
	
	aquifer = ReadAquifer(); 								// aquifer properties, returned from reader function (struct, not a class)
	writeln("Read aquifer properties.");

	basin = ReadMounds(aquifer); 							// recharge basin(s) properties (class), including recharge histories (struct)
	writeln("Read recharge basin(s) characteristics.");

	well = ReadWells(aquifer); 								// pumping well(s) properties (class), including pumping histories (struct)
	writeln("Read well(s) characteristics.");	
	

	// call potentiometric superposition model, as required for output(s)
	
	if (model.gridder==true){ 								// set up grid class; calculate heads and write output
		writeln("Processing grid.");	
		Grid grid = new Grid(model); 							
		grid.Surface(model, aquifer, basin, well);
		grid.WriteOutput();
		}
	
	if (model.monitor == true){ 							// conduct time series calculations at monitoring points
		writeln("Processing time series at monitor points.");
		Monitor monitor = new Monitor();
		monitor.WaterLevel(model, aquifer, basin, well);
		}
	
	if (model.particle_trace == true){ 						// conduct particle tracking, if applicable
		writeln("Tracking particles ...");
		Tracker trace = new Tracker();
		trace.Track(model, aquifer, basin, well);
		}

	writeln("Done.");
	
	}