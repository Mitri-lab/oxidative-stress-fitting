#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
//#include <chrono>
#include <ctime>
#include <math.h>
#include <cstdlib>
#include <stdlib.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <fstream>
#include <sstream>
#include <string>
#include <string.h>
#include <array>
#include <vector>
#include <algorithm>
#include <cfloat>

using namespace std;

bool is_inf_or_nan(double x) 
{
    return !(x <= DBL_MAX && x >= -DBL_MAX); 
}    


int equationsnoros (double t, const double y[], double f[], void *params)
{

    double* p = (double *) params;
      
    
    double rn1 = p[0];
    double rn2 = p[1];
    double Kn1 = p[2];
    double Kn2 = p[3];
    double Yn1 = p[4];
    double Yn2 = p[5];
    double rc1 = p[6];
    double rc2 = p[7];
    double Kc1 = p[8];
    double Kc2 = p[9];
    double Yc1 = p[10];
    double Yc2 = p[11];
    double beta1 = p[12];
    double beta2 = p[13];
    double gamma1 = p[14];
    double gamma2 = p[15];
	
	  // y[0] -> Ct, y[1] -> At, y[2] -> unknown MM nutrient, y[3]-> LA
	
	
	f[0] = rn1*y[2]/(y[2] + Kn1)*y[0] + rc1*y[3]/(y[3] + Kc1)*y[0] - (beta1 + gamma1*t)*y[0]*y[3];
	
	f[1] = rn2*y[2]/(y[2] + Kn2)*y[1] + rc2*y[3]/(y[3] + Kc2)*y[1] - (beta2 + gamma2*t)*y[1]*y[3];
	
	f[2] = - 1/Yn1*rn1*y[2]/(y[2] + Kn1)*y[0] - 1/Yn2*rn2*y[2]/(y[2] + Kn2)*y[1];
		
	f[3] =  -1/Yc1*rc1*y[3]/(y[3] + Kc1)*y[0] - 1/Yc2*rc2*y[3]/(y[3] + Kc2)*y[1];
	
	
	return GSL_SUCCESS;
	
	
}





int main(int argc, char *argv[])
{

    
    
    double rc1 = 2.378423;
    double Kc1 = 0.0005943006;
    double beta1 = 0.2994206;
    double Yc1 = 1382887081;
    double gamma1 = 3.014306e-06;
    double rn1 = 346.5335;
    double Kn1 = 0.100519;
    double Yn1 = 1277901049;
    
    double rc2 = 1.941358;
    double Kc2 = 0.001000012;
    double Yc2 =2380679472;
    double beta2 = 3.359035;
    double gamma2 = 0.5627613;
    double rn2 = 891.433;
    double Kn2 = 1.013203;
    double Yn2 = 313918816;

     vector<double> constparam;
	
	 constparam.push_back(rn1);
	 constparam.push_back(rn2);
	 constparam.push_back(Kn1);
	 constparam.push_back(Kn2);	 
	 constparam.push_back(Yn1);
	 constparam.push_back(Yn2);	 
	 constparam.push_back(rc1);
	 constparam.push_back(rc2);
	 constparam.push_back(Kc1);
	 constparam.push_back(Kc2);
	 constparam.push_back(Yc1);
	 constparam.push_back(Yc2);
	 constparam.push_back(beta1);
	 constparam.push_back(beta2);
	 constparam.push_back(gamma1);
	 constparam.push_back(gamma2);
	
	 

	 
	vector<double> vecparams = constparam;
	int tailleparams = vecparams.size();
	double params[tailleparams]; 
	for (int i = 0; i<tailleparams; ++i) {
		params[i] = vecparams[i];
	}


	
	double bact1init = 1e5; // Ct
    double bact2init = 1e5; //At
    double LAinit = 0.75;   // LA condition (0.01 or 0.75)
    double compini = 0.01;
    
    unsigned int nbequ = 4;
    
    
     const gsl_odeiv_step_type * T = gsl_odeiv_step_rkf45;
  
    gsl_odeiv_step * st = gsl_odeiv_step_alloc (T, nbequ); // creation an instance of the stepping function
	gsl_odeiv_control * control = gsl_odeiv_control_y_new (1e-8, 0.0); //
	gsl_odeiv_evolve * evol = gsl_odeiv_evolve_alloc (nbequ); // instance of an evolution function for a system of 1 dimensions
	gsl_odeiv_system sys = {equationsnoros, NULL, nbequ, &params}; // the Jacobian is useless with this method : use NULL instead
	
	double tbatch = 3.0;
	double tinteg = tbatch;
	int nbtransfer = 5; 
	double t = 0.0; // time span of integration
	double hfabs = 1e-10; // fabsolute accuracy;
	double y[nbequ]; 
	
	string typ;
	if (bact1init > 0 && bact2init > 0) {typ = "cocult";}
	if (bact1init > 0 && bact2init <= 0) {typ = "monoCt";}
	if (bact1init <= 0 && bact2init > 0) {typ = "monoAt";}
		
	
	string coextxt = "Model1_" + typ + "_" + to_string(nbtransfer) + ".txt";
		
		
	ofstream fout(coextxt.c_str(), ios::out);
	fout << "LA_initialconc"  <<"\t" << "dilu" << "\t"   << "nbtransfer" << "\t" << "finaltime" << "\t" << "Ct" << "\t" << "At" << "\t"  << "MMnu" << "\t"<< "LA" << "\t" << "finstate" <<   "\t" << "exttimeat" << endl;


	bool valnan = 0;
	
	
			
	int pasenregi = 100;
	

	
	double dilumin;
	double dilumax; 
	vector<double> vecdilu;	
	double pasdilu;	

	
	dilumin = 2;
	dilumax = 200;
	pasdilu= 2;
	double dilu = dilumin;
	
	while (dilu < dilumax + pasdilu) {
		vecdilu.push_back(dilu);
		dilu = dilu + pasdilu;	
	}
	

	

	double lamin;
	double lamax; 
	vector<double> vecla;	
	double pasla;	
	
	lamin = 0.01;
	lamax = 1;
	pasla= 0.01;
	double laini = lamin;
	
	while (laini < lamax + pasla) {
		vecla.push_back(laini);
		laini = laini + pasla;	
	}
	

cout << vecdilu.size() << endl;
cout << vecla.size() << endl;

	
	for (int pbif = 0; pbif != vecdilu.size(); ++pbif) { 
		for (int pbif2 = 0; pbif2 != vecla.size(); ++pbif2) { 
			
				
			LAinit = vecla[pbif2];
			dilu = vecdilu[pbif];
			
			cout << LAinit << endl;
			cout << dilu << endl;
			
			y[0] = bact1init;
			y[1] = bact2init;
			y[2] = compini;
			y[3] = LAinit;
			int tpoint = 0;		
			int transfer=0;
			t = 0.0;
			tinteg = tbatch;

			int succev = 0;
			int nbsuccev = 5;
			double eqthr = 1; 
			
			double Atfinal = 0.0;
			double Ctfinal = 0.0;	
			double Atprev = 0.0;
			double Ctprev = 0.0;
			
			
			double Atfinalfin = 0.0;
			double Ctfinalfin = 0.0;
			double MMnufinalfin = 0.0;
			double LAfinalfin = 0.0;
			
	
			bool notyetext = 1;
			int exttime = 0;
	
	
			while (transfer < nbtransfer and succev < nbsuccev) {
			
				
				Ctprev = y[0];
				Atprev = y[1];
				
				while (t < tinteg)		{
					
					int status = gsl_odeiv_evolve_apply (evol, control, st, &sys, &t, tinteg, &hfabs, y); // integration
					if (status != GSL_SUCCESS) {break;}
					
					for (int i = 0; i!=4; ++i) {
						if (y[i] < 1e-8) {y[i] = 0;}
						
				}
						
					
					if (is_inf_or_nan(y[0]) or is_inf_or_nan(y[1])  or is_inf_or_nan(y[2]) or is_inf_or_nan(y[3])  or is_inf_or_nan(y[4]) ) {
						valnan=1;
						cout << "NaN!" << endl;
						break;}        
							 
					tpoint++;
				} 
				
				Atfinalfin = y[1];
				Ctfinalfin = y[0];
				MMnufinalfin = y[2];
				LAfinalfin =y[3];
				
				y[0] = y[0] /dilu;
				y[1] = y[1] / dilu;
				y[2] = compini;
				y[3] = LAinit;
				
				Ctfinal = y[0];
				Atfinal = y[1];	
				
				if (Atfinal < eqthr && notyetext) {
					exttime = transfer; 
					notyetext = 0;}
					
				if ( abs(Ctfinal - Ctprev) < eqthr and abs(Atfinal - Atprev) < eqthr ) { succev++;}	

			
				
				transfer++;
				tinteg = tinteg + tbatch;
				
			}

				if (notyetext) {exttime = nbtransfer;} 

				string finstate;
				


			 if (Atfinalfin > eqthr and Ctfinalfin > eqthr) { finstate = "coex";}
			 if (Atfinalfin < eqthr and Ctfinalfin > eqthr) { finstate = "Ct";}
			 if (Atfinalfin > eqthr and Ctfinalfin < eqthr) { finstate = "At";}
			 if (Atfinalfin < eqthr and Ctfinalfin < eqthr) { finstate = "ext";}
				
			
				cout << t << endl;
				fout << LAinit  <<"\t" << dilu  << "\t" << transfer << "\t" << t << "\t" << Ctfinalfin << "\t" << Atfinalfin << "\t"  << MMnufinalfin << "\t"<< LAfinalfin << "\t" << finstate  << "\t" << exttime <<endl;

			}			
		}
				
				
	
	gsl_odeiv_evolve_free (evol);
	gsl_odeiv_control_free (control);
	gsl_odeiv_step_free (st);
	
	fout.close();
return 0;
}



	
	

     
