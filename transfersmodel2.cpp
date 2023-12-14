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


int equationsros (double t, const double y[], double f[], void *params)
{
    double* p = (double *) params;
      
    double m = p[0];
    double d = p[1];
    double e = p[2];
    double l = p[3];
    double rn1 = p[4];
    double rn2 = p[5];
    double Kn1 = p[6];
    double Kn2 = p[7];
    double Yn1 = p[8];
    double Yn2 = p[9];
    double rc1 = p[10];
    double rc2 = p[11];
    double Kc1 = p[12];
    double Kc2 = p[13];
    double Yc1 = p[14];
    double Yc2 = p[15];
    double beta1 = p[16];
    double beta2 = p[17];
    double alpha1 = p[18];
    double alpha2 = p[19];
 	
	
	// y[0] = Ct, y[1] = At, y[2] = MM nutrient, y[3] = LA, y[4] = ROS
	
	
	f[0] = rn1*y[2]/(y[2] + Kn1)*y[0] + rc1*y[3]/(y[3] + Kc1)*y[0] - beta1*y[0]*y[4];
	
	f[1] = rn2*y[2]/(y[2] + Kn2)*y[1] + rc2*y[3]/(y[3] + Kc2)*y[1] - beta2*y[1]*y[4];
	
	f[2] = - 1/Yn1*rn1*y[2]/(y[2] + Kn1)*y[0] - 1/Yn2*rn2*y[2]/(y[2] + Kn2)*y[1];
		
	f[3] =  -1/Yc1*rc1*y[3]/(y[3] + Kc1)*y[0] -1/Yc2*rc2*y[3]/(y[3] + Kc2)*y[1] - 1/m*(d*y[3] + e*y[4]*y[3]);
	
	f[4] = d*y[3] + e*y[4]*y[3] - l*y[4] - alpha1*y[0]*y[4] - alpha2*y[1]*y[4] ;
	
	return GSL_SUCCESS;
	
	
}




int main(int argc, char *argv[])
{


    double m = 0.879171;
    double d = 0.1139162;
    double e = 1.941401;
    double  l = 0.346206;
    
    double rn1 = 5.518182;
    double rn2 = 126.4633;
    double Kn1 = 0.00103794;
    double Kn2 = 0.3570926;
    double Yn1 = 407993404;
    double Yn2 = 326741920;
    
    
    double rc1 = 3.494828;
    double rc2 = 2.368324;
    double Kc1 = 0.04394823;
    double Kc2 = 0.04856385;
    double Yc1 = 5788433873;
    double Yc2 = 1e10;
    double beta1 = 8.580369;
    double beta2 = 20.0305;
    double alpha1 = 4.612671e-06;
    double alpha2 = 0.0;
   

	
	 vector<double> constparam;
	 constparam.push_back(m);
	 constparam.push_back(d);
	 constparam.push_back(e);
	 constparam.push_back(l);
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
	 constparam.push_back(alpha1);
	 constparam.push_back(alpha2);
	 

	 
	vector<double> vecparams = constparam;
	int tailleparams = vecparams.size();
	double params[tailleparams]; 
	for (int i = 0; i<tailleparams; ++i) {
		params[i] = vecparams[i];
	}


	
	double bact1init = 1e5; // Ct
    double bact2init = 1e5; //At
    double rosini = 0.01; 
    double LAinit = 0.75;   
    double compini = 0.01;
    
    unsigned int nbequ = 5;
        
    const gsl_odeiv_step_type * T = gsl_odeiv_step_rkf45;
  
    gsl_odeiv_step * st = gsl_odeiv_step_alloc (T, nbequ); // creation an instance of the stepping function
	gsl_odeiv_control * control = gsl_odeiv_control_y_new (1e-8, 0.0); //
	gsl_odeiv_evolve * evol = gsl_odeiv_evolve_alloc (nbequ); // instance of an evolution function for a system of 1 dimensions
	gsl_odeiv_system sys = {equationsros, NULL, nbequ, &params}; // the Jacobian is useless with this method : use NULL instead
	
	double tbatch = 3.0; //number of days in which the bacteria grow before being diluted
	double tinteg = tbatch;
	int nbtransfer = 5; //number of transfers: 5, 500...
	double t = 0.0; // time span of integration
	double hfabs = 1e-10; // fabsolute accuracy;
	double y[nbequ]; 
	
	string typ; //type of culture: monoculture, coculture
	if (bact1init > 0 && bact2init > 0) {typ = "cocult";} 
	if (bact1init > 0 && bact2init <= 0) {typ = "monoCt";}
	if (bact1init <= 0 && bact2init > 0) {typ = "monoAt";}
	
	stringstream randIDpar;
	randIDpar << rand();

	
	string coextxt = "Model2_" + typ + "_" + to_string(nbtransfer) +".txt";
		
	ofstream fout(coextxt.c_str(), ios::out);
	fout << "LA_initialconc"  <<"\t" << "dilu" << "\t"   << "nbtransfer" << "\t" << "finaltime" << "\t" << "Ct" << "\t" << "At" << "\t"  << "MMnu" << "\t"<< "LA" << "\t" << "ROS" << "\t" << "finstate" <<  "\t" << "exttimeat" << endl;


	bool valnan = 0;
	
	int pasenregi = 100;
	

	// creating the parameter grid
	
	double dilumin;
	double dilumax; 
	vector<double> vecdilu;	
	double pasdilu;	
	
	
	dilumin = 10;
	dilumax = 200;
	pasdilu= 10;
	double dilu = dilumin;
	
	while (dilu < dilumax + pasdilu) {
		vecdilu.push_back(dilu);
		dilu = dilu + pasdilu;	
	}
	

	double lamin;
	double lamax; 
	vector<double> vecla;	
	double pasla;	
	
	lamin = 0.05;
	lamax = 1;
	pasla= 0.05;
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
			y[4]= rosini;
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
			double Rosfinalfin = 0.0;
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
					
					
					
					if (is_inf_or_nan(y[0]) or is_inf_or_nan(y[1])  or is_inf_or_nan(y[2]) or is_inf_or_nan(y[3])  or is_inf_or_nan(y[4]) ) {
						valnan=1;
						cout << "NaN!" << endl;
						break;}        
							 
					tpoint++;
				} 
				
				Atfinalfin = y[1];
				Ctfinalfin = y[0];
				Rosfinalfin = y[4];
				MMnufinalfin = y[2];
				LAfinalfin =y[3];
				
				y[0] = y[0] /dilu;
				y[1] = y[1] / dilu;
				y[2] = compini;
				y[3] = LAinit;
				y[4] = rosini;
				
				Ctfinal = y[0];
				Atfinal = y[1];	
				
				
				if (Atfinal < eqthr && notyetext) {
					exttime = transfer; 
					notyetext = 0;}
				
				if ( abs(Ctfinal - Ctprev) < eqthr and abs(Atfinal - Atprev) < eqthr ) { succev++;}	
				
				transfer++;
				tinteg = tinteg + tbatch;
				
			}

				string finstate;
			if (notyetext) {exttime = nbtransfer;} 


			 if (Atfinalfin > eqthr and Ctfinalfin > eqthr) { finstate = "coex";}
			 if (Atfinalfin < eqthr and Ctfinalfin > eqthr) { finstate = "Ct";}
			 if (Atfinalfin > eqthr and Ctfinalfin < eqthr) { finstate = "At";}
			 if (Atfinalfin < eqthr and Ctfinalfin < eqthr) { finstate = "ext";}
				
			
				
				cout << t << endl;
				fout << LAinit  <<"\t" << dilu  << "\t" << transfer << "\t" << t << "\t" << Ctfinalfin << "\t" << Atfinalfin << "\t"  << MMnufinalfin << "\t"<< LAfinalfin<< "\t" << Rosfinalfin << "\t" << finstate << "\t" << exttime <<endl;

			}			
		}
				
				
	
	gsl_odeiv_evolve_free (evol);
	gsl_odeiv_control_free (control);
	gsl_odeiv_step_free (st);
	
	fout.close();
return 0;
}



	
	

     
