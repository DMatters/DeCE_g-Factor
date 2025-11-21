/******************************************************************************/
/**     DeCE Maxwellian Average and Westcott Factor                          **/
/******************************************************************************/

#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <sstream>
#include <vector>
#include <utility>

using namespace std;

#include "dece.h"
#include "terminate.h"

const int    subdiv    = 100;
const double PI        = 3.14159265358979323846;
const double BOLTZMANN = 8.6173324e-5; // Boltzman Constant [eV/K]
const double Ethermal  = 0.02526;      // Thermal energy [eV]
const double eps       = 1.0e-99;

static double specaverage (const int, double *, double *, const double, const double, const double, const bool);
static double arbspecaverage (const int, double *, double *, vector<double>, vector<double>, const double, const double, const bool);
static inline double maxwellian (const double, const double);

static bool firstcall = true;


/******************************************************/
/*      Import neutron spectrum from .csv file        */
/******************************************************/

//Function takes two-column .csv file for spectrum, skipping header line
vector< pair<double, double> > spectrumCSV(const string& filename)
{
    vector< pair<double, double> > data;
    ifstream file(filename);
    
    if (!file.is_open()) {
        cerr << "Failed to open file: " << filename << endl;
        return data; // Return empty vector if file cannot be opened
    }

    // Skip the header row
    string header;
    if (!getline(file, header)) {
        cerr << "File is empty or cannot read header." << endl;
        return data;
    }

    string line;
    while (getline(file, line)) {
        istringstream iss(line);
        string col1, col2;
        
        if (getline(iss, col1, ',') && getline(iss, col2, ',')) {
            try {
                double val1 = stod(col1);
                double val2 = stod(col2);
                data.emplace_back(val1, val2);
            } catch (const invalid_argument& e) {
                cerr << "Invalid argument in line: " << line << endl;
            } catch (const out_of_range& e) {
                cerr << "Out of range error in line: " << line << endl;
            }
        } else {
            cerr << "Line does not contain two columns: " << line << endl;
        }
    }

    file.close();
    return data;
}

/*************************************************/
/*    Interpolate to define flux at any energy   */
/*************************************************/

double fluxInterpolate(const vector<double>& x, const vector<double>& y, double x_val) {
    if (x.empty() || y.empty() || x.size() != y.size()) {
        return 0; // Return 0 if input is invalid
    }

    // Find the position where x_val would fit in the sorted x array
    auto it = lower_bound(x.begin(), x.end(), x_val);

    if (it == x.end() || (it == x.begin() && *it != x_val)) {
        return 0; // Return 0 if x_val is out of bounds
    }

    if (it == x.begin()) {
        return y[0]; // x_val is less than or equal to the first element
    }

    if (it == x.end()) {
        return y.back(); // x_val is greater than the last element
    }

    size_t index = it - x.begin() - 1;
    double x0 = x[index];
    double y0 = y[index];
    double x1 = x[index + 1];
    double y1 = y[index + 1];

    // Perform linear interpolation
    return y0 + (y1 - y0) * (x_val - x0) / (x1 - x0);
}

/*****************************************************/
/*      Calculate Westcott g Factor                  */
/*****************************************************/
void DeceMaxwellian(ENDFDict *dict, ENDF *lib[], const double temperature, string ope)
{
  int k0 = dict->getID(3,102);
  if(k0 < 0){ message << "Capture channel not found"; TerminateCode("DeceWestcottFactor"); }

  bool westcott = (ope == "westcottfactor");

  double t = temperature;
  /*** when negative, the temperature is in K */
  if(temperature < 0.0) t *= -BOLTZMANN;

  /*** copy pointwise cross sections */
  int     n = lib[k0]->rdata[0].n2;
  double *x = new double [n];
  double *y = new double [n];
  for(int i=0 ; i<n ; i++){
    x[i] = lib[k0]->xdata[2*i  ];
    y[i] = lib[k0]->xdata[2*i+1];
  }

  /*** thermal cross section */
  double sig0 = ENDFInterpolation(lib[k0], Ethermal, false, 0);

  Record head = lib[k0]->getENDFhead();
  double awr  = head.c2;
  double kcms = awr/(awr+1.0); // convert energy into CMS for MACS
  if(westcott) kcms = 1.0;

  /*** print header part */
  cout << setprecision(7) << setiosflags(ios::scientific);
  if(firstcall){
    if(westcott){
      if(temperature <  0.0) cout << "# MacsTemp[K]    Westcott-g" << endl;
      else                   cout << "# MacsTemp[eV]   Westcott-g" << endl;
    }
    else{
      if(temperature <  0.0) cout << "# MacsTemp[K]    MACS[b]" << endl;
      else                   cout << "# MacsTem[eV]    MACS[b]" << endl;
    }
    firstcall = false;
  }

  /*** when temperature is not given, calculate at default grid 
  if(temperature == 0.0){
    if(westcott){
      double t0 =  10.0; // from 10 K to 500 K
      double t1 = 500.0;
      double dt =  10.0;

      for(int i=0 ; ; i++){
        double tx = t0 + i * dt;
        double gfac = specaverage(n,x,y,tx*BOLTZMANN,sig0,kcms,westcott);
        cout << setw(15) << tx << setw(15) << gfac << endl;
        if(tx >= t1) break;
      }
    }
    else{
      double e0 = 1e+4; // from 10 keV to 1 MeV
      double e1 = 1e+6;
      double e2 = 1e+5;
      double de = 1e+4;

      double ex = e0;
      for(int i=0 ; ; i++){
	i = i // quiet warning during compile
        double macs = specaverage(n,x,y,ex,sig0,kcms,westcott);
        cout << setw(15) << ex << setw(15) << macs << endl;
        if(ex >= e1) break;
        if(ex >= e2) de = 5e+4;
        ex += de;
      }
    }
  }
  */

  /*** When temperature is not given, go with user-defined neutron energy spectrum */
  if(temperature == 0.0){
    string filename;
    cout << "You have elected to apply a user-defined neutron energy spectrum instead of a Maxwellian. Enter the name of the CSV file for the spectrum: " <<endl;
    cin >> filename;
    
    vector< pair<double, double> > data = spectrumCSV(filename);

    // Determine the length of the pair structure
    size_t length = data.size();
    vector<double> x_spectrum(length);
    vector<double> y_spectrum(length);

    // Copy data into x_spectrum and y_spectrum
    for (size_t i = 0; i < length; ++i) {
        x_spectrum[i] = data[i].first;
        y_spectrum[i] = data[i].second;
    }
    
    double macs = arbspecaverage(n,x,y,x_spectrum,y_spectrum,sig0,kcms,westcott);
    cout << setw(15) << "Westcott g-factor, user-defined neutron energy spectrum:" << setw(15) << macs <<endl;
  }
  
  /*** MACS or Westcott g-factor at given Maxwellian temperature */
  else{
    double macs = specaverage(n,x,y,t,sig0,kcms,westcott);
    if(temperature < 0.0){
      cout << setw(15) << -temperature << setw(15) << macs <<endl;
    }
    else{
      cout << setw(15) <<  temperature << setw(15) << macs <<endl;
    }
  }

  delete [] x;
  delete [] y;
}

/**********************************************************/
/*      User-defined Arbitrary Spectrum Average           */
/**********************************************************/
double arbspecaverage(const int n, double *x, double *y, vector<double> x_spec, vector<double> y_spec, const double s0, const double c1, const bool westcott)
{
  double c0 = 1.0 / (s0 * sqrt(Ethermal));

  double s1 = 0.0, s2 = 0.0;
  for(int i=0 ; i<n-1 ; i++){
    double d = (x[i+1] - x[i]) * c1; // energy width
    if(d == 0.0) continue;

    /*** sub-divide each of energy bins */
    double d1 = d/(double)subdiv;
    double y1 = (y[i+1] - y[i])/(double)subdiv;

    double e1 =  x[i]       * c1;
    double e2 = (x[i] + d1) * c1;

    double u1 = y[i];
    double u2 = y[i] + y1;

    /*** sqrt(E) x sigma */
    if(westcott){
      u1 *= sqrt(e1);
      u2 *= sqrt(e2);
    }
    
    /*** User-defined (interpolated) spectrum at both bin sides */
    double w1 = fluxInterpolate(x_spec,y_spec,e1);
    double w2 = fluxInterpolate(x_spec,y_spec,e2);

    double p1 = 0.0; // sum of cross section x spectrum
    double p2 = 0.0; // sum of spectrum

    for(int j=1 ; j<=subdiv ; j++){

      p1 += (u1*w1 + u2*w2) * d1 * 0.5;
      p2 += (   w1 +    w2) * d1 * 0.5;

      e1 = e2; e2 = x[i] + d1*(j+1);
      u1 = u2; u2 = y[i] + y1*(j+1); if(westcott){ u2 *= sqrt(e2); }
      w1 = w2; w2 = fluxInterpolate(x_spec,y_spec,e2);

      if(w1 == 0.0 && w2 == 0.0) break;
    }

    s1 += p1;
    s2 += p2;
  }

  double macs = 0.0;
  if(s2 > 0.0){
    macs = s1 / s2;
    if(westcott) macs *= c0;
    else         macs *= 2.0/sqrt(PI);
  }

  return(macs);
}

  
/***********************************************/
/*      Maxwellian Average                     */
/***********************************************/
double specaverage(const int n, double *x, double *y, const double t, const double s0, const double c1, const bool westcott)
{
  double c0 = 1.0 / (s0 * sqrt(Ethermal));

  double s1 = 0.0, s2 = 0.0;
  for(int i=0 ; i<n-1 ; i++){
    double d = (x[i+1] - x[i]) * c1; // energy width
    if(d == 0.0) continue;

    /*** sub-divide each of energy bins */
    double d1 = d/(double)subdiv;
    double y1 = (y[i+1] - y[i])/(double)subdiv;

    double e1 =  x[i]       * c1;
    double e2 = (x[i] + d1) * c1;

    double u1 = y[i];
    double u2 = y[i] + y1;

    /*** sqrt(E) x sigma */
    if(westcott){
      u1 *= sqrt(e1);
      u2 *= sqrt(e2);
    }
    
    /*** Maxwellian spectrum at both bin sides */
    double w1 = maxwellian(e1,t);
    double w2 = maxwellian(e2,t);

    double p1 = 0.0; // sum of cross section x spectrum
    double p2 = 0.0; // sum of spectrum

    for(int j=1 ; j<=subdiv ; j++){

      p1 += (u1*w1 + u2*w2) * d1 * 0.5;
      p2 += (   w1 +    w2) * d1 * 0.5;

      e1 = e2; e2 = x[i] + d1*(j+1);
      u1 = u2; u2 = y[i] + y1*(j+1); if(westcott){ u2 *= sqrt(e2); }
      w1 = w2; w2 = maxwellian(e2,t);

      if(w1 == 0.0 && w2 == 0.0) break;
    }

    s1 += p1;
    s2 += p2;

    if(p1 < eps && p2 < eps) break;
/*
    cout << setprecision(5) << setiosflags(ios::scientific);
    cout << setw(5) << i;
    cout << setw(14) << x[i];
    cout << setw(14) << y[i];
    cout << setw(14) << s0 * sqrt(Ethermal / x[i]);
    cout << setw(14) << maxwellian(x[i]*c1,t);
    cout << setw(14) << p1;
    cout << setw(14) << p2;
    if(s2 > 0.0) cout << setw(12) << s1/s2 * c0;
    cout << endl;
*/
  }

  double macs = 0.0;
  if(s2 > 0.0){
    macs = s1 / s2;
    if(westcott) macs *= c0;
    else         macs *= 2.0/sqrt(PI);
  }

  return(macs);
}


inline double maxwellian(const double e, const double t)
{
  double w = e * exp(-e/t);
  if(w < eps) w = 0.0;
  return(w);
}
