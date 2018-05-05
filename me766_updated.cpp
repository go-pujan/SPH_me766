#include <cstdlib>
#include <algorithm>
#include <functional>
#include <math.h>
#include <random>
#include <time.h>
#include <stdio.h>
#include <string>
#include <vector>
#include <random>
#include <fstream>
#include <iostream>
//#include "cl_helpers.h"

using namespace std;
float length(vector<float> a){
	float A = pow(a[0],2);
	float B = pow(a[1],2);
	return pow((A+B),0.5);
}

float distance(vector <float>a, vector<float>b)
{
	float dist;
	float A = pow((a[0]-b[0]),2);
	float B = pow((a[1]-b[1]),2);
	dist = pow((A+B),0.5);
	return dist;
}
float dot(vector <float>a, vector<float>b){
	return a[0]*b[0]+a[1]*b[1];
}

float time_step(float cfl, float h0, float c0)
{
    float dt=cfl*h0/c0 ;  
    return dt;
}

float kernel_cubic(vector<float> xi, vector<float> xj, float h, float dist)
{
	float pi = 3.14159;
    float q = dist/h;
    float W = 0;
    if (q <= 1.&& q >= 0)
        W += 10 / (7 * pi * h * h) * (1. - 3 / 2 * q * q * (1 - q / 2));
    if (q > 1. && q < 2.)
        W += 10 / (28 * pi * h * h) * pow((2 - q),3);
    return W;
}

vector<float> kernel_derivative(vector<float> xi, vector<float> xj, float h, float dist)
{	
	float pi = 3.14159;
    float q = dist / h;
    float dwdq = 0;
    if (q <= 1.)
        dwdq = (9 / 4 * q - 3) * 10 / (7 * pi * h * h);
    if (q > 1. && q < 2.)
        dwdq = -7.5 * (2 - q) * (2-q) / (7 * pi * q * h * h);

    vector<float> dW (2);
    dW[0] = dwdq * (xi[0] - xj[0]) / h / h;
    dW[1] = dwdq * (xi[1] - xj[1]) / h / h;

    return dW;
}

float art_visc(vector<float> x_i, vector<float> x_j, float r_i, float r_j, vector<float> v_i, vector<float> v_j, float p_i, float p_j, float h)
{
    float alpha = 1;
    float beta = 1;
    vector<float> x(2);
    x[0] = (x_i[0] - x_j[0]);
    x[1] = (x_i[1] - x_j[1]);
    float neta = 0.01 * h;
    float pia = 0;
    vector <float> v(2);
    v[0] = (v_i[0] - v_j[0]);
    v[1] = (v_i[1] - v_j[1]);

    if (dot(x, v) <= 0)
    {
        float ca = (sqrt(1.4 * p_i / r_i) + sqrt(1.4 * p_j / r_j)) / 2;
        float ra = (r_i + r_j) / 2;
        float mu = h * dot(v, x) / (pow(length(x), 2) + neta*neta);
        pia = (-alpha * ca * mu + beta * mu * mu) / ra;
    }
    return pia;
}


 void UPDATE_POS(int iterator, vector<vector<float> >& x, vector<vector<float> > &v, vector<float> &r, float m, float h, int N, float dt)
{

    vector<float> tmp(2);
    tmp[0]=0.0;
    tmp[1]=0.0;
    float e_con = 0.5;
    int i = iterator;
    for(int j = 0; j < N; j++)
    {   
        float dist = distance(x[i], x[j]);
        if (dist < 2.0)
        {
            float W = kernel_cubic(x[i], x[j], h, dist);
            tmp[0] += e_con * (m * 2 / (r[j] + r[i])) * (v[j][0] - v[i][0]) * W;
            tmp[1] += e_con * (m * 2 / (r[j] + r[i])) * (v[j][1] - v[i][1]) * W;
        }
    }

    x[i][0] += (tmp[0] + v[i][0]) * dt;
    x[i][1] += (tmp[1] + v[i][1]) * dt;
}

// N = number of fluid particles. Nw = wall particles. Launch one kernel per fluid particle
 void SUMDEN(int iterator, vector<vector<float> > &x, vector<vector<float> > &xw, vector<float> &r, float m, float h, int N, int Nw)
{
	int i = iterator;
    r[i] = 0;

    for (int j=0; j<N; j++){

        float dist = distance(x[i], x[j]);
        if(dist < 2.0){
            r[i] += m * kernel_cubic(x[i], x[j], h, dist);
        }   
    }

    for (int j=0; j<N; j++){

        float dist = distance(x[i], xw[j]);
        if(dist < 2.0){
            r[i] += m * kernel_cubic(x[i], xw[j], h, dist);
        }   
    }
}

 void DEN(int iterator, vector<vector<float> > &x, vector<vector<float> > &xw, vector<vector<float> > &v, vector<vector<float> > &vw, vector<float> &r, vector<float> &rw, float m, float h, float dt, int N, int Nw)
{
	int i = iterator;
    float dr = 0;
    vector <float> vec(2);
    vec[0]=0.0;
    vec[1]=0.0;
    for (int j=0; j<N; j++)
    {   
        float dist = distance(x[i], x[j]);
        if(dist < 2.0)
        {
    	   vec[0] = v[i][0] - v[j][0];
    	   vec[1] = v[i][1] - v[j][1];
            dr += 1 / r[j] * dot(vec, kernel_derivative(x[i], x[j], h, dist));
        }
    }
    vector <float> vect(2);
    vect[0]=0.0;
    vect[1]=0.0;
    for (int j=0; j<Nw; j++)
    {   
        float dist = distance(x[i], xw[j]);
        if(dist < 2.0)
        {
    	   vect[0] = v[i][0] - vw[j][0];
    	   vect[1] = v[i][1] - vw[j][1];
            dr += 1 / rw[j] * dot(vect, kernel_derivative(x[i], xw[j], h, dist));
        }
    }
    r[i] += dr * m * r[i] * dt;
}

// Launch one kernel per FLUID particle
 void INCOMP_P(int iterator, vector<float> &r,vector<float> &p, float c0, float rho0)
{	
	int i = iterator;
	float gamma = 1.4;
    float B = rho0 * c0 * c0 / gamma;
    double tmp = pow((double)r[i] / rho0, (double)gamma) - 1;
    p[i] = B * (tmp) ;
    
}


// Launch one kernel per FLUID particle
 void UPDATE_VEL(int iterator, vector<vector<float> > &x, vector<vector<float> > &xw, vector<float> &p, vector<float> &pw, vector<vector<float> > &v, vector<vector<float> > &vw,vector<float> &r, vector<float> &rw, float m, int N, int Nw, float dt, float h)
{
	int i = iterator;
    vector<float> dv (2);
	dv[0]=0.0;
    dv[1]=0.0;

    for (int j=0; j<N; j++)
    {
        if (distance(x[i], x[j]) < 2*h)
        {
            float dist = distance(x[i], xw[j]);
            float av = art_visc(x[i], x[j], r[i], r[j], v[i], v[j], p[i], p[j], h);
            vector<float> dW = kernel_derivative(x[i], x[j], h, dist);
        
            float calc = (p[j] / r[j] / r[j] + p[i] / r[i] / r[i] + av);

            dv[0] += - m * calc * dW[0];
            dv[1] += - m * calc * dW[1];
        }
    }

    for (int j=0; j<Nw; j++)
    {
        if (distance(x[i], xw[j]) < 2*h)
        {
            float dist = distance(x[i], xw[j]);
            float av = art_visc(x[i], xw[j], r[i], rw[j], v[i], vw[j], p[i], pw[j], h);
            vector<float> dW = kernel_derivative(x[i], xw[j], h, dist);
        
            float calc = (pw[j] / rw[j] / rw[j] + p[i] / r[i] / r[i] + av);

            dv[0] += - m * calc * dW[0];
            dv[1] += - m * calc * dW[1];
        }
    }
    v[i][0] += dv[0] * dt;
    v[i][1] += dv[1] * dt;
}

// Launch one kernel per wall particle
 void WALL(int iterator, vector<vector<float> > &x, vector<vector<float> > &xw, vector<vector<float> > &v, vector<vector<float> > &vw, vector<float> &p, vector<float> &pw, vector<float> &rw, float h, float rho0, float c0, int N)
{	
	int i = iterator;
	float gamma = 1.4;
    vector<float> num_v_w (2);
    num_v_w[0]=0.0;
    num_v_w[1]=0.0;
    float num_p_w =0;
    float den_w = 0;
	
    for (int j=0; j<N ; j++)
    {
        float dist = distance(x[i], xw[j]);
        if(dist < 2.0)
        {
            num_v_w[0] += v[j][0]*kernel_cubic(xw[i], x[j], h, dist);
            num_v_w[1] += v[j][1]*kernel_cubic(xw[i], x[j], h, dist);
            num_p_w += p[j]*kernel_cubic(xw[i], x[j], h, dist);
            den_w += kernel_cubic(xw[i], x[j], h, dist) + 1e-10;
        }
    }

    vw[i][0] = - num_v_w[0]/den_w;
    vw[i][1] = - num_v_w[1]/den_w;
    pw[i] = num_p_w/den_w;
	
    float B = rho0 * c0 * c0 / gamma;
    double tmp = pw[i]/B + 1; 
    double tmp1 = 1/gamma;
    double tmp2 =  pow((double) tmp, (double) tmp1);
    rw[i] = rho0 * (tmp2) ;
}
// ______________________________________________________________End of Kernel Functions____________________________________________//


void readParams(string filename, string &output_dir, int* numpts, int*Nw, float* box_size_x, float* box_size_y, float* density, float* viscosity, float* velocity, float* total_t, int* local_size, float *dt, float *h, int* saveFreq)
{

    int nf = 8; //Number of floating point parameters
    int nd = 4; //Number of integer parameters
    
    int *d_vars[] = { numpts, Nw, local_size, saveFreq };
    float *f_vars[] = { box_size_x,  box_size_y, density, viscosity, velocity, total_t, dt, h };

    string d_params[] = {"number of fluid particles : ", "number of wall particles : ",
                         "threads per block : ", "checkpoint save frequency : "};
    string f_params[] = {"box size x : ", "box size y : ", "density : ",
                         "viscosity : ", "velocity : ", "simulation time : ",
                         "time step : ", "kernel size : "};

    ifstream f(filename);
    string line;

    if(f.is_open())
    {
        while(getline (f,line) )
        {
             int found;
             char* pEnd;
             transform(line.begin(), line.end(), line.begin(), ::tolower);

             for (int i=0; i < nd; i++){
                 if (line.find(d_params[i]) != string::npos)
                     *(d_vars[i]) = (int)strtof(line.substr(d_params[i].size()).c_str(), &pEnd);
             }
             for (int i=0; i < nf; i++){
                 if (line.find(f_params[i]) != string::npos)
                     *(f_vars[i]) = strtof(line.substr(f_params[i].size()).c_str(), &pEnd);
             }
             if (line.find("output directory") != string::npos)
                 output_dir = line.substr(19);
        }

        f.close();
    }
    else
    {
        cout << filename << " does not exist. Exiting\n";
        exit(1);
    }
}

void set_ic(vector<vector<float> > &x, vector<vector<float> > &xw, vector<vector<float> > &v,
            vector<vector<float> > &vw, vector<float> &r, vector<float> &rw,
            float dx, float box_size_x, float box_size_y)
{
    uniform_real_distribution <double> randx(0,box_size_x);
    uniform_real_distribution <double> randy(0,box_size_y);
    uniform_real_distribution <double> rvx(1, -1);
    uniform_real_distribution <double> rvy(1, -1);
    default_random_engine rex, rey, revx, revy;
    rex.seed(time(NULL));
    rey.seed(time(NULL));
    revx.seed(time(NULL));
    revy.seed(time(NULL));

    for (int i=0; i < x.size(); i++)
    {
        /* x[i].s[0] = i * dx; */
        /* x[i].s[1] = i * dx; */
        v[i][0] = 0;
        v[i][1] = 0;
        x[i][0] = randx(rex);
        x[i][1] = randy(rey);
        /* v[i].s[0] = rvx(revx); */
        /* v[i].s[1] = rvy(revy); */
        r[i] = 1000;
    }
    for (int i=0; i < xw.size(); i++)
    {
        xw[i][0] = 1000;
        xw[i][1] = 1000;
        vw[i][0] = 1;
        vw[i][1] = 0;
        rw[i] = 100000;
    }
}
void saveCheckpoint(string i, string output_dir, vector<vector<float>> &x, vector<vector<float>> &v, vector<float> &r, vector<float> &p)
{
    int numpts = x.size();
    string command_str = "mkdir -p " + output_dir;
    const char* command = command_str.c_str();
    system(command);
    string filename = output_dir + "/" + i + ".csv";
    ofstream f(filename); 
    f << "Particle Number,X pos,Y pos,X vel,Y vel,Density,Pressure\n";
    //cout<<i<<endl;
    
    for (int i=0; i < x.size(); i++)
    {
        f << i << "," << x[i][0] << "," << x[i][1] << "," << v[i][0];
        f << "," << v[i][1] << "," << r[i] << "," << p[i] << endl;
    }
    f.close();
}
int main()
{   
    int N=1024;
	string input_params, output_dir="Fuck_Yeah";
	int error;

	int numpts, Nw, local_size, saveFreq;
    float box_size_x, box_size_y, rho0, viscosity, velocity, total_t, dt, h;

    /*readParams(input_params, output_dir, &numpts, &Nw, &box_size_x, &box_size_y,
               &rho0, &viscosity, &velocity, &total_t,
               &local_size, &dt, &h, &saveFreq);*/
    numpts = 1024;
    box_size_x = 1; box_size_y=1; rho0=100;total_t=0.1;dt=0.0005;h=0.01;saveFreq=1;local_size=256;
    velocity = 100; viscosity=1000; Nw=1024;


    float c0 = 10;
    float dx = h/1.1, m = rho0 * dx * dx;
    vector<vector<float> > x(numpts, vector<float> (2.0)), xw(Nw, vector<float> (2.0)), vw(Nw, vector<float> (2.0)), v(numpts, vector<float> (2.0));
    vector<float> r(numpts), rw(Nw), p(numpts), pw(Nw);
    set_ic(x, xw, v, vw, r, rw, dx, box_size_x, box_size_y);

    //Call functions directly
    for (int i=0; i<numpts; i++){
    	UPDATE_POS(i, x, v, r, m, h, numpts, dt);
    	SUMDEN(i, x, xw, r, m, h, numpts, Nw);
    	DEN(i, x, xw, v, vw, r, rw, m, h, dt, numpts, Nw);
    	INCOMP_P(i, r, p, c0, rho0);
    	UPDATE_VEL(i, x, xw, p, pw, v, vw, r, rw, m, numpts, Nw, dt, h);
    	WALL(i, x, xw, v, vw, p, pw, rw, h, rho0, c0, Nw);
        saveCheckpoint(to_string(i), output_dir, x, v, r, p);
	}
	saveCheckpoint("final", output_dir, x, v, r, p);


}
