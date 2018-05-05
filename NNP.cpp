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
#include <ctime>
#include <iostream>

using namespace std;

float distance(float a[], float b[])
{
    float dist = 0;
    float A = pow((a[0]-b[0]),2);
    float B = pow((a[1]-b[1]),2);
    dist = pow((A+B),0.5);
    return dist;
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


//______________________________KD Tree______________________________________// 
//A structure to represent node of kd tree
struct Node
{
    float data[2]; // To store k dimensional point; // To store k dimensional data
    Node *left, *right;
};
 
// A method to create a node of K D tree
struct Node* newNode(float input[])
{
    struct Node* temp = new Node;
 
    for (int i=0; i<2; i++)
       temp->data[i] = input[i];
 
    temp->left = temp->right = NULL;
    return temp;
}

// Inserts a new node and returns root of modified tree
// The parameter depth is used to decide axis of comparison
Node *insertRec(Node *root, float data[], unsigned depth)
{
    // Tree is empty?
    if (root == NULL)
       return newNode(data);
 
    // Calculate current dimension (cd) of comparison
    unsigned cd = depth % 2;
 
    // Compare the new data with root on current dimension 'cd'
    // and decide the left or right subtree
    if (data[cd] < (root->data[cd]))
        root->left  = insertRec(root->left, data, depth + 1);
    else
        root->right = insertRec(root->right, data, depth + 1);
 
    return root;
}
 
// Function to insert a new data with given data in
// KD Tree and return new root. It mainly uses above recursive
// function "insertRec()"
Node* insert(Node *root, float data[])
{
	float forward[2];
	forward[0]=data[0];
	forward[1]=data[1];
    return insertRec(root, forward, 0);
}

// A utility method to determine if two datas are same
// in K Dimensional space
bool aredatasSame(float data1[], float data2[])
{
    // Compare individual datainate values
    for (int i = 0; i < 2; ++i)
        if (data1[i] != data2[i])
            return false;
 
    return true;
}

vector<vector<float> > run(Node* root, float data[], unsigned depth, float radius, vector<vector<float> > &final)
{
	
	if(root == NULL)
	{
		return final;
	}

	if(distance(root->data, data)> radius || aredatasSame(root->data, data))
	{
		unsigned cd = depth % 2;
 
    // Compare data with root with respect to cd (Current dimension)
    	if (data[cd] < root->data[cd])
    	{
    		return run(root->left, data, depth+1, radius, final);
    	}

    	return run(root->right, data, depth+1, radius, final);
	}

	if(distance(root->data, data)<radius)
	{
		vector<float> push(2);
		push[0]= root->data[0];
		push[1]= root->data[1];

		final.push_back(push);
		unsigned cd = depth % 2;
 
    // Compare data with root with respect to cd (Current dimension)
    	if (data[cd] < root->data[cd])
    	{
    		return run(root->left, data, depth+1, radius, final);
    	}

    	return run(root->right, data, depth+1, radius, final);

	}
}

int main()
{
    int N=1024;
    string input_params, output_dir = "Time";
    int error;

    int numpts, Nw, local_size, saveFreq;
    float box_size_x, box_size_y, rho0, viscosity, velocity, total_t, dt, h;

    /*readParams(input_params, output_dir, &numpts, &Nw, &box_size_x, &box_size_y,
               &rho0, &viscosity, &velocity, &total_t,
               &local_size, &dt, &h, &saveFreq);*/
    numpts = 1024;
    box_size_x = 1; box_size_y=1; rho0=100;total_t=0.1;dt=0.01;h=0.01;saveFreq=1;local_size=256;
    velocity = 100; viscosity=0.1; Nw=1024;
    int numrep = (total_t/dt);
    float c0 = 10;
    float dx = h/1.1, m = rho0 * dx * dx;
    vector<vector<float> > x(numpts,vector<float> (2.0)), xw(Nw,vector<float> (2.0)), vw(Nw,vector<float> (2.0)), v(numpts,vector<float> (2.0));
    vector<float> r(numpts), rw(Nw), p(numpts), pw(Nw);
    set_ic(x, xw, v, vw, r, rw, dx, box_size_x, box_size_y);

	vector<vector<float> > final(0,vector<float>(2));//The vector of nearest neighbours

	struct Node *root = NULL;
    float points[x.size()][2];

    for (int i = 0; i < x.size(); ++i)
    {
        points[i][0]=x[i][0];
        points[i][1]=x[i][1];
    }

    int n = 7;
    float radius = 2;
    for (int i=0; i<n; i++)
       root = insert(root, points[i]);

 	run(root, points[265], 0, radius, final);
 	
 	int end = final.size();
 	if(end == 0)
 	{
 		cout<<"No such point";
 		return 0;//Can be used as a condition for checking if no nearest neighbours
 	}
 	for(int i=0; i<end; i++)
 	{
 		cout<<final[i][0]<<", "<<final[i][1]<<endl;
 	}
 
    return 1;
}