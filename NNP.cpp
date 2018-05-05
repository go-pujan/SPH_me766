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
 
// A structure to represent node of kd tree
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
	vector<vector<float> > final(0,vector<float>(2));

	struct Node *root = NULL;
    float points[7][2] = {{3.25, 6.64}, {17.18, 15.58}, {13.14, 15.96}, {6.36, 12.58}, {9.41, 1.00}, {2.58, 7.2}, {10.12, 19.78}};
 
    int n = 7;
    float radius = 6;
    for (int i=0; i<n; i++)
       root = insert(root, points[i]);

    float point1[] = {3.25, 6.64};
 	run(root, point1, 0, radius, final);
 	
 	int end = final.size();
 	if(end == 0)
 	{
 		cout<<"No such point";
 		return 0;
 	}
 	for(int i=0; i<end; i++)
 	{
 		cout<<final[i][0]<<endl;
 		cout<<final[i][1]<<endl;
 	}
 
    return 1;
}