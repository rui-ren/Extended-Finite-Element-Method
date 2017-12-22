#include"Geometry.h"
#include<iostream>
#include<vector>
#include<list>
#include<map>

using namespace std;


class Mesh
{
	Mesh();	 //constructor
	~Mesh(); // destructor
public:
	vector<int> Elem_Num;
	vector<int> Elem_Node;
};


double Elem_Node(vector<int>& ElemNum)
{
	vector<int> Elem_Node;
	double whole_Elem;
	whole_Elem = ElemNum[0] * ElemNum[1];
	vector<int> Elem_Node0, Elem_Node1, Elem_Node2, Elem_Node3;

	/////////////////////////////////// Element and Node Association/////////////////////

	for (int i = 0; i < whole_Elem; i++)
	{
		vector<int> Elem_Num.push_back(i);
		for (int j = 0; j < ElemNum[0]; j++)
		{
			Elem_Node0.push_back(i + 1 + 1);
			Elem_Node1.push_back(i + 1 + ElemNum[0] * j + 1);
			Elem_Node2.push_back(i + 1 + ElemNum[0] * j);
			Elem_Node3.push_back(i + 1);
		}
	}
};

