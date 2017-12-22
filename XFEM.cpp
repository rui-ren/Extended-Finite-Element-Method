#include"Geometry.h"
#include"XFEM.h"
#include<vector>
#include"Mesh.h"
#include<list>
#include<cstdlib>

#include<iostream>

using namespace std;

int main()
{
	/////////////////////////////////// Plot the geometry of the domain///////////////////////////////////

	////////////////////////////////// Input the in formation of Mesh////////////////////////////////////
	vector<int> ElemNum;
	cout << "Input the Element Number in X direction" << endl;
	cin >> ElemNum[0];			// input the number of element in X direction
	cout << "Input the element number in Y direction" << endl;
	cin>>ElemNum[1] ;			// input the number of element in Y direction

	system("Pause");
	return 0;
}
