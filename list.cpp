// List.cpp file   Rui. Dec-1-2017

#include "list.h"
#include<iostream>
#include<string>
#include<stdlib.h>

using namespace std;

List::List(int s)   // constructor creates a list of size s.
{
	register int i;       // register i which is faster than memory to access
	FEMComponent** p;
	if (size)
	{
		values = new FEMComponent*[size];
		p = values;
		i = size;
		while (i--)        
			*p++ = NULL;         //add the compounts(nodes, materials, element into list)
	}
	else
		values = NULL;
}

List::~List()		// Destructor
{
	int i = size;
	if (size)
	{
		while (i--)         // delete the components in the list.
		{
			delete values[i];
			delete[] values;
		}
	}
};

void List::growTo(int newSize)
{
		register int i;
		FEMComponent **newValues, **p1, **p2;
		if (newSize <= size)   
		{
			cout << newSize << endl;   // output the newSize
			throw logic_error("New list size is not larger enough");
		}

		// expand the list p1;
		newValues = new FEMComponent*[newSize];
		p1 = values;
		p2 = newValues;
		for (int i = 0; i < size; i++)
			*p2++ = *p1++;     // give the origin value to p2
		for (int j = size; j < newSize; j++)
			*p2++ = NULL;     // Null the new value
		if (values)
			delete[] values; // delete memory

		// set the value and size of the list
		values = newValues;
		size = newSize;
}


// check i th components in the list.
bool List::includes(int i)
	{
		if (i > size)			// check whether the index within the list 
			return false;
		else
			return true;       // within in the list
	}

// output the element/ nodes/ materials in the list
	void List::printYourself()
	{
		register int i;
		for (int i = 0; i < size; i++)
		{
			if (values[i - 1] == NULL)
				cout << i << endl;
			else
				cout << (long int)values[i - 1] << endl;
		}
	}

// push_back material, nodes, element in the list.
	void List::put(int i, FEMComponent* anObject)
		//stores anObject at position i, enlarge the receiver if too small
	{
		if (size < i)
		{
			this->growTo(i);
			values[i - 1] = anObject;        // use the push_back method
		}
	}
