//List.h   Rui.   Dec-1-2017

#pragma once
#ifndef _LIST_H_
#define _LIST_H_
#include<stdlib.h>
#include<string>
#include<stddef.h>
#include<iostream>

using namespace std;

class FEMComponent; // supclass

class List
{
	// define an array for containing elements, nodes, materials
	// loads  and load-time functions
typedef	enum { FALSE };
protected:
	int size;
	FEMComponent** values;

public:
	List(int){ };
	~List() { };

	FEMComponent* at(int i) { return values[i - 1]; }   // return the value of buffer i
	int			giveSize() { return size; }             // check the size of the list
	void        growTo(int );                           // increase the storage capacity
	int         isEmpty() { return (size == 0); }       // check the empty list
	bool		includes(int);                          // check include
	int			isNotEmpty() { return (size != 0); }    // check not empty list
	void		printYourself();                        // print out the number in the list
	void		put(int, FEMComponent*);                // push_back the (nodes, materials, elements) into the list  
};
#endif // !_STORAGE_H_