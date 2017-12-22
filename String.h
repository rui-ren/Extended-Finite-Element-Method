// This file is used for comparision, concatenation, conversion

#ifndef STRING_H
#define STIRNG_H

#include<string>
#include "mathfem.h"
#include<stdio.h>
#include<stdlib.h>
#include<ctype.h>
using namespace std;

char* concatenate(char*, int, char*);
char* ctos(char, char*);
bool isInteger(char*);
bool isNumber(char*);
char* itos(int, char*);

#endif

char* concatenate(char* s , int i, char* target)
{
	char s2[8];

	target[0] = '\0';
	strcat(target, s);			//string concatenation
	strcat(target, itos(i, s2));
	return target;
}

char* ctos(char c, char* s)
{
	s[0] = c;
	s[1] = '\0';
	return s;
}

bool isInteger(char* s)
{
	// find the integer in s function, else return false.
	int  len, i;
	char c;

	len = strlen(s);
	if (!len)
		return false;

	c = s[0];
	if (c == '+' || c == '-')
		i = 1;
	else if (isdigit(c))
		i = 0;
	else
		return false;

	for (; i<len; i++)
	if (!isdigit(s[i]))
		return false;
	return true;
}

bool isNumber(char* s)
{
	int  len, i;
	char c;

	len = strlen(s);
	if (!len)
		return false;

	c = s[0];                             // look for a sign
	if (c == '+' || c == '-')
		i = 1;
	else if (isdigit(c))
		i = 0;
	else
		return false;


	for (i = i + 1; i<len; i++)              // process integer part
	if (!isdigit(c = s[i])) {
		if (c == '.')
			goto decimal;
		else if (c == 'e' || c == 'E')
			goto exponent;
		else
			return false;
	}
	return true;

decimal:
	for (i = i + 1; i<len; i++)
	if (!isdigit(c = s[i])) {
		if (c == 'e' || c == 'E')
			goto exponent;
		else
			return false;
	}
	return true;

exponent:
	c = s[i + 1];
	if (c == '+' || c == '-')
		i++;
	if (i == len - 1)
		return false;
	for (i = i + 1; i<len; i++)
	if (!isdigit(s[i]))
		return false;
	return true;
}

char*  itos(int n, char* s)
// A non member function, copied from Kernighan & Ritchie, p64. Implemen-
// ted for portability (the PC function itoa(int,char*,int) is not ANSI
// C). Returns s, the string conversion of n.
{
	int  i, j, sign;
	char c;

	if (n == 1)
		return strcpy(s, "1");
	if (n == 2)
		return strcpy(s, "2");

	if ((sign = n) < 0)
		n = -n;
	i = 0;
	do
	s[i++] = n % 10 + '0';
	while ((n /= 10) > 0);
	if (sign < 0)
		s[i++] = '-';
	s[i] = '\0';

	// reverse 's'
	j = strlen(s) - 1;
	for (i = 0; i<j; i++, j--) {
		c = s[i];
		s[i] = s[j];
		s[j] = c;
	}

	return s;
}