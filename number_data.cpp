// random_numbers.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <cstdlib>
#include "number_data.h"
using namespace std;

void fill_random_numbers (int arr [], const unsigned int SIZE)
{
	// Pre:  SIZE is no more than a million, arr has been declared to be
	//       an array of at least SIZE elements.
	// Post: The first SIZE elements of the array arr have been populated with
	//       random integers between 0 and one million.

	for (int i = 0; i < SIZE; i++)
	{
		arr[i] = round(rand() / double(RAND_MAX) * 1000000);
	}
}

void fill_sorted_numbers (int asc [], int dsc [], const unsigned int SIZE)
{
	// Pre:  SIZE is no more than a million, and asc and dsc have been declared
	//       to be arrays of at least SIZE elements.
	// Post: The first SIZE elements of the array asc have been populated with
	//       integers between 0 and one million in ascending order, and dsc has
	//       the same numbers in descending order.

	int step = 1000000 / SIZE;

	for (int i = 0; i < SIZE; i++)
	{
		asc[i] = dsc [SIZE-1-i] = step*i;
	}
}

/*
Example for usage:

const unsigned int SIZE = 8000;

int main()
{
	int numbers [SIZE];
	int increasing [SIZE];
	int decreasing [SIZE];

	fill_random_numbers (numbers, SIZE);
	fill_sorted_numbers (increasing, decreasing, SIZE);

	// NOW apply sorting algorithms to a copy of each array or copy of a slice of the array.

}
*/