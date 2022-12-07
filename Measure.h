#ifndef MEASURE_H
#define MEASURE_H

#include <mkl.h>

class Measure
{
public:
	explicit Measure(int d);
	~Measure();
	void measurePI(int d2);
	void measureIPI(int d1, int d2);
	void measureIP(int d1);
	MKL_Complex8 *GetValues();
	int *GetCol();
	int *GetPB();
	int *GetPE();
private:
	int d1;
	int d2;
	MKL_Complex8 *valuesPtr;
	int *colPtr;
	int *pBPtr;
	int *pEPtr;
};
#endif


