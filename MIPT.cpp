#include <iostream>
#include <iomanip>
#include <mkl.h>
#include "kron.h"
#include <time.h>
#include <stdio.h>
#include <cmath>
#include "Measure.h"
#include <fstream>
#include <cstdlib>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

using namespace std;

int main(int argc, char *argv[])
{

	int Initial_Nqubit = atoi(argv[1]);
	int Final_Nqubit = atoi(argv[2]);
	int Incre_Nqubit = atoi(argv[3]);
	int Nqubit;
	int Nsample = atoi(argv[4]);
	float initalRate = atof(argv[5]);
	float finalRate = atof(argv[6]);
	float rateIncrease = atof(argv[7]);
	//cout << "Initial_Nqubit";
	//cin >> Initial_Nqubit;
	//cout << "Final_Nqubit";
	//cin >> Final_Nqubit;
	//cout << "Incre_Nqubit";
	//cin >> Incre_Nqubit;
	//cout << "Nsample";
	//cin >> Nsample;
	//cout << "initalRate";
	//cin >> initalRate;
	//cout << "finalRate";
	//cin >> finalRate;
	//cout << "rateIncrease";
	//cin >> rateIncrease;

	//Some Parameters
	double start, end, rateStart, rateEnd;
	int d2;
	int d1;
	float scalar;
	MKL_Complex8 aph = { (float)1,0 };
	MKL_Complex8 bet = { 0,0 };

	//File
	char efilename[110];
	char cfilename[110];
	sprintf(efilename, "%s%d~%d(%d)%s%d%s%.2f%s%.2f%s%.2f", "Nqubit=", Initial_Nqubit, Final_Nqubit, Incre_Nqubit, "Nsample=", Nsample, "IniRate=", initalRate
		, "FinalRate=", finalRate, "RateIncre=", rateIncrease);
	sprintf(cfilename, "%s%d~%d(%d)%s%d%s%.2f%s%.2f%s%.2f", "cri_Nqubit=", Initial_Nqubit, Final_Nqubit, Incre_Nqubit, "Nsample=", Nsample, "IniRate=", initalRate
		, "FinalRate=", finalRate, "RateIncre=", rateIncrease);

	ofstream entropy(efilename, ios::out);
	ofstream crtic(cfilename, ios::out);

	if (!entropy || !crtic)
	{
		cerr << "bd";
		exit(EXIT_FAILURE);
	}

	//Random
	srand(time(NULL));

	start = clock();
	for (Nqubit = Initial_Nqubit; Nqubit <= Final_Nqubit; Nqubit = Nqubit + Incre_Nqubit)
	{
		int EElength = ((int)((finalRate - initalRate) / rateIncrease) + 1)*Nsample*(Nqubit + 1);
		int AElength = ((int)((finalRate - initalRate) / rateIncrease) + 1)*(Nqubit + 1);
		kron kr((int)pow(2, Nqubit - 2), EElength, AElength);
		Measure mr((int)pow(2, Nqubit - 1));

		sparse_matrix_t A;
		struct matrix_descr descrA;
		descrA.type = SPARSE_MATRIX_TYPE_GENERAL;

		sparse_matrix_t B;
		struct matrix_descr descrB;
		descrB.type = SPARSE_MATRIX_TYPE_DIAGONAL;
		descrB.diag = SPARSE_DIAG_NON_UNIT;

		entropy << "Nqubit : " << Nqubit << "\n\n";
		cout << "Nqubit : " << Nqubit << "\n\n";
		for (int rate = 0; rate < ((int)((finalRate - initalRate) / rateIncrease + 1)); rate++)
		{
			rateStart = clock();
			entropy << "Rate : " << initalRate + rateIncrease * rate << "\n\n";
			cout << "Rate : " << initalRate + rateIncrease * rate << "\n\n";
			for (int sample = 0; sample < Nsample; sample++)
			{
				entropy << sample + 1 << " , " << setw(10) << 0 << ",";
				cout << Nqubit << " : ";
				cout << initalRate + rateIncrease * rate << " : ";
				cout << sample + 1 << " : " << setw(20) << 0 << ",";
				kr.ProductState(Nqubit);
				kr.GetEE()[rate*(Nsample+1)*Nqubit + sample * (Nqubit+1) + 0] = 0;
				for (int layer = 1; layer < Nqubit + 1; layer++)
				{
					if (Nqubit / 2 % 2 == 0)
					{
						//U1
						for (int j = 0; j < Nqubit / 2; j++)
						{

							if ((j + 1) == 1)
							{
								d1 = 1;
								d2 = (int)pow(2, Nqubit - 2);
								kr.Haar();
								kr.kronUI(d2);
								mkl_sparse_c_create_csr(&A, SPARSE_INDEX_BASE_ZERO, 4 * (int)pow(2, Nqubit - 2), 4 * (int)pow(2, Nqubit - 2), kr.GetPB(), kr.GetPE(), kr.GetCol(), kr.GetValues());
								mkl_sparse_c_mv(SPARSE_OPERATION_NON_TRANSPOSE, aph, A, descrA, kr.Getv(), bet, kr.Getvv());
								mkl_sparse_destroy(A);
								cblas_ccopy((int)pow(2, Nqubit), kr.Getvv(), 1, kr.Getv(), 1);
								
							}
							else if ((j + 1) == Nqubit / 2)
							{
								d2 = 1;
								d1 = (int)pow(2, Nqubit - 2);
								kr.Haar();
								kr.kronIU(d1);
								mkl_sparse_c_create_csr(&A, SPARSE_INDEX_BASE_ZERO, 4 * (int)pow(2, Nqubit - 2), 4 * (int)pow(2, Nqubit - 2), kr.GetPB(), kr.GetPE(), kr.GetCol(), kr.GetValues());
								mkl_sparse_c_mv(SPARSE_OPERATION_NON_TRANSPOSE, aph, A, descrA, kr.Getv(), bet, kr.Getvv());
								mkl_sparse_destroy(A);
								cblas_ccopy((int)pow(2, Nqubit), kr.Getvv(), 1, kr.Getv(), 1);

							}
							else
							{

								d1 = (int)pow(2, 2 * j);
								d2 = (int)pow(2, Nqubit - 2 - 2 * j);
								kr.Haar();
								kr.kronIUI(d1, d2);
								mkl_sparse_c_create_csr(&A, SPARSE_INDEX_BASE_ZERO, 4 * (int)pow(2, Nqubit - 2), 4 * (int)pow(2, Nqubit - 2), kr.GetPB(), kr.GetPE(), kr.GetCol(), kr.GetValues());
								mkl_sparse_c_mv(SPARSE_OPERATION_NON_TRANSPOSE, aph, A, descrA, kr.Getv(), bet, kr.Getvv());
								mkl_sparse_destroy(A);
								cblas_ccopy((int)pow(2, Nqubit), kr.Getvv(), 1, kr.Getv(), 1);

							}
						}
						//Measure1
						for (int j = 0; j < Nqubit; j++)
						{
							if (j + 1 == 1)
							{
								if ((float)rand() / RAND_MAX > initalRate + rateIncrease * rate)
								{

								}
								else
								{
									d1 = (int)pow(2, j);
									d2 = (int)pow(2, Nqubit - j - 1);
									mr.measurePI(d2);
									mkl_sparse_c_create_csr(&B, SPARSE_INDEX_BASE_ZERO, (int)pow(2, Nqubit), (int)pow(2, Nqubit), mr.GetPB(), mr.GetPE(), mr.GetCol(), mr.GetValues());
									mkl_sparse_c_mv(SPARSE_OPERATION_NON_TRANSPOSE, aph, B, descrB, kr.Getv(), bet, kr.Getvv());
									mkl_sparse_destroy(B);
									scalar = 1 / cblas_scnrm2((int)pow(2, Nqubit), kr.Getvv(), 1);
									cblas_csscal((int)pow(2, Nqubit), scalar, kr.Getvv(), 1);
									cblas_ccopy((int)pow(2, Nqubit), kr.Getvv(), 1, kr.Getv(), 1);
								}
							}
							else if (j + 1 == Nqubit)
							{
								if ((float)rand() / RAND_MAX > initalRate + rateIncrease * rate)
								{

								}
								else
								{
									d1 = (int)pow(2, j);
									d2 = (int)pow(2, Nqubit - j - 1);
									mr.measureIP(d1);
									mkl_sparse_c_create_csr(&B, SPARSE_INDEX_BASE_ZERO, (int)pow(2, Nqubit), (int)pow(2, Nqubit), mr.GetPB(), mr.GetPE(), mr.GetCol(), mr.GetValues());
									mkl_sparse_c_mv(SPARSE_OPERATION_NON_TRANSPOSE, aph, B, descrB, kr.Getv(), bet, kr.Getvv());
									mkl_sparse_destroy(B);
									scalar = 1 / cblas_scnrm2((int)pow(2, Nqubit), kr.Getvv(), 1);
									cblas_csscal((int)pow(2, Nqubit), scalar, kr.Getvv(), 1);
									cblas_ccopy((int)pow(2, Nqubit), kr.Getvv(), 1, kr.Getv(), 1);
								}
							}
							else
							{
								if ((float)rand() / RAND_MAX > initalRate + rateIncrease * rate)
								{

								}
								else
								{
									d1 = (int)pow(2, j);
									d2 = (int)pow(2, Nqubit - j - 1);
									mr.measureIPI(d1, d2);
									mkl_sparse_c_create_csr(&B, SPARSE_INDEX_BASE_ZERO, (int)pow(2, Nqubit), (int)pow(2, Nqubit), mr.GetPB(), mr.GetPE(), mr.GetCol(), mr.GetValues());
									mkl_sparse_c_mv(SPARSE_OPERATION_NON_TRANSPOSE, aph, B, descrB, kr.Getv(), bet, kr.Getvv());
									mkl_sparse_destroy(B);
									scalar = 1 / cblas_scnrm2((int)pow(2, Nqubit), kr.Getvv(), 1);
									cblas_csscal((int)pow(2, Nqubit), scalar, kr.Getvv(), 1);
									cblas_ccopy((int)pow(2, Nqubit), kr.Getvv(), 1, kr.Getv(), 1);
								}

							}
						}

						//U2
						for (int j = 0; j < Nqubit / 2 - 1; j++)
						{


							d1 = (int)pow(2, 2 * j + 1);
							d2 = (int)pow(2, Nqubit - 2 - 2 - 2 * j + 1);
							kr.Haar();
							kr.kronIUI(d1, d2);
							mkl_sparse_c_create_csr(&A, SPARSE_INDEX_BASE_ZERO, 4 * (int)pow(2, Nqubit - 2), 4 * (int)pow(2, Nqubit - 2), kr.GetPB(), kr.GetPE(), kr.GetCol(), kr.GetValues());
							mkl_sparse_c_mv(SPARSE_OPERATION_NON_TRANSPOSE, aph, A, descrA, kr.Getv(), bet, kr.Getvv());
							mkl_sparse_destroy(A);
							cblas_ccopy((int)pow(2, Nqubit), kr.Getvv(), 1, kr.Getv(), 1);


						}
						//Measure2
						for (int j = 1; j < Nqubit - 1; j++)
						{

							if ((float)rand() / RAND_MAX > initalRate + rateIncrease * rate)
							{

							}
							else
							{
								d1 = (int)pow(2, j);
								d2 = (int)pow(2, Nqubit - j - 1);
								mr.measureIPI(d1, d2);
								mkl_sparse_c_create_csr(&B, SPARSE_INDEX_BASE_ZERO, (int)pow(2, Nqubit), (int)pow(2, Nqubit), mr.GetPB(), mr.GetPE(), mr.GetCol(), mr.GetValues());
								mkl_sparse_c_mv(SPARSE_OPERATION_NON_TRANSPOSE, aph, B, descrB, kr.Getv(), bet, kr.Getvv());
								mkl_sparse_destroy(B);
								scalar = 1 / cblas_scnrm2((int)pow(2, Nqubit), kr.Getvv(), 1);
								cblas_csscal((int)pow(2, Nqubit), scalar, kr.Getvv(), 1);
								cblas_ccopy((int)pow(2, Nqubit), kr.Getvv(), 1, kr.Getv(), 1);
							}
						}
					}
					else
					{
						//U2
						for (int j = 0; j < Nqubit / 2 - 1; j++)
						{


							d1 = (int)pow(2, 2 * j + 1);
							d2 = (int)pow(2, Nqubit - 2 - 2 - 2 * j + 1);
							kr.Haar();
							kr.kronIUI(d1, d2);
							mkl_sparse_c_create_csr(&A, SPARSE_INDEX_BASE_ZERO, 4 * (int)pow(2, Nqubit - 2), 4 * (int)pow(2, Nqubit - 2), kr.GetPB(), kr.GetPE(), kr.GetCol(), kr.GetValues());
							mkl_sparse_c_mv(SPARSE_OPERATION_NON_TRANSPOSE, aph, A, descrA, kr.Getv(), bet, kr.Getvv());
							mkl_sparse_destroy(A);
							cblas_ccopy((int)pow(2, Nqubit), kr.Getvv(), 1, kr.Getv(), 1);


						}
						//Measure1
						for (int j = 0; j < Nqubit; j++)
						{
							if (j + 1 == 1)
							{
								if ((float)rand() / RAND_MAX > initalRate + rateIncrease * rate)
								{

								}
								else
								{
									d1 = (int)pow(2, j);
									d2 = (int)pow(2, Nqubit - j - 1);
									mr.measurePI(d2);
									mkl_sparse_c_create_csr(&B, SPARSE_INDEX_BASE_ZERO, (int)pow(2, Nqubit), (int)pow(2, Nqubit), mr.GetPB(), mr.GetPE(), mr.GetCol(), mr.GetValues());
									mkl_sparse_c_mv(SPARSE_OPERATION_NON_TRANSPOSE, aph, B, descrB, kr.Getv(), bet, kr.Getvv());
									mkl_sparse_destroy(B);
									scalar = 1 / cblas_scnrm2((int)pow(2, Nqubit), kr.Getvv(), 1);
									cblas_csscal((int)pow(2, Nqubit), scalar, kr.Getvv(), 1);
									cblas_ccopy((int)pow(2, Nqubit), kr.Getvv(), 1, kr.Getv(), 1);
								}
							}
							else if (j + 1 == Nqubit)
							{
								if ((float)rand() / RAND_MAX > initalRate + rateIncrease * rate)
								{

								}
								else
								{
									d1 = (int)pow(2, j);
									d2 = (int)pow(2, Nqubit - j - 1);
									mr.measureIP(d1);
									mkl_sparse_c_create_csr(&B, SPARSE_INDEX_BASE_ZERO, (int)pow(2, Nqubit), (int)pow(2, Nqubit), mr.GetPB(), mr.GetPE(), mr.GetCol(), mr.GetValues());
									mkl_sparse_c_mv(SPARSE_OPERATION_NON_TRANSPOSE, aph, B, descrB, kr.Getv(), bet, kr.Getvv());
									mkl_sparse_destroy(B);
									scalar = 1 / cblas_scnrm2((int)pow(2, Nqubit), kr.Getvv(), 1);
									cblas_csscal((int)pow(2, Nqubit), scalar, kr.Getvv(), 1);
									cblas_ccopy((int)pow(2, Nqubit), kr.Getvv(), 1, kr.Getv(), 1);
								}
							}
							else
							{
								if ((float)rand() / RAND_MAX > initalRate + rateIncrease * rate)
								{

								}
								else
								{
									d1 = (int)pow(2, j);
									d2 = (int)pow(2, Nqubit - j - 1);
									mr.measureIPI(d1, d2);
									mkl_sparse_c_create_csr(&B, SPARSE_INDEX_BASE_ZERO, (int)pow(2, Nqubit), (int)pow(2, Nqubit), mr.GetPB(), mr.GetPE(), mr.GetCol(), mr.GetValues());
									mkl_sparse_c_mv(SPARSE_OPERATION_NON_TRANSPOSE, aph, B, descrB, kr.Getv(), bet, kr.Getvv());
									mkl_sparse_destroy(B);
									scalar = 1 / cblas_scnrm2((int)pow(2, Nqubit), kr.Getvv(), 1);
									cblas_csscal((int)pow(2, Nqubit), scalar, kr.Getvv(), 1);
									cblas_ccopy((int)pow(2, Nqubit), kr.Getvv(), 1, kr.Getv(), 1);
								}

							}
						}
						//U1
						for (int j = 0; j < Nqubit / 2; j++)
						{

							if ((j + 1) == 1)
							{
								d1 = 1;
								d2 = (int)pow(2, Nqubit - 2);
								kr.Haar();
								kr.kronUI(d2);
								mkl_sparse_c_create_csr(&A, SPARSE_INDEX_BASE_ZERO, 4 * (int)pow(2, Nqubit - 2), 4 * (int)pow(2, Nqubit - 2), kr.GetPB(), kr.GetPE(), kr.GetCol(), kr.GetValues());
								mkl_sparse_c_mv(SPARSE_OPERATION_NON_TRANSPOSE, aph, A, descrA, kr.Getv(), bet, kr.Getvv());
								mkl_sparse_destroy(A);
								cblas_ccopy((int)pow(2, Nqubit), kr.Getvv(), 1, kr.Getv(), 1);

							}
							else if ((j + 1) == Nqubit / 2)
							{
								d2 = 1;
								d1 = (int)pow(2, Nqubit - 2);
								kr.Haar();
								kr.kronIU(d1);
								mkl_sparse_c_create_csr(&A, SPARSE_INDEX_BASE_ZERO, 4 * (int)pow(2, Nqubit - 2), 4 * (int)pow(2, Nqubit - 2), kr.GetPB(), kr.GetPE(), kr.GetCol(), kr.GetValues());
								mkl_sparse_c_mv(SPARSE_OPERATION_NON_TRANSPOSE, aph, A, descrA, kr.Getv(), bet, kr.Getvv());
								mkl_sparse_destroy(A);
								cblas_ccopy((int)pow(2, Nqubit), kr.Getvv(), 1, kr.Getv(), 1);

							}
							else
							{

								d1 = (int)pow(2, 2 * j);
								d2 = (int)pow(2, Nqubit - 2 - 2 * j);
								kr.Haar();
								kr.kronIUI(d1, d2);
								mkl_sparse_c_create_csr(&A, SPARSE_INDEX_BASE_ZERO, 4 * (int)pow(2, Nqubit - 2), 4 * (int)pow(2, Nqubit - 2), kr.GetPB(), kr.GetPE(), kr.GetCol(), kr.GetValues());
								mkl_sparse_c_mv(SPARSE_OPERATION_NON_TRANSPOSE, aph, A, descrA, kr.Getv(), bet, kr.Getvv());
								mkl_sparse_destroy(A);
								cblas_ccopy((int)pow(2, Nqubit), kr.Getvv(), 1, kr.Getv(), 1);

							}
						}
						//Measure2
						for (int j = 1; j < Nqubit - 1; j++)
						{

							if ((float)rand() / RAND_MAX > initalRate + rateIncrease * rate)
							{

							}
							else
							{
								d1 = (int)pow(2, j);
								d2 = (int)pow(2, Nqubit - j - 1);
								mr.measureIPI(d1, d2);
								mkl_sparse_c_create_csr(&B, SPARSE_INDEX_BASE_ZERO, (int)pow(2, Nqubit), (int)pow(2, Nqubit), mr.GetPB(), mr.GetPE(), mr.GetCol(), mr.GetValues());
								mkl_sparse_c_mv(SPARSE_OPERATION_NON_TRANSPOSE, aph, B, descrB, kr.Getv(), bet, kr.Getvv());
								mkl_sparse_destroy(B);
								scalar = 1 / cblas_scnrm2((int)pow(2, Nqubit), kr.Getvv(), 1);
								cblas_csscal((int)pow(2, Nqubit), scalar, kr.Getvv(), 1);
								cblas_ccopy((int)pow(2, Nqubit), kr.Getvv(), 1, kr.Getv(), 1);
							}
						}
					}
					kr.schmidtDecomposition((int)pow(2, Nqubit), (int)pow(2, Nqubit / 2));
					kr.EECal(Nqubit, rate, Nsample, sample, layer);
					cout << setw(20) << kr.GetEE()[rate*Nsample*(Nqubit + 1) + sample * (Nqubit + 1) + layer] << ",";				
				}
				cout << "\n";
				for (int i = 1; i < Nqubit + 1; i++)
				{
					entropy << setw(10) << kr.GetEE()[rate*Nsample*(Nqubit + 1) + sample * (Nqubit + 1) + i] << ",";
				}
				entropy << "\n";
			}
			kr.AECal(Nqubit, rate, Nsample);
			cout << Nsample + 1 << " : " << setw(20) << 0 << ",";
			for (int i = 1; i < Nqubit + 1; i++)
			{
				cout << setw(20) << kr.GetAE()[rate*(Nqubit + 1) + i] << ",";
			}
			cout << "\n";
			entropy << Nsample + 1 << " , " << setw(20) << 0 << ",";
			for (int i = 1; i < Nqubit + 1; i++)
			{
				entropy << setw(20) << kr.GetAE()[rate*(Nqubit + 1) + i] << ",";
			}
			rateEnd = clock();
			cout << "Time Cost : " << (rateEnd - rateStart) / CLOCKS_PER_SEC;
			cout << "\n\n";
			entropy << "\n";
			entropy << "Time Cost : " << (rateEnd - rateStart) / CLOCKS_PER_SEC;
			entropy << "\n\n";
		}
		cout << "\n\n";
		crtic << Nqubit << ",";
		for (int rate = 0; rate < ((int)((finalRate - initalRate) / rateIncrease + 1)); rate++)
		{
			crtic << setw(20) << kr.GetAE()[rate*(Nqubit + 1) + Nqubit] << ",";
		}
		crtic << "\n";
		entropy << "\n\n";
	}
	end = clock();
	cout << "Total Time Cost : " << (end - start) / CLOCKS_PER_SEC;
	entropy << "Total Time Cost : " << (end - start) / CLOCKS_PER_SEC;
}

