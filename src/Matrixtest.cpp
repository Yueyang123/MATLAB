/*
*版权所有 (C)2019,Yueyang.
*文件名称：main.cpp
*内容摘要：MATLAB矩阵类测试函数
*当前版本：V1.0
*作者：杨越
*完成日期：20191020
*/

#include "Random.h"		//随机算法头文件
#include "common.h"
#include <iostream>			//输入输出流头文件			
#include <time.h>			//时间头文件
#include <complex>
#include "Matrix.h"

using namespace std;		//名字空间

//矩阵基本运算函数
int MATEST1(void)
{	

	//使用构造函数一定义两矩阵对象
	matrix<int> matS(2, 1);		//定义整型2*1矩阵matS
	matrix<int> matR(1, 2);		//定义整型1*2矩阵matR
	
	//数组元素赋值
	matS(0, 0) = 1;
	matS(1, 0) = 2;
	matR(0, 0) = 3;
	matR(0, 1) = 4;

	//输出两矩阵
	cout << endl << "matS : " << endl;
	
	MatrixLinePrint(matS);			//按行输出矩阵matS
	cout << endl << "matR : " << endl;
	
	MatrixLinePrint(matR);			//按行输出矩阵matR

	const double dma[4][4] = {	{  1.2,  2.6,  3.7,  4.8 },
								{  5.3,  6.0,  7.0,  8.0 },
								{  9.4, 10.0, 11.6, 12.0 },
								{  3.5, 14.0, 15.0, 16.8 }};

	const double dmb[4][4] = {	{  3.0, -3.0, -2.0,  4.0 },
								{  5.0, -5.0,  1.0,  8.0 },
								{  1.0,  8.0,  5.0, -7.0 },
								{  5.0, -1.0, -3.0, -1.0 } };

	//使用构造函数二定义两矩阵对象matA, matB
	const matrix<double> matA(&dma[0][0], 4, 4);
	const matrix<double> matB(&dmb[0][0], 4, 4);
	
	//使用构造函数三定义两矩阵对象matC
	matrix<double> matC(matA);		
	
	//输出两矩阵matA, matb
	cout << endl << "matA : " << endl;
	MatrixLinePrint(matA);			//按行输出矩阵matA

	cout << endl << "matB : " << endl;
	MatrixLinePrint(matB);			//按行输出矩阵matB

	cout << endl << "matC : " << endl;
	MatrixLinePrint(matC);			//按行输出矩阵matC

	cout << endl << "Row of matA : " << matA.GetRowNum();	//matA行数
	cout << "\t Column of matA : " << matA.GetRowNum() << endl;	//matA列数
	
	matC += matA;					//矩阵自加矩阵
	cout << endl << "matC += matA : " << endl;
	MatrixLinePrint(matC);			//按行输出矩阵matC

	matC *= matA;					//矩阵自乘矩阵
	cout << endl << "matC  *= matA : " << endl;
	MatrixLinePrint(matC);			//按行输出矩阵matC

	matC -= 12.3;					//矩阵自-数
	cout << endl << "matC -= 12.3 : " << endl;
	MatrixLinePrint(matC);			//按行输出矩阵matC

	matC = -matC;					//矩阵赋-号
	cout << endl << "matC = -matC : " << endl;
	MatrixLinePrint(matC);			//按行输出矩阵matC

	matC = matA * matB;					//矩阵乘矩阵
	cout << endl << "matC = matA * matB : " << endl;
	MatrixLinePrint(matC);			//按行输出矩阵matC

	matC = matC / 2.0;					//矩阵除以数
	cout << endl << "matC = matC / 2.0 : " << endl;
	MatrixLinePrint(matC);			//按行输出矩阵matC

	matC = 0.5 * matC;					//数乘以矩阵
	cout << endl << "matC = 0.5 * matC : " << endl;
	MatrixLinePrint(matC);			//按行输出矩阵matC

	matC = matA + matC;					//矩阵加矩阵
	cout << endl << "matC = matA + matC : " << endl;
	MatrixLinePrint(matC);			//按行输出矩阵matC

	//比较两矩阵是否相同不相同
	cout << endl << "if (matA == matC) : " << (matA == matC) << endl;	//是否matA与matC相同
	cout << endl << "if (matA != matB) : " << (matA != matB) << endl;	//是否matA与matB不相同
	
	matrix< complex<float> > clhs(2, 1);	//定义复数浮点型2*1数组clhs
	matrix< complex<float> > crhs(1, 2);	//定义复数浮点型1*2数组crhs
	
	//给复数型矩阵元素赋
	clhs(0, 0) = complex<float> (1, 2);		//给复数型矩阵元素赋值
	clhs(1, 0) = complex<float> (3, 4);
	crhs(0, 0) = complex<float> (5, 6);
	crhs(0, 1) = complex<float> (7, 8);
	
	cout << endl << "ComplexMatrix clhs :" << endl;
	MatrixLinePrint(clhs);		//输出复矩阵clhs

	cout << endl << "ComplexMatrix crhs :" << endl;
	MatrixLinePrint(crhs);		//输出复矩阵crhs

	matrix< complex<float> > cmm(clhs * crhs);	//cmm=clhs*crhs
	cout << endl << "ComplexMatrix * ComplexMatrix :" << endl;
	MatrixLinePrint(cmm);	//输出复矩阵cmm=clhs*crhs

	matrix< complex<float> > cmp(cmm);		//生成cmp，由cmm初始化
	cmp += cmm;				//cmp += cmm
	cout << endl << "ComplexMatrix * ComplexMatrix :" << endl;
	MatrixLinePrint(cmp);	//输出复矩阵cmp+=cmm
    return 0;
}

//对矩阵进行LU分解测试
int MATLUTEST(void)
{
	void DislayMatrixInfo(int iRet, matrix<double>& A, 
							matrix<double>& L, matrix<double>& U);
	int iRet;

	const double dma[4][4] = {	{  1.0,  2.0,  3.0,  4.0 },
								{  5.0,  6.0,  7.0,  8.0 },
								{  9.0, 10.0, 11.0, 12.0 },
								{ 13.0, 14.0, 15.0, 16.0 }};

	const double dmb[4][4] = {	{  3.0, -3.0, -2.0, 2.0 },
								{ -3.0,  8.0,  4.0, 1.0 },
								{ -2.0,  4.0,  6.0, 3.0 }, 
								{  2.0,  1.0,  3.0, 9.0 } };

	const double dmc[4][4] = {	{2.0, 4.0,  4.0, 2.0},
								{3.0, 3.0, 12.0, 6.0},
								{2.0, 4.0, -1.0, 2.0},
								{4.0, 2.0,  1.0, 1.0}	};

	matrix<double> matA(&dma[0][0], 4, 4);
	matrix<double> matB(&dmb[0][0], 4, 4);
	matrix<double> matC(&dmc[0][0], 4, 4);

	matrix<double> matA_L(4,4);
	matrix<double> matA_U(4,4);
	
	iRet = MatrixLU(matA, matA_L, matA_U);

	cout << "matA : " << endl;
	DislayMatrixInfo(iRet, matA, matA_L, matA_U); 

	iRet = MatrixLU(matB, matA_L, matA_U);
	cout << endl << "matB : " << endl ;

	DislayMatrixInfo(iRet, matB, matA_L, matA_U); 

	iRet = MatrixLU(matC, matA_L, matA_U);
	cout << endl << "matC : " << endl;

	DislayMatrixInfo(iRet, matC, matA_L, matA_U); 

}
void DislayMatrixInfo(int iRet, matrix<double>& A, 
						matrix<double>& L, matrix<double>& U)
{
	if(iRet==-1)
	{	
		cout << "矩阵不是方阵！" << endl;
	}
	else
		if(iRet==0)
			cout << "矩阵分解失败！" << endl;
		else
		{
			MatrixLinePrint(A);
			cout << endl << "L : " << endl;
			MatrixLinePrint(L);
			cout << endl << "U : " << endl;
			MatrixLinePrint(U);
		}
}

