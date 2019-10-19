/*
*版权所有 (C)2019,Yueyang.
*文件名称：common.h
*内容摘要：MATLAB工作空间顶层头文件
*当前版本：V1.0
*作者：杨越
*完成日期：20191018
*/

#ifndef COMMON_H
#define COMMON_H

#include <complex> 
#include <valarray>
#include <iostream>

#define __API__       //用户函数标记
#define __KERNEl__    //内部函数标记

using namespace std;

//浮点数精确度，用来做浮点数的等值判断
const float FLOATERROR = 1.0e-6F;
const float DOUBLEERROR =1.0e-15;

const long double LONGDOUBLEERROR = 1.0e-30;

const double GOLDNO =0.61903339;//黄金分割率

 
using namespace std;

//返回任意类型数值的绝对值
template <class T>
inline long double		//判断绝对值
__API__ Abs(const T& x)
{
	complex<long double> cld(x);
	long double ldAbs = abs(cld);
	return(ldAbs);
}

//取符号函数
template <class T>
inline T			//取x符号，+-0
__API__ Sgn(const T& x)
{
	return x < T(0) ? T(-1) : (x > T(0) ? T(1) : T(0));
}

//浮点型等值判断函数
//这里采用函数重载实现不同类型浮点数的判断
inline bool			//判断float浮点数相等
__API__ FloatEqual(float lhs, float rhs)
{
	if (Abs(lhs - rhs) < FLOATERROR)
		return true;
	else
		return false;
}

inline bool			//判断float浮点数不相等
__API__ FloatNotEqual(float lhs, float rhs)
{
	if (Abs(lhs - rhs) >= FLOATERROR)
		return true;
	else
		return false;
}

inline bool			//判断double浮点数相等
__API__ FloatEqual(double lhs, double rhs)
{
	if (Abs(lhs - rhs) < DOUBLEERROR)
		return true;
	else
		return false;
}

inline bool			//判断double浮点数不相等
__API__ FloatNotEqual(double lhs, double rhs)
{
	if (Abs(lhs - rhs) >= DOUBLEERROR)
		return true;
	else
		return false;
}

inline bool				//比较两long double浮点数相等
__API__ FloatEqual(long double lhs, long double rhs)
{
	if (Abs(lhs - rhs) < LONGDOUBLEERROR)
		return true;
	else
		return false;
}

inline bool				//比较两long double浮点数不相等
__API__ FloatNotEqual(long double lhs, long double rhs)
{
	if (Abs(lhs - rhs) >= LONGDOUBLEERROR)
		return true;
	else
		return false;
}
 
 //取最小值
template <class T>
inline T
__API__ Min(const T& x, const T& y)			//比较x与y，返回小者
{
	if(x < y)
		return x;
	else
		return y;
}

//返回最大值
template <class T>
__API__ T Max(const T& x, const T& y)		//求x与y的最大值，返回大者
{
	if(x > y)
		return x;
	else
		return y;
}

//打印数组(向量)所有元素值
template <class T>
__API__ void ValarrayPrint(const valarray<T>& vOut)
{
	size_t vsize = vOut.size();		//取数组元素的个数
	for(size_t st=0; st<vsize; st++)
	{
		cout.width(15);				//元素对齐，让每个元素占15列
		cout << vOut[st] << ' ';
	}
	cout << endl;
}

//打印某个指定数组(向量)元素值
template <class T>
__API__ void ValarrayPrint(const valarray<T>& vOut, size_t sPosition)
{
	cout.width(15);					//元素对齐，让每个元素占15列
	cout << vOut[sPosition] << endl;
}
#endif // !COMMON_H
