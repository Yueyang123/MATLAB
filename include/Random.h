/*
*版权所有 (C)2019,Yueyang.
*文件名称：Random.h
*内容摘要：MATLAB随机数生成头文件
*功能描述：提供多种生成随机数的API
*当前版本：V1.0
*作者：杨越
*完成日期：20191018
*/


#ifndef _RANDOM_H
#define _RANDOM_H

#define RandCoef2053 2053	
#define RandCoef13849 13849
#define MODUL65536 65536

#include "common.h"		//公共部分头文件

//产生一个[0,1]区间内均匀分布伪随机数
double rand_01_One(double& seed);

//产生多个[0,1]区间内均匀分布伪随机数
inline void 
rand_01_Series(double& seed, valarray<double>& dp, const size_t stCount);

//产生任意[a,b]区间内一个均匀分布伪随机整数
inline size_t 
rand_ab_One(size_t a, size_t b, size_t& seed);

//产生任意[a,b]区间内均匀分布伪随机整数序列
inline void 
rand_ab_Series(size_t a, size_t b, size_t& seed, valarray<size_t>& sp, size_t stCount);

//产生一个任意均值与方差的正态分布随机数
inline double 
rand_NormalDistributing_One(double mu, double ro, double& seed);

//产生任意均值与方差的正态分布随机数序列
inline void 
rand_NormalDistributing_Series(double mu, double ro, double seed, valarray<double>& dp, size_t stCount);

#endif //_RANDOM_H

