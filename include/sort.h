/*
*版权所有 (C)2019,Yueyang.
*文件名称：sort.h
*内容摘要：MATLAB快速排序算法
*当前版本：V1.0
*作者：杨越
*完成日期：20191018
*/

#ifndef _SORT_H
#define _SORT_H

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <vector>

using namespace std;

    void BubbleSort(vector<int> &a);
    void BubbleSort(vector<float> &a);
    void CockSort(vector<int> &arr);
    void CockSort(vector<float> &arr);

    // void QuickSort(vector<int> &a);
    // int Partition(vector<int> &a, int left, int right);
    void InsertSort(vector<int> &a);
    void SelectSort(vector<int> &a);
    void InsertSort(vector<float> &a);
    void SelectSort(vector<float> &a);

#endif 