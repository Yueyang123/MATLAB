/*
*版权所有 (C)2019,Yueyang.
*文件名称：sort.cpp
*内容摘要：MATLAB快速排序算法
*当前版本：V1.0
*作者：杨越
*完成日期：20191018
*/

#include "sort.h"
#include <iostream>

using namespace std;

/****************************************
* function          冒泡排序             *
* param     a[]     待排序的数组          *
* param     n       数组长度              *
* return            无                   *
* good time         O(n)                *
* avg time          O(n^2)              *
* bad time          O(n^2)              *
* space             O(1)                *
* stable            yes                 *
*****************************************/
/*主要实现步骤：
1.从0开始将每一个元素与他后面的元素进行比较，将符合条件的进行交换，
一次循环下来，最符合条件的将被放置在最前面
2.重复进行，第二符合条件的将被放在第二个位置
最后排序完成
*/
void BubbleSort(vector<int> &a)
{
    bool changed ;
    for (int i = 0; i < a.size(); i++)
    {
        changed=false;
        for (int j = 0; j < a.size(); j++)
        {
            if (a[j] < a[j+1])
                {
                changed=true;
                swap(a[j], a[j+1]);
                }
        }
    }
}
void BubbleSort(vector<float> &a)
{
    bool changed ;
    for (int i = 0; i < a.size(); i++)
    {
        changed=false;
        for (int j = 0; j < a.size(); j++)
        {
            if (a[j] < a[j+1])
            {
                swap(a[j], a[j+1]);
                changed=true;
            }
        }
    }
}

//衍生算法：鸡尾酒排序
//相当于将原数组分为两半，分别冒泡排序
void CockSort(vector<int> &arr)
{
        bool changed;
        for(int i = 0;i < arr.size()/2;i++){
            changed = false;
            //从前往后,比较相邻两个数,把大的放在后边.之前已放置成功的可以不再参与比较
            for(int j = i;j < arr.size()-1-i;j++){
 
                if(arr[j]>arr[j+1]) {
                    swap(arr[j], arr[j+1]);
                    changed =true;
                }
            }
            if(!changed){
                break;
            }
            for(int j = arr.size()-1-i; j > i; j--){
 
                if(arr[j]<arr[j-1]) {
                    swap(arr[j], arr[j-1]);
                    changed = true;
                }
            }
            if(!changed){
                break;
            }
        }
}

void CockSort(vector<float> &arr)
{
        bool changed;
        for(int i = 0;i < arr.size()/2;i++){
            changed = false;
            //从前往后,比较相邻两个数,把大的放在后边.之前已放置成功的可以不再参与比较
            for(int j = i;j < arr.size()-1-i;j++){
 
                if(arr[j]>arr[j+1]) {
                    swap(arr[j], arr[j+1]);
                    changed =true;
                }
            }
            if(!changed){
                break;
            }
            for(int j = arr.size()-1-i; j > i; j--){
 
                if(arr[j]<arr[j-1]) {
                    swap(arr[j], arr[j-1]);
                    changed = true;
                }
            }
            if(!changed){
                break;
            }
        }
}
/****************************************
* function          选择排序法           *
* param     a[]     待排序的数组         *
* param     n       数组长度             *
* return            无                  *
* good time         O(n^2)              *
* avg time          O(n^2)              *
* bad time          O(n^2)              *
* space             O(1)                *
* stable            no                  *
*****************************************/
/*
实现步骤：
1.从数组中找到最小的元素，和第一个位置的元素互换。
2.从第二个位置开始，找到最小的元素，和第二个位置的元素互换。
 直到选出array.length-1个较小元素，剩下的最大的元素自动排在最后一位。
*/

void SelectSort(vector<int> &a)
{
    int min = 0;
    for (int i = 0; i < a.size(); i++)
    {
        min = i;
        for (int j = i + 1; j < a.size(); j++)
        {
            if (a[j] < a[min])
                min = j;
        }
        if (min != i)
            swap(a[min], a[i]);
    }
}

void SelectSort(vector<float> &a)
{
    int min = 0;
    for (int i = 0; i < a.size(); i++)
    {
        min = i;
        for (int j = i + 1; j < a.size(); j++)
        {
            if (a[j] < a[min])
                min = j;
        }
        if (min != i)
            swap(a[min], a[i]);
    }
}
/****************************************
* function          插入排序法           *
* param     a[]     待排序的数组          *
* param     n       数组长度                *
* return            无                   *
* good time         O(n)                *
* avg time          O(n^2)              *
* bad time          O(n^2)              *
* space             O(1)                *
* stable            yes                 *
*****************************************/
/*实现步骤
从第二个元素开始，将当前元素插入到前面对应位置，使当前元素i
和之前元素形成有序数组。 
*/
void InsertSort(vector<int> &a)
{
    int temp = 0;
    for (int i = 1,j = 1; i < a.size(); i++)
    {
        j = i;
        temp = a[i];  //待插入的数，假设第一个已经排好序
        while (j > 0 && temp < a[j - 1])a[j] = a[j - 1],--j; //向前逐个比较，若该位置插入不了，则往后移
        a[j] = temp; //插入
    }
}

/****************************************
* function          快速排序法--分区       *
* param     a[]     待分区的数组          *
* param     left    区间左             *
* param     right   区间右             *
* return            返回新的基准keyindex  *
*****************************************/
int Partition(int a[], int left, int right)
{
    int midIndex = left;
    while (left < right)
    {       
        //从右往左扫描，若找到比基准key小的数，则与a[midIndex]交换
        while (left < right && a[right] >= a[midIndex])
            --right;
        if (left < right)
        {
            swap(a[right], a[midIndex]);
            midIndex = right;
        }   
        //从左往右扫描，若找到比基准key大的数，则与a[midIndex]交换
        while (left < right && a[left] <= a[midIndex])
            ++left;
        if (left < right)
        {
            swap(a[left], a[midIndex]);
            midIndex = left;
        }
    }
    return midIndex;
}

/****************************************
* function          快速排序法           *
* param     a[]     待分区的数组          *
* param     left    区间左             *
* param     right   区间右             *
* return            无                   *
* good time         O(N*logN)           *
* avg time          O(N*logN)           *
* bad time          O(N^2)              *
* space             O(NlogN)            *
* stable            no                  *
*****************************************/
void Quick_Sort(int a[], int left, int right)
{
    if (left < right)
    {
        //1.随机取基准值，然后交换到left那里
//      srand(GetTickCount());
//      int m = (rand() % (right - left)) + left;
        //2.取前中后的中值，然后交换到left那里
//      int m = Mid(left, (left + right / 2), right);
//      swap(a[m], a[left]);
        int midIndex = Partition(a, left, right); //获取新的基准keyindex
        Quick_Sort(a, left, midIndex - 1);  //左半部分排序
        Quick_Sort(a, midIndex + 1, right); //右半部分排序
    }
}

void QuickSort(int a[], int n)
{
    Quick_Sort(a, 0, n - 1);
}
