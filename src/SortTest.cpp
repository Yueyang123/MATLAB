#include "sort.h"
#include "common.h"
#include <iostream>		
#include <complex>
#include <vector>

using namespace std;		//名字空间

template <class T>
void printsort(vector<T>& val)
{
    for(int i=0;i<=val.size()-1;i++)
    {
    cout<<val[i]<<endl;
    }
}

int SortTest(void)
{
    cout<<"*******************************"<<endl;
    cout<<"请选择排序类型："<<endl;
    cout<<"1.冒泡排序"<<endl;
    cout<<"2.选择排序"<<endl;
    cout<<"3.快速排序"<<endl;
    cout<<"4.插入排序"<<endl;
    cout<<"*******************************"<<endl;
    vector<int>sort;
    u8 choose;
    choose=getchar();
    cout<<"请输入排序数组："<<endl;
    u8 k;
    while(k<=10)
    {   
        k++;
        int temp;
        cin>>temp;
        sort.push_back(temp);
    }
    cout<<"要排序的数组为："<<endl;
    for(u8 i=0;i<sort.size();i++)
    cout<<sort[i]<<endl;
    
    switch (choose)
    {
    case '1':
        cout<<"冒泡排序："<<endl;
        BubbleSort(sort);
        break;

    case '2':
        cout<<"选择排序："<<endl;
        SelectSort(sort);
        break;

    case '3':
        cout<<"快速排序："<<endl;
        QuickSort(sort);
        break;

    case '4':
        cout<<"插入排序："<<endl;
        InsertSort(sort);
        break;            
            
    default:
        break;
    }
    
    printsort(sort);
    return 0;
}