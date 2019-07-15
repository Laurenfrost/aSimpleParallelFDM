// Soute_Difussion.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include <iostream>
#include <vector>
#include <ctime>
using namespace std;

//const int n = 81;
//const int m = 81;
//float fConten_pre[n][m];//存储溶质浓度
//float fConten_now[n][m];

float fDomainX;		//计算域边长
float fDomainY;
int nDom_X;			//计算域X方向节点数
int nDom_Y;
float fDelta_Xh;	//空间步长
float fDelta_Yh;	

float fDelta_Time;	//时间步长

int nStep;			//迭代总次数
int nstep;			//经历的迭代次数

int Center_X;		//计算域中心点坐标
int Center_Y;

float fDiff_Coeffient;//扩散系数

vector<vector<float>> fConten_pre;	//存储数据的二维数组，上一时刻
vector<vector<float>> fConten_now;	//存储数据的二维数组，当前时刻


void Intial_Cacul_Domain()
{
	nDom_X = 300;		//计算域X方向节点数
	nDom_Y = 300;
	fDelta_Xh = 6e-5f;	//空间步长
	fDelta_Yh = 6e-5f;	
	fDomainX = fDelta_Xh * nDom_X;	//计算域尺寸  
	fDomainY = fDelta_Yh * nDom_Y;

	//初始化存储数据的二维数组
	fConten_pre.clear();
	fConten_pre.resize(nDom_X);
	fConten_now.clear();
	fConten_now.resize(nDom_X);
	for(int i = 0; i < nDom_X; i++)
	{
		fConten_pre[i].resize(nDom_Y);
		fConten_now[i].resize(nDom_Y);
	}


	//将计算域所有节点值初始化，低的值
	for(int i = 0; i < nDom_X; i++)
	{
		for(int j = 0; j < nDom_Y; j++)
		{
			fConten_pre[i][j] = fConten_now[i][j] = 0.6f;
		}
	}

	nStep = 4000;  
	nstep = 0;
}


void Intial_Ncleus()
{
	Center_X = nDom_X / 2;
	Center_Y = nDom_Y / 2;
	fConten_pre[Center_X][Center_Y] = 20.f;//计算域中心的节点设置高的温度值
}


void Intial_Physical_Variation()//设置一个扩散系数
{
	fDiff_Coeffient = 2e-009f;
}


void Caculation()
{
	//绝热边界条件,B为坐标轴原点
	for(int i = 0; i < nDom_X ; i++)
	{
		for(int j = 0; j < nDom_Y ; j++)
		{
			int nLeft = 0;
			int nRight = nDom_X - 1;
			int nBottom = 0;
			int nTop = nDom_Y - 1;
			float deltaX = 0.f;
			float deltaY = 0.f;

			if((i == nLeft) && (j == nTop)) // A
			{
				deltaX = (2.f * fConten_pre[i + 1][j] - 2.f * fConten_pre[i][j]) / (fDelta_Xh * fDelta_Xh);
				deltaY = (2.f * fConten_pre[i][j - 1] - 2.f * fConten_pre[i][j]) / (fDelta_Yh * fDelta_Yh);
				fConten_now[i][j] = fConten_pre[i][j] + fDelta_Time * fDiff_Coeffient * (deltaX + deltaY);
			}
			else if((i == nLeft) && (j == nBottom)) // B	
			{
				deltaX = (2.f * fConten_pre[i + 1][j] - 2.f * fConten_pre[i][j]) / (fDelta_Xh * fDelta_Xh);
				deltaY = (2.f * fConten_pre[i][j + 1] - 2.f * fConten_pre[i][j]) / (fDelta_Yh * fDelta_Yh);
				fConten_now[i][j] = fConten_pre[i][j] + fDelta_Time * fDiff_Coeffient * (deltaX + deltaY);
			}
			else if((i == nRight) && (j == nBottom)) // C
			{
				deltaX = (2.f * fConten_pre[i - 1][j] - 2.f * fConten_pre[i][j]) / (fDelta_Xh * fDelta_Xh);
				deltaY = (2.f * fConten_pre[i][j + 1] - 2.f * fConten_pre[i][j]) / (fDelta_Yh * fDelta_Yh);
				fConten_now[i][j] = fConten_pre[i][j] + fDelta_Time * fDiff_Coeffient * (deltaX + deltaY);
			}
			else if((i == nRight) && (j == nTop)) // D
			{
				deltaX = (2.f * fConten_pre[i - 1][j] - 2.f * fConten_pre[i][j]) / (fDelta_Xh * fDelta_Xh);
				deltaY = (2.f * fConten_pre[i][j - 1] - 2.f * fConten_pre[i][j]) / (fDelta_Yh * fDelta_Yh);
				fConten_now[i][j] = fConten_pre[i][j] + fDelta_Time * fDiff_Coeffient * (deltaX + deltaY);
			}
			else if((i == nLeft) && ((j > nBottom) && (j < nTop))) //AB	
			{
				deltaX = (2.f * fConten_pre[i + 1][j] - 2.f * fConten_pre[i][j]) / (fDelta_Xh * fDelta_Xh);
				deltaY = (fConten_pre[i][j - 1] + fConten_pre[i][j + 1] - 2.f * fConten_pre[i][j]) / (fDelta_Yh * fDelta_Yh);
				fConten_now[i][j] = fConten_pre[i][j] + fDelta_Time * fDiff_Coeffient * (deltaX + deltaY);
			}
			else if(((i > nLeft) && (i < nRight)) &&  (j == nBottom)) //BC
			{
				deltaX = (fConten_pre[i + 1][j] + fConten_pre[i - 1][j] - 2.f * fConten_pre[i][j]) / (fDelta_Xh * fDelta_Xh);
				deltaY = (2.f * fConten_pre[i][j + 1] - 2.f * fConten_pre[i][j]) / (fDelta_Yh * fDelta_Yh);
				fConten_now[i][j] = fConten_pre[i][j] + fDelta_Time * fDiff_Coeffient * (deltaX + deltaY);
			}
			else if((i == nRight) && ((j > nBottom) && (j < nTop))) // CD
			{
				deltaX = (2.f * fConten_pre[i - 1][j] - 2.f * fConten_pre[i][j]) / (fDelta_Xh * fDelta_Xh);
				deltaY = (fConten_pre[i][j + 1] + fConten_pre[i][j - 1] - 2.f * fConten_pre[i][j]) / (fDelta_Yh * fDelta_Yh);
				fConten_now[i][j] = fConten_pre[i][j] + fDelta_Time * fDiff_Coeffient * (deltaX + deltaY);
			}
			else if(((i > nLeft) && (i < nRight)) && (j == nTop)) // AD
			{
				deltaX = (fConten_pre[i - 1][j] + fConten_pre[i + 1][j] - 2.f * fConten_pre[i][j]) / (fDelta_Xh * fDelta_Xh);
				deltaY = (2.f * fConten_pre[i][j - 1] - 2.f * fConten_pre[i][j]) / (fDelta_Yh * fDelta_Yh);
				fConten_now[i][j] = fConten_pre[i][j] + fDelta_Time * fDiff_Coeffient * (deltaX + deltaY);
			}
			else //中心
			{
				deltaX = (fConten_pre[i + 1][j] + fConten_pre[i - 1][j] - 2.f * fConten_pre[i][j]) / (fDelta_Xh * fDelta_Xh);
				deltaY = (fConten_pre[i][j + 1] + fConten_pre[i][j - 1] - 2.f * fConten_pre[i][j]) / (fDelta_Yh * fDelta_Yh);
				fConten_now[i][j] = fConten_pre[i][j] + fDelta_Time * fDiff_Coeffient * (deltaX + deltaY);
			}
		}
	}

	//整个计算域经过差分计算后更新节点数据
	for(int i = 0; i < nDom_X; i++)
	{
		for(int j = 0; j < nDom_Y; j++)
		{
			fConten_pre[i][j] = fConten_now[i][j];
		}
	}	
}


void Save_Caculation_Data()
{
	FILE* file = fopen("result-Soute_Diffusion.txt", "w+");
	//fprintf(file,"i, j,fConten_now\n");

// 	for(int i = 0; i < nDom_X; i++)
// 		for(int j = 0; j < nDom_Y; j++)
// 			fprintf(file, "%03d    %03d      %9.7f\n", i, j,fConten_now[i][j]);


	for(int j = nDom_Y - 1; j >= 0; j--)
	{
		for(int i = 0; i < nDom_X; i++)
		{
			fprintf(file, "%9.7f ", fConten_now[i][j]);
		}
		fprintf(file, "\n");
	}	
	fclose(file);
}


int _tmain(int argc, _TCHAR* argv[])
{
	Intial_Cacul_Domain();
	Intial_Ncleus();
	Intial_Physical_Variation();
	//时间步长设置，为了保证计算收敛性
	fDelta_Time = (fDelta_Xh * fDelta_Xh) / (4.5f * fDiff_Coeffient);

	while (nstep < nStep)
	{
		Caculation();
		nstep++;
		printf("第%d次迭代！\n",nstep);
	}	

	Save_Caculation_Data();
	cout << "Total time:" <<  clock() / CLOCKS_PER_SEC << "s" << endl;
	return 0;
}

/*
计算判断：
1) 4个角点, A, B, C, D
2) 4条边, AB, AC, BD, CD
3) 中心区域
                 Q4
   A *************************** D
   *                             *
   *                             *
   *                             *
Q1 *          中心区域         * Q3
   *                             *
   *                             *
   *                             *
   B *************************** C  
                  Q2
*/
