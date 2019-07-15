/* ******************************************************* */
//大连理工大学凌云超算学生工作室ASCDUT出品
//串行版程序原作者：魏晶晶
//并行版程序贡献者：叶林伟
//
//Version 1.0
//本程序说明：
//本程序使用MPI实现，但数据划分类似与“十字形切蛋糕”
//如下图，4个mpi进程在计算域中的数据划分情况
//	 ———————
//	| 2 | 3 |
//	|———|———|	
//	| 0 | 1 |
//	 ———————
//值得一提的是这种划分方式，纵向的数据重叠部分的数据交换需要用小数组倒个手
//这制约了数据的交换性能
/* ******************************************************* */
#include "stdio.h"
#include "iostream"
#include "ctime"
#include "mpi.h"
using namespace std;

int nStep = 5000;							            //迭代总次数
int nstep = 0;								            //经历的迭代次数

const int nDom_X = 200;						        //真实计算域大小
const int nDom_Y = 200;

const int holo_x = 1;						          //真实计算域外光环，辅助边界计算
const int holo_y = 1;

const int global_x = nDom_X + 2 * holo_x;	//加环计算域大小
const int global_y = nDom_Y + 2 * holo_y;

const int global_x_start = holo_x;			  //加环后真实计算域范围
const int global_x_end   = holo_x + nDom_X;
const int global_y_start = holo_y;
const int global_y_end   = holo_y + nDom_Y;

const int local_x = nDom_X / 2;				    //划分给子进程的计算域大小
const int local_y = nDom_Y / 2;

int local_x_start;							          //加环后子进程的计算域范围
int local_x_end;
int local_y_start;
int local_y_end;

int size;//size == 4						          //mpi总进程数
int rank;									                //本进程进程号
MPI_Status status;

float buffer[2][global_x][global_y];		//buffer[0]即fConten_pre存储数据的二维数组，上一时刻
											//buffer[1]即fConten_now存储数据的二维数组，当前时刻
float buf[local_x][local_y];

const float fDelta_Xh = 6e-5f;				//空间步长
const float fDelta_Yh = 6e-5f;	

float fDelta_Time;							//时间步长

int Center_X;								//计算域中心点坐标
int Center_Y;

float fDiff_Coeffient;						//扩散系数


void Intial_Cacul_Domain()
{
	//初始化加环后子进程的计算域范围
	local_x_start = (rank % 2) * local_x + holo_x;
	local_x_end   = local_x_start + local_x;
	local_y_start = (rank / 2) * local_y + holo_y;
	local_y_end   = local_y_start + local_y;

	//将计算域所有节点值初始化，低的值
	for(int k = 0; k < 2; k ++)
		for(int i = 0; i < global_x; i ++)
			for(int j = 0; j < global_y; j ++)
				buffer[k][i][j] = 0.6f;
}


void Intial_Ncleus()
{
	Center_X = global_x / 2;
	Center_Y = global_y / 2;
	buffer[0][Center_X][Center_Y] = 20.f;	//计算域中心的节点设置高的温度值
}


void Intial_Physical_Variation()//设置一个扩散系数
{
	fDiff_Coeffient = 2e-009f;
}


void Caculation()
{
	float deltaX = 0.f;
	float deltaY = 0.f;

	int now = nstep % 2;		//通过nstep的奇偶变化交换pre和now的指针
	int pre = ( nstep + 1) % 2; 
//MPI并行实现
	//绝热边界条件,B为坐标轴原点
	for(int i = local_x_start; i < local_x_end ; i++)
	{
		for(int j = local_y_start; j < local_y_end ; j++)
		{
			deltaX = (buffer[now][i + 1][j] + buffer[now][i - 1][j] - 2.f * buffer[now][i][j]) / (fDelta_Xh * fDelta_Xh);
			deltaY = (buffer[now][i][j + 1] + buffer[now][i][j - 1] - 2.f * buffer[now][i][j]) / (fDelta_Yh * fDelta_Yh);
			buffer[pre][i][j] = buffer[now][i][j] + fDelta_Time * fDiff_Coeffient * (deltaX + deltaY);
		}
	}
/*
	//交换子进程计算域之间的数据
	float up[local_x];	//向上传递
	float down[local_x];//向下传递
		//horizontal
	if (rank % 2 == 0) 	//rank == 0 or 2
	{
		MPI_Send (&buffer[pre][local_x_end - 1][local_y_start], local_y, MPI_FLOAT, rank + 1, 12, MPI_COMM_WORLD);
		MPI_Recv (&buffer[pre][local_x_end    ][local_y_start], local_y, MPI_FLOAT, rank + 1, 21, MPI_COMM_WORLD, &status);
	}
	else				//rank == 1 or 3
	{
		MPI_Send (&buffer[pre][local_x_start    ][local_y_start], local_y, MPI_FLOAT, rank - 1, 21, MPI_COMM_WORLD);
		MPI_Recv (&buffer[pre][local_x_start - 1][local_y_start], local_y, MPI_FLOAT, rank - 1, 12, MPI_COMM_WORLD, &status);		
	}
		//vertical
	if (rank / 2 == 0)	//rank == 0 or 1
	{
		for (int i = local_x_start; i < local_x_end; i ++)
			up[i] = buffer[pre][i][local_y_end - 1];
		MPI_Send (up,   local_x, MPI_FLOAT, rank + 2, 34, MPI_COMM_WORLD);
		MPI_Recv (down, local_x, MPI_FLOAT, rank + 2, 34, MPI_COMM_WORLD, &status);
		for (int i = local_x_start; i < local_x_end; i ++)
			buffer[pre][i][local_y_end] = down[i];
	}
	else				//rank == 2 or 3
	{
		for (int i = local_x_start; i < local_x_end; i ++)
			down[i] = buffer[pre][i][local_y_start];
		MPI_Send (down, local_x, MPI_FLOAT, rank - 2, 34, MPI_COMM_WORLD);
		MPI_Recv (up,   local_x, MPI_FLOAT, rank - 2, 34, MPI_COMM_WORLD, &status);
		for (int i = local_x_start; i < local_x_end; i ++)
			buffer[pre][i][local_y_start - 1] = up[i];
	}
*/
	//更新辅助边界计算用的光环
	//AB & CD
	for(int j = global_y_start; j < global_y_end ; j++)
	{
		buffer[pre][0][j] = buffer[pre][global_x_start][j];
		buffer[pre][global_x_end][j] = buffer[pre][global_x_end - 1][j];
	}

	for(int i = global_x_start; i < global_x_end ; i++)
	{
		buffer[pre][i][0] = buffer[pre][i][global_x_start];
		buffer[pre][i][global_y_end] = buffer[pre][i][global_y_end - 1];
	}
/* 
//串行代码
	//绝热边界条件,B为坐标轴原点
	for(int i = global_x_start; i < global_x_end ; i++)
	{
		for(int j = global_y_start; j < global_y_end ; j++)
		{
			deltaX = (buffer[now][i + 1][j] + buffer[now][i - 1][j] - 2.f * buffer[now][i][j]) / (fDelta_Xh * fDelta_Xh);
			deltaY = (buffer[now][i][j + 1] + buffer[now][i][j - 1] - 2.f * buffer[now][i][j]) / (fDelta_Yh * fDelta_Yh);
			buffer[pre][i][j] = buffer[now][i][j] + fDelta_Time * fDiff_Coeffient * (deltaX + deltaY);
		}
	}

	//更新辅助边界计算用的光环
	//AB & CD
	for(int j = global_y_start; j < global_y_end ; j++)
	{
		buffer[pre][0][j] = buffer[pre][global_x_start][j];
		buffer[pre][global_x_end][j] = buffer[pre][global_x_end - 1][j];
	}

	for(int i = global_x_start; i < global_x_end ; i++)
	{
		buffer[pre][i][0] = buffer[pre][i][global_x_start];
		buffer[pre][i][global_y_end] = buffer[pre][i][global_y_end - 1];
	}
*/
	
}


void Save_Caculation_Data()
{
	FILE* file = fopen("result-Soute_Diffusion.txt", "w+");
	//fprintf(file,"i, j,fConten_now\n");

// 	for(int i = 0; i < nDom_X; i++)
// 		for(int j = 0; j < nDom_Y; j++)
// 			fprintf(file, "%03d    %03d      %9.7f\n", i, j,fConten_now[i][j]);


//输出加环计算域
/*
	for(int j = 0; j < global_y; j++)
	{
		for(int i = 0; i < global_x; i++)
		{
			fprintf(file, "%9.7f ", buffer[1][i][j]);
		}
		fprintf(file, "\n");
	}
*/

//输出真实计算域
	for(int j = global_y_end - 1; j >= global_y_start; j--)
	{
		for(int i = global_x_start; i < global_x_end; i++)
		{
			fprintf(file, "%9.7f ", buffer[1][i][j]);
		}
		fprintf(file, "\n");
	}
	
	fclose(file);
}


int main(void)
{
	//MPI初始化
	MPI_Init(NULL, NULL);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	//计算域及相关参数初始化
	Intial_Cacul_Domain();
	Intial_Ncleus();
	Intial_Physical_Variation();
	//时间步长设置，为了保证计算收敛性
	fDelta_Time = (fDelta_Xh * fDelta_Xh) / (4.5f * fDiff_Coeffient);
	//开始迭代
	while (nstep < nStep)
	{
		Caculation();
		nstep++;
		if (rank == 0)
			printf("The %dth iteration accomplished\n", nstep);
	}
	//整合计算结果到rank 0
/* 
	if (rank == 0)
	{
		local_x_start = (1 % 2) * local_x + holo_x;
		local_x_end   = local_x_start + local_x;
		local_y_start = (1 / 2) * local_y + holo_y;
		local_y_end   = local_y_start + local_y;
		MPI_Recv (&buf, local_x*local_y, MPI_FLOAT, 1, 1, MPI_COMM_WORLD, &status);
		buffer[1][local_x+holo_x:2*local_x+holo_x-1][local_y+holo_y:2*local_y+holo_y-1] = buf[local_x+holo_x:2*local_x+holo_x-1][local_y+holo_y:2*local_y+holo_y-1];
		MPI_Recv (&buf, local_x*local_y, MPI_FLOAT, 2, 2, MPI_COMM_WORLD, &status);
		buffer[1][local_x+holo_x:2*local_x+holo_x-1][local_y+holo_y:2*local_y+holo_y-1] = buf[local_x+holo_x:2*local_x+holo_x-1][local_y+holo_y:2*local_y+holo_y-1];
		MPI_Recv (&buf, local_x*local_y, MPI_FLOAT, 3, 3, MPI_COMM_WORLD, &status);
		buffer[1][local_x+holo_x:2*local_x+holo_x-1][local_y+holo_y:2*local_y+holo_y-1] = buf[local_x+holo_x:2*local_x+holo_x-1][local_y+holo_y:2*local_y+holo_y-1];
	}
	else
	{
		buf[local_x_start:local_x_end-1][local_y_start:local_y_end-1] = buffer[1][local_x_start:local_x_end-1][local_y_start:local_y_end-1];
		MPI_Send (&buf, local_x*local_y, MPI_FLOAT, 0, rank, MPI_COMM_WORLD);
	}
*/	
	//raank 0输出计算结果
	if (rank == 0)
	{
		Save_Caculation_Data();
		cout << "Total time:" <<  clock() / CLOCKS_PER_SEC << "s" << endl;
	}
	//MPI结束
	MPI_Finalize();
	return 0;
}

