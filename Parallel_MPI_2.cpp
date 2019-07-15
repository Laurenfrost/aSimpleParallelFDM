/* ******************************************************* */
//大连理工大学凌云超算学生工作室ASCDUT出品
//串行版程序原作者：魏晶晶
//并行版程序贡献者：叶林伟
//
//Version 2.0
//本程序说明：
//2号程序依然用MPI实现，但数据划分方式不同
//区别与1的“十字形切蛋糕”，本程序横向划分数据，即每个进程分得相同行数的数据
//横向划分也是考虑到C语言数组横向存储的特点，希望可以藉此提升访存效率
//如下图，在mpi size == 4时在计算域中的数据划分情况
//	0———————————————y
//	|       0       |
//	|———————————————|
//	|       1       |
//	|———————————————|
//	|       2       |
//	|———————————————|
//	|       3       |
//	x———————————————
/* ******************************************************* */
#include "stdio.h"
#include "iostream"
#include "ctime"
#include "mpi.h"
using namespace std;

const int nStep = 5000;						//迭代总次数
const int nstep = 0;						//经历的迭代次数

const int nDom_X = 200;						//真实计算域大小
const int nDom_Y = 200;

const int holo_x = 1;						//在真实计算域外添加的辅助计算环的厚度
const int holo_y = 1;

const int global_x = nDom_X + 2 * holo_x;	//加计算辅助环后整体计算域大小
const int global_y = nDom_Y + 2 * holo_y;

const int global_x_start = holo_x;			//加环后真实计算域范围，因为C语言的数组从0开始，故-1
const int global_x_end   = holo_x + nDom_X - 1;
const int global_y_start = holo_y;
const int global_y_end   = holo_y + nDom_Y - 1;

const int lay = 1;							//分割后计算域条带之间重叠区域的厚度

const int local_x = nDom_X / 4;		 		//划分给子进程的计算域大小，因为是横向下刀“切蛋糕”，所以y不变，x均分
const int local_y = nDom_Y;

int local_x_start;							//加环后子进程的计算域范围
int local_x_end;
int local_y_start;
int local_y_end;

int size;//size == 4						//mpi总进程数
int rank;									//本进程进程号
MPI_Status status;

float buffer[2][global_x][global_y];		//buffer[0]即fConten_pre存储数据的二维数组，上一时刻
											//buffer[1]即fConten_now存储数据的二维数组，当前时刻
float buf[2][local_x][local_y];				

const float fDelta_Xh = 6e-5f;				//空间步长
const float fDelta_Yh = 6e-5f;	

float fDelta_Time;							//时间步长

int Center_X;								//计算域中心点坐标
int Center_Y;

float fDiff_Coeffient;						//扩散系数


void Intial_Cacul_Domain()
{
	//初始化加环后子进程的计算域范围
	local_y_start = global_y_start;
	local_y_end   = global_y_end;
	local_x_start = rank * local_x + lay - 1;
	local_x_end   = (rank + 1) * local_x + lay - 1;

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


void Intial_Physical_Variation()			//设置一个扩散系数
{
	fDiff_Coeffient = 2e-009f;
}


void Caculation()
{
	float deltaX = 0.f;
	float deltaY = 0.f;

	int now = nstep % 2;					//通过nstep的奇偶变化交换pre和now的指针
	int pre = ( nstep + 1) % 2; 
//MPI并行实现
	for(int i = local_x_start; i <= local_x_end ; i++)
	{
		for(int j = local_y_start; j <= local_y_end ; j++)
		{
			deltaX = (buffer[now][i + 1][j] + buffer[now][i - 1][j] - 2.f * buffer[now][i][j]) / (fDelta_Xh * fDelta_Xh);
			deltaY = (buffer[now][i][j + 1] + buffer[now][i][j - 1] - 2.f * buffer[now][i][j]) / (fDelta_Yh * fDelta_Yh);
			buffer[pre][i][j] = buffer[now][i][j] + fDelta_Time * fDiff_Coeffient * (deltaX + deltaY);
		}
	}

	//交换子进程计算域之间的数据
	if(rank != 0)	//上
	{
		MPI_Send (&buffer[pre][local_x_start  ][local_y_start], local_y, MPI_FLOAT, rank-1, rank  , MPI_COMM_WORLD);
		MPI_Recv (&buffer[pre][local_x_start-1][local_y_start], local_y, MPI_FLOAT, rank-1, rank-1, MPI_COMM_WORLD, &status);
	}
	if(rank != 3)	//下
	{
		MPI_Send (&buffer[pre][local_x_end  ][local_y_start], local_y, MPI_FLOAT, rank+1, rank  , MPI_COMM_WORLD);
		MPI_Recv (&buffer[pre][local_x_end+1][local_y_start], local_y, MPI_FLOAT, rank+1, rank+1, MPI_COMM_WORLD, &status);
	}

	//更新辅助边界计算用的光环
	for(int j = global_y_start; j <= global_y_end ; j++)
	{
		buffer[pre][0][j] = buffer[pre][global_x_start][j];
		buffer[pre][global_x_end][j] = buffer[pre][global_x_end - 1][j];
	}

	for(int i = global_x_start; i <= global_x_end ; i++)
	{
		buffer[pre][i][0] = buffer[pre][i][global_x_start];
		buffer[pre][i][global_y_end] = buffer[pre][i][global_y_end - 1];
	}

/*
//串行代码
	for(int i = global_x_start; i <= global_x_end ; i++)
	{
		for(int j = global_y_start; j <= global_y_end ; j++)
		{
			deltaX = (buffer[now][i + 1][j] + buffer[now][i - 1][j] - 2.f * buffer[now][i][j]) / (fDelta_Xh * fDelta_Xh);
			deltaY = (buffer[now][i][j + 1] + buffer[now][i][j - 1] - 2.f * buffer[now][i][j]) / (fDelta_Yh * fDelta_Yh);
			buffer[pre][i][j] = buffer[now][i][j] + fDelta_Time * fDiff_Coeffient * (deltaX + deltaY);
		}
	}

	//更新辅助边界计算用的光环
	for(int j = global_y_start; j <= global_y_end ; j++)
	{
		buffer[pre][0][j] = buffer[pre][global_x_start][j];
		buffer[pre][global_x_end][j] = buffer[pre][global_x_end - 1][j];
	}

	for(int i = global_x_start; i <= global_x_end ; i++)
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


//输出加环后的完整计算域
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

//输出完整的真实计算域
	for(int j = global_y_end; j >= global_y_start; j--)
	{
		for(int i = global_x_start; i <= global_x_end; i++)
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

	if (rank == 0)
	{
		for (int i = 1; i < 4; i ++)
			MPI_Recv (&buffer[nStep%2 + 1][local_x_start][0], local_x * (local_y + holo_y * 2), MPI_FLOAT, i, i, MPI_COMM_WORLD, &status);
	}
	else
	{
		MPI_Send (&buffer[nStep%2 + 1][local_x_start][0], local_x * (local_y + holo_y * 2), MPI_FLOAT, 0, rank, MPI_COMM_WORLD);
	}

	//rank 0输出计算结果
	if (rank == 0)
	{
		Save_Caculation_Data();
		cout << "Total time:" <<  clock() / CLOCKS_PER_SEC << "s" << endl;
	}
	//MPI结束
	MPI_Finalize();
	return 0;
}
