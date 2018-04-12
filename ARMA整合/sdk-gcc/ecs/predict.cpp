#include "predict.h"
#include "distribution.h"
#include "prediction.h"

#include <iostream>
#include <sstream>
#include <map>
#include <vector>

#define MAX_PERIOD 30
#include "Date.h"

//你要完成的功能总入口
void predict_server(char * info[MAX_INFO_NUM], char * data[MAX_DATA_NUM], int data_num, char * filename)
{

	std::stringstream ss;						// 字符串流
	std::string serial, flaName, time, target;	// 虚拟机ID，虚拟机规格，创建时间，优化资源维度名
	std::string str, strResult;					// 临时str，结果字符串
	Date date, dateFirstTrans, dateLastTrans,	// 当前日期，训练的最早日期，训练的最晚日期
		dateFirstPre, dateLastPre;				// 预测的最早日期，预测的最晚日期
	int sumCPU, sumMEM, sumHD,					// 物理服务器CPU核数，内存大小（GB），硬盘大小（GB）
		numFla, period, sumDate, indxDate,		// Flavor数，周期长度，总日期数，日期序号
		vCPU, vMEM, numPHY, numFlaPre,			// CPU核数，内存大小（MB），物理服务器数目,有效Flavor数
		interval;								// 预测开始与训练结束之间的空余天数

	std::map<std::string, int> mapFlaIndx;		// Flavor名map
	std::string arrFlaName[MAX_FLAVOR];			// Flavor名数组
	int arrFlaCPU[MAX_FLAVOR];					// 各Flavor对应CPU核数数组
	int arrFlaMEM[MAX_FLAVOR];					// 各Flavor对应内存数组
	int arrFlaPre[MAX_FLAVOR];					// 各Flavor预测结果
	int res[MAX_PHY][MAX_FLAVOR];				// 分配结果数组

	// 读取info
	ss << info[0];
	ss >> sumCPU >> sumMEM >> sumHD;
	ss << info[2];
	ss >> numFla;
	for (int i = 0; i < numFla; i++)
	{
		ss << info[i + 3];
		ss >> flaName >> vCPU >> vMEM;
		arrFlaName[i] = flaName;
		arrFlaCPU[i] = vCPU;
		arrFlaMEM[i] = vMEM / 1024;
		mapFlaIndx.insert(std::make_pair(flaName, i));
	}
	ss << info[4 + numFla];
	ss >> target;
	ss << info[6 + numFla];
	ss >> dateFirstPre >> time;
	ss << info[7 + numFla];
	ss >> dateLastPre >> time;
	ss.clear();
	period = dateLastPre - dateFirstPre;

	// 分析data中的日期
	ss << data[0];
	ss >> serial >> flaName >> dateFirstTrans >> time;
	ss << data[data_num - 1];
	ss >> serial >> flaName >> dateLastTrans >> time;
	ss.clear();
	sumDate = dateLastTrans - dateFirstTrans + 1;
	interval = dateFirstPre - dateLastTrans - 1;

	// 二维数据向量vecData
	std::vector<std::vector<double>> vecData(numFla, std::vector<double>(sumDate, 0));

	// 统计各flavor每周数目
	indxDate = 0;
	for (int i = 0; i < data_num; i++)
	{
		ss << data[i];
		ss >> serial >> flaName >> date >> time;
		if (mapFlaIndx.find(flaName) == mapFlaIndx.end())
			continue;
		indxDate = date - dateFirstTrans;
		vecData[mapFlaIndx[flaName]][indxDate]++;
	}
	ss.clear();

	// 预测 input: 二维数据向量, 预测周期长度, 预测开始与训练结束之间间隔日期, 虚拟服务器种类数, 虚拟服务器预测结果向量 
	PredictAll(vecData, period, interval, numFla, arrFlaPre);

	// 计算有效Flavor数
	numFlaPre = 0;
	for (int i = 0; i < numFla; i++)
		numFlaPre += arrFlaPre[i];

	// 分配
	for (int i = 0; i < MAX_PHY; i++)
		for (int j = 0; j < numFla; j++)
			res[i][j] = 0;
	numPHY = distribution(sumCPU, sumMEM, numFla, target, arrFlaCPU, arrFlaMEM, arrFlaPre, res);

	// 减去过小服务器
	int numLastPHY = 0;
	for (int i = 0; i < numFla; i++)
		numLastPHY += res[numPHY - 1][i];
	if (numLastPHY <= 5)
	{
		numPHY--;
		numFlaPre -= numLastPHY;
		for (int i = 0; i < numFla; i++)
			arrFlaPre[i] -= res[numPHY][i];
	}

	// 生成结果字符串
	strResult += std::to_string(numFlaPre);
	strResult += '\n';
	for (int i = 0; i < numFla; i++)
	{
		strResult += arrFlaName[i];
		strResult += ' ';
		strResult += std::to_string(arrFlaPre[i]);
		strResult += '\n';
	}
	strResult += '\n';
	strResult += std::to_string(numPHY);
	strResult += '\n';
	for (int i = 0; i < numPHY; i++)
	{
		strResult += std::to_string(i + 1);
		for (int j = 0; j < numFla; j++)
			if (res[i][j] != 0)
			{
				strResult += ' ';
				strResult += arrFlaName[j];
				strResult += ' ';
				strResult += std::to_string(res[i][j]);
			}
		strResult += '\n';
	}

	//// 需要输出的内容
	//char * result_file = (char *)"17\n\n0 8 0 20";

	// 直接调用输出文件的方法输出到指定文件中(ps请注意格式的正确性，如果有解，第一行只有一个数据；第二行为空；第三行开始才是具体的数据，数据之间用一个空格分隔开)
	write_result(strResult.c_str(), filename);

	return;
}
