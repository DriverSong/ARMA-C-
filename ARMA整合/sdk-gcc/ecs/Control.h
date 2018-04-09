

#ifndef _MYARMA_COMTROL_H_
#define _MYARMA_COMTROL_H_

std::vector<int> predictByAIC(std::vector<double> array, int maxP, int maxQ, std::vector<std::vector<int>> usedPQ, int flag0);
std::vector<double> preDiff(std::vector<double> array, int gap);
std::vector<double> contrastDiff(std::vector<double> array, int gap);
double gaussrand1();

#endif