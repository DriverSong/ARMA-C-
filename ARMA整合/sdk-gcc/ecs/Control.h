

#ifndef _MYARMA_COMTROL_H_
#define _MYARMA_COMTROL_H_

std::vector<int> predictByAIC(std::vector<double> array, int maxP, int maxQ);
std::vector<double> preDiff(std::vector<double> array, int gap);
std::vector<double> contrastDiff(std::vector<double> array, int gap);
double gaussrand();

#endif