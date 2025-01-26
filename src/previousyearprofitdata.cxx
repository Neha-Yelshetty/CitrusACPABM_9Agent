#include "../headers/previousyearprofitdata.hpp"
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include <iomanip>
#include <math.h>
#include <algorithm>
#include <random>
#include <boost/algorithm/string.hpp>
#include <boost/math/distributions.hpp>
#include <boost/math/distributions/normal.hpp>



using namespace std;

previousyearprofitdata::previousyearprofitdata() {  
}

/*previousyearprofitdata::previousyearprofitdata(int time,double cummulativeprofit,double annualprofit,string stratergyname,string stratergyparameter) {
    this->time = 0;
    this->cummulativeprofit = 0;
    this->annualprofit = 0;
    this->stratergyname = "";
    this->stratergyparameter = "";
}*/


void previousyearprofitdata::setPreviousyearprofitdata(int time,double cummulativeprofit,double annualprofit,string stratergyname,string stratergyparameter,double cummulative5yearprofit,double cummulative5yearhlbzeroprofit,double hlbzeroannualprofit) {
        this->time = time;
        this->cummulativeprofit = cummulativeprofit;
        this->annualprofit = annualprofit;
        this->stratergyname = stratergyname;
        this->stratergyparameter = stratergyparameter;
        this->cummulative5yearprofit = cummulative5yearprofit;
        this->cummulative5yearhlbzeroprofit = cummulative5yearhlbzeroprofit;
        this->hlbzeroannualprofit = hlbzeroannualprofit;
}

void previousyearprofitdata::ReadPreviousData(previousyearprofitdata pdata[],int currentyear)
{
    //if(pdata[0].getPreviousyeartime() != currentyear)
    //{
        ifstream inFile;
        inFile.open("./YearlyProfit_NA_SP_RG_RS/YearlyProfit_NA_SP_RG_RS_OTC.csv");

        std::vector<std::string> tokens;
        std::string token;
        string readString;
        int i = 0;
        while (!inFile.eof()) {
            readString = "";
            tokens.clear() ;
            getline(inFile,readString);
            std::stringstream tempstr; 
            
            tempstr<< readString;
            while (std::getline(tempstr, token, ',')) {
                tokens.push_back(token);
            }
            if(tokens[0] != "t")
            {
                int y = stoi(tokens[0])/365;
                pdata[i].setPreviousyearprofitdata(stoi(tokens[0])/365,stod(tokens[1]),stod(tokens[4]),tokens[2],tokens[3],stod(tokens[5]),stod(tokens[6]),stod(tokens[7]));
                i++;
            }
            
        }
        inFile.close();
   // }
}

void previousyearprofitdata::ReadPreviousDataCommonlibrary(previousyearprofitdata pdata[])
{
    string filename = "./YearlyProfit_NA_SP_RG_RS/YearlyProfit_NA_SP_RG_RS.csv";
    fstream infile(filename);
    char buffer[65536];
    infile.rdbuf()->pubsetbuf(buffer, sizeof(buffer));
    string line;
    vector<string> splittedString;
    /*while (getline(infile, line)) {
        splittedString.clear();
        size_t last = 0, pos = 0;
        while ((pos = line.find('|', last)) != std::string::npos) {
            splittedString.emplace_back(line, last, pos - last);
            last = pos + 1;
            cout<<splittedString<<endl;
        }
        if (last)
            splittedString.emplace_back(line, last);
        int a = stoi(splittedString[0]);
    }*/
}