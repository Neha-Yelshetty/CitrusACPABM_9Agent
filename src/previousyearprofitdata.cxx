#include "../headers/previousyearprofitdata.hpp"
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
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

/*Previousyearprofitdata::Previousyearprofitdata(int time,double cummulativeprofit,double annualprofit,string stratergyname,string stratergyparameter) {
    this->time = time;
    this->cummulativeprofit = cummulativeprofit;
    this->annualprofit = annualprofit;
    this->stratergyname = stratergyname;
    this->stratergyparameter = stratergyparameter;
}*/


void previousyearprofitdata::setPreviousyearprofitdata(int time,double cummulativeprofit,double annualprofit,string stratergyname,string stratergyparameter) {
        this->time = time;
        this->cummulativeprofit = cummulativeprofit;
        this->annualprofit = annualprofit;
        this->stratergyname = stratergyname;
        this->stratergyparameter = stratergyparameter;
}

void previousyearprofitdata::ReadPreviousData(previousyearprofitdata pdata[])
{
    ifstream inFile;
    inFile.open("./YearlyProfit_NA_SP_RG_RS/YearlyProfit_NA_SP_RG_RS.csv");

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
            pdata[i].setPreviousyearprofitdata(stoi(tokens[0])/365,stod(tokens[1]),stod(tokens[4]),tokens[2],tokens[3]);
            i++;
            cout<<"Excel Data~~"<<stoi(tokens[0])/365<<"~~"<<stod(tokens[1])<<"~~"<<stod(tokens[4])<<"~~"<<tokens[2]<<"~~"<<tokens[3]<<endl;
        }
        
    }

   


    inFile.close();
}