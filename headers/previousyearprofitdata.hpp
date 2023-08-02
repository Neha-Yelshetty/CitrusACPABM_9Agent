#ifndef PREVIOUSYEARPROFITDATA_HPP
#define PREVIOUSYEARPROFITDATA_HPP
#include <vector>
#include <sstream>
using namespace std;

// This class is used to keep the information of all the grovers

class previousyearprofitdata;

class previousyearprofitdata {

public:

  int time;
  double cummulativeprofit;
  double annualprofit;
  string stratergyname;
  string stratergyparameter;

    previousyearprofitdata();
    //Previousyearprofitdata(int time,double cummulativeprofit,double annualprofit,string stratergyname,string stratergyparameter);

    double getPreviousyearcummulativeprofit() { return this->cummulativeprofit; }
    double getPreviousyearannualprofit() { return this->annualprofit; }
    string getPreviousyearstratergyname() { return this->stratergyname; }
    string getPreviousyearstratergyparameter() { return this->stratergyparameter; }
    int getPreviousyeartime() { return this->time; }

    void setPreviousyearprofitdata(int,double,double,string,string);
    void ReadPreviousData(previousyearprofitdata[]);



    string getpreviousprofitdata() {
        stringstream ss;
        ss << time << "," << cummulativeprofit << "," << annualprofit << "," << stratergyname << "," << stratergyparameter;
        return ss.str();
        //return stratergyparameter;
    }



};


#endif