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
  double cummulative5yearprofit;
  double cummulative5yearhlbzeroprofit;
  double annualprofit;
  double hlbzeroannualprofit;
  string stratergyname;
  string stratergyparameter;

    previousyearprofitdata();
    //previousyearprofitdata(int time,double cummulativeprofit,double annualprofit,string stratergyname,string stratergyparameter);

    double getPreviousyearcummulativeprofit() { return this->cummulativeprofit; }
    double getPreviousyearcummulative5yearprofit() { return this->cummulative5yearprofit; }
    double getPreviousyearcummulative5yearhlbzeroprofit() { return this->cummulative5yearhlbzeroprofit; }
    double getPreviousyearannualprofit() { return this->annualprofit; }
    double getPreviousyearhlbzeroannualprofit() { return this->hlbzeroannualprofit; }
    string getPreviousyearstratergyname() { return this->stratergyname; }
    string getPreviousyearstratergyparameter() { return this->stratergyparameter; }
    int getPreviousyeartime() { return this->time; }

    void setPreviousyearprofitdata(int,double,double,string,string,double,double,double);
    void ReadPreviousData(previousyearprofitdata[],int);

     void ReadPreviousDataCommonlibrary(previousyearprofitdata[]);



    string getpreviousprofitdata() {
        stringstream ss;
        ss << time << "," << cummulativeprofit << "," << annualprofit << "," << stratergyname << "," << stratergyparameter<<","<<cummulative5yearprofit<<","<<cummulative5yearhlbzeroprofit<<","<<hlbzeroannualprofit;
        return ss.str();
        //return stratergyparameter;
    }



};


#endif