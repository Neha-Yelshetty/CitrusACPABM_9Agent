#ifndef GROVERSBANK_HPP
#define GROVERSBANK_HPP
#include <vector>
#include <sstream>
using namespace std;

// This class is used to keep the information of all the grovers

class Groversbank;


class Groversbank {
private:
  double profit;
  double hlbseverity;
  string behaviortype;
  string strategyParameters;
  int withhlbseverityyearcount;
  int nohlbseverityyearcount ;
  double previou_year_profit;

public:

    Groversbank();
    Groversbank(double a,double hlbseverity,string behaviortype,string strategyParameters,int withhlbseverityyearcount,int nohlbseverityyearcount,double previou_year_profit);

    double getgroversbankprofit() { return this->profit; }
    double getgroversbankhlbseverity() { return this->hlbseverity; }
    string getgroversbankbehaviortype() { return this->behaviortype; }
    string getgroversbankstrategyParameters() { return this->strategyParameters; }
    int getgroversbankwithhlbseverityyearcount() { return this->withhlbseverityyearcount; }
    int getgroversbanknohlbseverityyearcount() { return this->nohlbseverityyearcount; }
    double getgroversbankpreviouyearprofit() { return this->previou_year_profit; }
    string getgroversbankinformation() {
      stringstream ss;
        ss << profit << "~" << behaviortype << "~" << strategyParameters ;
       // cout<<thresholdseverity<<endl;
        return ss.str();
        }

    void setgroverbankparameters(double,double,string);
    void setgroverbankprofit(double);
    void setgroverbankhlbseverity(double);
    void setgroverbankbehaviortype(string);
    void setgroversbankstrategyParameters(string);
    void setgroversbankwithhlbseverityyearcount(int);
    void setgroversbanknohlbseverityyearcount(int);
    void setgroversbankpreviouyearprofit(double );

};


#endif