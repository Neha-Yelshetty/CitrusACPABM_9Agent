#include "../headers/groversbank.hpp"


Groversbank::Groversbank() {
   
}



Groversbank::Groversbank(double profit,double hlbseverity,string behaviortype,string strategyParameters,
int withhlbseverityyearcount,int nohlbseverityyearcount,double previou_year_profit) {
    this->profit = profit;
    this->hlbseverity = hlbseverity;
    this->behaviortype = behaviortype;
    this->strategyParameters = strategyParameters;
    this->withhlbseverityyearcount = withhlbseverityyearcount;
    this->nohlbseverityyearcount = nohlbseverityyearcount;
    this->previou_year_profit = previou_year_profit;
}


void Groversbank::setgroverbankparameters(double profit,double hlbseverity,string behaviortype) {
        this->profit = profit;
        this->hlbseverity = hlbseverity;
        this->behaviortype = behaviortype;
}

void Groversbank::setgroverbankprofit(double profit) {
        this->profit = profit;
}


void Groversbank::setgroverbankhlbseverity(double hlbseverity) {
        this->hlbseverity = hlbseverity;
}

void Groversbank::setgroverbankbehaviortype(string behaviortype) {
        this->behaviortype = behaviortype;
}

void Groversbank::setgroversbankstrategyParameters(string strategyParameters) {
        this->strategyParameters = strategyParameters;
}

void Groversbank::setgroversbankwithhlbseverityyearcount(int withhlbseverityyearcount){
        this->withhlbseverityyearcount = withhlbseverityyearcount;
}

void Groversbank::setgroversbanknohlbseverityyearcount(int nohlbseverityyearcount){
        this->nohlbseverityyearcount = nohlbseverityyearcount;
}

void Groversbank::setgroversbankpreviouyearprofit(double previou_year_profit){
        this->previou_year_profit = previou_year_profit;
}

void Groversbank::setgroversbankcost(double cost){
        this->cost = cost;
}

void Groversbank::setgroversbankreturnvaluet(double return_value){
        this->return_value = return_value;
}

void Groversbank::setgroversbankcummulativeprofit(double cummulative_profit){
        this->cummulative_profit = cummulative_profit;
}

