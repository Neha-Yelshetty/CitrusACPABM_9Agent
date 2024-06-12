#include "../headers/groversbank.hpp"


Groversbank::Groversbank() {
   
}

Groversbank::Groversbank(double profit,double hlbseverity,string behaviortype,string strategyParameters,int withhlbseverityyearcount,int nohlbseverityyearcount) {
    this->profit = profit;
    this->hlbseverity = hlbseverity;
    this->behaviortype = behaviortype;
    this->strategyParameters = strategyParameters;
    this->withhlbseverityyearcount = withhlbseverityyearcount;
    this->nohlbseverityyearcount = nohlbseverityyearcount;
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

