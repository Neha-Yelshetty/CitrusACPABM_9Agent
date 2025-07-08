#include <iostream>
using namespace std;
#include "../headers/grove.hpp"
#include "../headers/commodity.hpp"
#include "../headers/behavior.hpp"
#include "../headers/coord.hpp"
#include "../headers/parameterSet.hpp"
#include "../headers/bioABM.h"
#include "../headers/groversbank.hpp"
#include "../headers/previousyearprofitdata.hpp"
#include "../headers/grovernetworkinfo.hpp"
//#include "../src/previousyearprofitdata.cxx"
#include <boost/algorithm/string.hpp>
#include <boost/math/distributions.hpp>
#include <boost/math/distributions/normal.hpp>
#include <math.h>
#include <algorithm>
#include <string>
#include <random>
#include <sstream>
#include <fstream>
#include <string>
#include <iomanip>

//#define SURVIVAL_DEBUG

#ifdef SURVIVAL_DEBUG
#define sdebug(x) cout << x
#else
#define sdebug(x) 
#endif

//Needed define-aheads
void Phase4();

/********************************************************************
* GLOBALS
*********************************************************************/
// Grove* agents;
Grove agents[ParameterSet::gridLength][ParameterSet::gridWidth];
Groversbank agentsinfo[ParameterSet::gridLength][ParameterSet::gridWidth];
Behavior* behaviorPatterns[3];
string harvestDays;
string yieldFilename;
double commodityPrice;
double commodityCost;
int commodityMaxAge;
ofstream outputFile;
string outputFilename;
string strategyParameters;
string strategyFlags; 
string agentsbankinfo;
string agencyFlags;
double memoryvalue;
double alphaplus ;
double alphaminus ;
double phiplus ;
double phiminus ;
double lambda ;
double referenceincome;
double wl;
double wh;
double typeofnetwork;
//string (*grovernetwrokinfo)[9];
std::vector<std::vector<string>> grovernetwrokinfo;
std::vector<std::vector<string>> groversprioritytoacessoxt;

int experimentID;
int memorylength = ParameterSet::parammemorylength;
double memorylengthyear[ParameterSet::parammemorylength];
Commodity cropvalue =Commodity();
previousyearprofitdata prevdata = previousyearprofitdata();
previousyearprofitdata pdata[6105];
int totalcountofpreviousdata = 6105;

Grovernetworkinfo gnw = Grovernetworkinfo();

boost::random::mt19937 econ_rng(std::time(0));
boost::random::uniform_01<> econ_gen;


/*****************************************************
 * String split
 * ***************************************************/
vector<string> split(string s, string delimiter) {
    size_t pos_start = 0, pos_end, delim_len = delimiter.length();
    string token;
    vector<string> res;

    while ((pos_end = s.find (delimiter, pos_start)) != string::npos) {
        token = s.substr (pos_start, pos_end - pos_start);
        pos_start = pos_end + delim_len;
        res.push_back (token);
    }

    res.push_back (s.substr (pos_start));
    return res;
}

/**************************************************************
 * Configure initial behaviors
 * parses behavior parameter and populates behavior vector
 * ************************************************************/
void configureInitialBehaviors(Grove *g, string behaviorString) {
     //Configure initial behaviors
    vector<string> behavioridxs;
    boost::split(behavioridxs, behaviorString, boost::is_any_of(","));
    for (int i = 0; i < behavioridxs.size(); i++) {
        int idx = stoi(behavioridxs[i]);
        Behavior* b = behaviorPatterns[idx];
        g->behaviorPatterns.push_back(b);
    }
}


/**************************************************************
* Initialize CHMA 
* Take in our grove grid and uniform crop and populate the grid
*
* STRATEGY PARAMS SPEC
* 0: Rogue cost
* 1: Rogue frequency
* 2: Rogue radius
* *************************************************************/
void InitialiseCHMA(Commodity crop) {
   // agents = new Grove[ParameterSet::gridLength*ParameterSet::gridWidth];
    vector<string> sParams = split(strategyParameters, ";");
    vector<string> sFlags = split(strategyFlags, ";");
    vector<string> agencyParams = split(agencyFlags, ";");
    vector<string> sagentsbankinfo = split(agentsbankinfo,";");
   
    int latticeRows = bioABM::getNumRows(); //75
    int latticeCols = bioABM::getRowLength(); //33
    int cropRowSize = latticeRows / ParameterSet::gridLength;  // 25
    int cropColSize = latticeCols / ParameterSet::gridWidth;  // 11
    int k =0;   
    for (int i = 0; i < ParameterSet::gridLength; i++) {

         for(int j=0; j<ParameterSet::gridWidth;j++){
        //Create grove
                /*int i_lb = (i / ParameterSet::gridLength) * cropRowSize;
                int i_ub = i_lb + cropRowSize;
                int j_lb = (i % ParameterSet::gridLength) * cropColSize;
                int j_ub = j_lb + cropColSize;*/
                int i_lb = j * cropRowSize;
                int i_ub = i_lb + cropRowSize;
                int j_lb = i * cropColSize;
                int j_ub = j_lb + cropColSize;
                //cout << i <<"~~" << j <<"~~"<< i_lb <<"~~"<< i_ub <<"~~"<< j_lb <<"~~"<< j_ub << endl;

                bool agency = (stoi(agencyParams[0]) == 1);
                 //Create and assign strategy vector
                vector<string> sFlags_agent = split(sFlags[k], ",");
                vector<string> sParams_agent = split(sParams[k], ",");
                vector<string> ssagentsbankinfo = split(sagentsbankinfo[k], ",");
                k = k+1;
              
                agents[i][j] = Grove(crop, agency, i_lb, i_ub, j_lb, j_ub);

                //Rogue trees
                if (stoi(sFlags_agent[0]) == 1) {
                    stringstream ss;
                    ss << stod(sParams_agent[2]) << ";" << stod(sParams_agent[3]) << ";" << stod(sParams_agent[0])<< ";" << stod(sParams_agent[1])<< ";" << stod(sParams_agent[4]);
                   agentsinfo[i][j] = Groversbank(stod(ssagentsbankinfo[0]),stod(ssagentsbankinfo[1]),"Rogue",ss.str(),0,0,0);
                    Behavior* rogue =  new RogueTrees(stod(sParams_agent[0]),
                                                    stod(sParams_agent[1]),
                                                    stod(sParams_agent[2]),
                                                    stod(sParams_agent[3]),
                                                    stod(sParams_agent[4]));
                                                                    
                                                        
                    agents[i][j].behaviorPatterns.push_back(rogue);
                }

                    
                //Spraying
                if (stoi(sFlags_agent[1]) == 1) {
                    stringstream ss;
                    ss << stod(sParams_agent[5]) << ";" << stod(sParams_agent[6]);
                   agentsinfo[i][j] = Groversbank(stod(ssagentsbankinfo[0]),stod(ssagentsbankinfo[1]),"Spray",ss.str(),0,0,0);
                    Behavior* spray = new SprayTrees(stod(sParams_agent[5]),
                                                    stod(sParams_agent[6]),
                                                    bioABM::getSpringStart(),
                                                    bioABM::getSummerStart(),
                                                    bioABM::getFallStart());
                    agents[i][j].behaviorPatterns.push_back(spray);
                }

                //Oxytetracycline
                if (stoi(sFlags_agent[2]) == 1) {
                    stringstream ss;
                    ss << stod(sParams_agent[7]) << ";" << stod(sParams_agent[8]);
                   agentsinfo[i][j] = Groversbank(stod(ssagentsbankinfo[0]),stod(ssagentsbankinfo[1]),"OTC",ss.str(),0,0,0);
                    Behavior* otc = new OTC(stod(sParams_agent[7]),
                                                    stod(sParams_agent[8]),
                                                    bioABM::getSpringStart(),
                                                    bioABM::getSummerStart(),
                                                    bioABM::getFallStart());
                    agents[i][j].behaviorPatterns.push_back(otc);
                }
                            
                //Denser planting
                if (stoi(sFlags_agent[3]) == 1) {
                    stringstream ss;
                    ss << stod(sParams_agent[9]) << ";" << stod(sParams_agent[10]);
                  agentsinfo[i][j] = Groversbank(stod(ssagentsbankinfo[0]),stod(ssagentsbankinfo[1]),"DensePlanting",ss.str(),0,0,0);
                    Behavior * dPlant = new DensePlanting(
                        stod(sParams_agent[7]),
                        stod(sParams_agent[8])
                    );
                    agents[i][j].behaviorPatterns.push_back(dPlant);
                }

                if (stoi(sFlags_agent[4]) == 1) {
                    stringstream ss;
                    ss << stod(sParams_agent[11]) << ";" << stod(sParams_agent[12]) << ";" << stod(sParams_agent[13])<< ";" << stod(sParams_agent[14])<< ";" << stod(sParams_agent[15])<< ";" << stod(sParams_agent[16]);
                   agentsinfo[i][j] = Groversbank(stod(ssagentsbankinfo[0]),stod(ssagentsbankinfo[1]),"RectangularRogue",ss.str(),0,0,0);
                    Behavior * wideRogue = new RectangularRogue(
                        stod(sParams_agent[11]),
                        stod(sParams_agent[12]),
                        stoi(sParams_agent[13]),
                        stoi(sParams_agent[14]),
                        stoi(sParams_agent[15]),
                        stod(sParams_agent[16])
                    );
                    agents[i][j].behaviorPatterns.push_back(wideRogue);
                }

                
                
              // cout << i << j << agentsinfo[i][j].getgroversbankprofit() << agentsinfo[i][j].getgroversbankhlbseverity() << agentsinfo[i][j].getgroversbankbehaviortype() << endl;

               // agentsinfo[i][j].setgroverbankparameters(stod(ssagentsbankinfo[0]),10,"Rogue");
        }
    }

    
      
    //Initial logging and planning
    for (int i = 0; i < ParameterSet::gridLength ; i++) { 
        for(int j=0; j<ParameterSet::gridWidth;j++){
            for (int k = 0; k < agents[i][j].behaviorPatterns.size(); k++) {
                agents[i][j].behaviorPatterns[k]->PlanActions();
            }
        }
    }

    
    
    return;
}




/*************************************************************
* Phase 1
* Psyllid growth, movement, and infection. All handled by the
* biological model
*************************************************************/
void Phase1() {
    bioABM::advanceBiologicalModel();
}

/*************************************************************
* Phase 1.1
* Change the grovers behavior type if the current practice of grover is not profitable
*************************************************************/
void ChangeGroverBehaviorType(Commodity crop,int i,int j,vector<string> sParams,string changestratergname,double changeprofit,double chnagehlbseverity) {
    int period_t = bioABM::getModelDay();

   // vector<string> sParams = split(strategyParameters, ";");
    vector<string> agencyParams = split(agencyFlags, ";");

    int latticeRows = bioABM::getNumRows();
    int latticeCols = bioABM::getRowLength();
    int cropRowSize = latticeRows / ParameterSet::gridLength;
    int cropColSize = latticeCols / ParameterSet::gridWidth;
    int k =0;   
  

   /* int i_lb = (i / ParameterSet::gridLength) * cropRowSize;
    int i_ub = i_lb + cropRowSize;
    int j_lb = (i % ParameterSet::gridLength) * cropColSize;
    int j_ub = j_lb + cropColSize;*/
    bool agency = (stoi(agencyParams[0]) == 1);


    //Create and assign strategy vector
    vector<string> sParams_agent = split(sParams[k], ";");
    k = k+1;

   // agents[i][j] = Grove(crop, agency, i_lb, i_ub, j_lb, j_ub);

   
    agents[i][j].behaviorPatterns.clear();
    //Rogue trees
    if (changestratergname == "RogueTrees" ) {
        stringstream ss;
        ss << stod(sParams_agent[2]) << ";" << stod(sParams_agent[3]) << ";" << stod(sParams_agent[0])<< ";" << stod(sParams_agent[1])<< ";" << stod(sParams_agent[4]);
    agentsinfo[i][j] = Groversbank(changeprofit,chnagehlbseverity,"Rogue",ss.str(),agentsinfo[i][j].getgroversbankwithhlbseverityyearcount(),agentsinfo[i][j].getgroversbanknohlbseverityyearcount(),agentsinfo[i][j].getgroversbankpreviouyearprofit());
       Behavior* rogue = new RogueTrees(stod(sParams_agent[2]),
                       stod(sParams_agent[3]),
                       stod(sParams_agent[0]),
                       stod(sParams_agent[1]),
                       stod(sParams_agent[4]));

        agents[i][j].behaviorPatterns.clear();
        agents[i][j].behaviorPatterns.push_back(rogue);
        agents[i][j].behaviorPatterns[0]->PlanActions();
    }
     
    //Spraying
    if (changestratergname == "SprayTrees") {
     stringstream ss;
        ss << stod(sParams_agent[0]) << ";" << stod(sParams_agent[1]);
     agentsinfo[i][j] = Groversbank(changeprofit,chnagehlbseverity,"Spray",ss.str(),agentsinfo[i][j].getgroversbankwithhlbseverityyearcount(),agentsinfo[i][j].getgroversbanknohlbseverityyearcount(),agentsinfo[i][j].getgroversbankpreviouyearprofit());
      Behavior* spray = new SprayTrees(stod(sParams_agent[0]),
                                        stod(sParams_agent[1]),
                                        bioABM::getSpringStart(),
                                        bioABM::getSummerStart(),
                                        bioABM::getFallStart());

        agents[i][j].behaviorPatterns.clear();
        agents[i][j].behaviorPatterns.push_back(spray); 
        agents[i][j].behaviorPatterns[0]->PlanActions();
        
    }

    //Oxytetracycline
    if (changestratergname == "NoActionOTC") {
     stringstream ss;
     //cout<<"Hello :" << sParams[k] << endl;
        ss << stod(sParams_agent[0]) << ";" << stod(sParams_agent[1]);
     agentsinfo[i][j] = Groversbank(changeprofit,chnagehlbseverity,"NoActionOTC",ss.str(),agentsinfo[i][j].getgroversbankwithhlbseverityyearcount(),agentsinfo[i][j].getgroversbanknohlbseverityyearcount(),agentsinfo[i][j].getgroversbankpreviouyearprofit());
      Behavior* otc = new OTC(stod(sParams_agent[0]),
                                        stod(sParams_agent[1]),
                                        bioABM::getSpringStart(),
                                        bioABM::getSummerStart(),
                                        bioABM::getFallStart());

        agents[i][j].behaviorPatterns.clear();
        agents[i][j].behaviorPatterns.push_back(otc); 
        agents[i][j].behaviorPatterns[0]->PlanActions();
        
    }

    if (changestratergname == "SprayTreesOTC") {
     stringstream ss;
      agentsinfo[i][j] = Groversbank(changeprofit,chnagehlbseverity,"SprayTreesOTC",ss.str(),agentsinfo[i][j].getgroversbankwithhlbseverityyearcount(),agentsinfo[i][j].getgroversbanknohlbseverityyearcount(),agentsinfo[i][j].getgroversbankpreviouyearprofit());
      Behavior* spray = new SprayTrees(stod(sParams_agent[0]),
                                        stod(sParams_agent[1]),
                                        bioABM::getSpringStart(),
                                        bioABM::getSummerStart(),
                                        bioABM::getFallStart());

    agents[i][j].behaviorPatterns.clear();
    agents[i][j].behaviorPatterns.push_back(spray); 
    //agents[i][j].behaviorPatterns[0]->PlanActions();
      
      Behavior* otc = new OTC(stod(sParams_agent[2]),
                                        stod(sParams_agent[3]),
                                        bioABM::getSpringStart(),
                                        bioABM::getSummerStart(),
                                        bioABM::getFallStart());

        //agents[i][j].behaviorPatterns.clear();
        agents[i][j].behaviorPatterns.push_back(otc);
        agents[i][j].behaviorPatterns[0]->PlanActions();
        
    }

    return;

}

/*************************************************************
* Phase 2
* Execution of planned actions
************************************************************/
void Phase2() {
    int period_t = bioABM::getModelDay();
    for (int i = 0; i < ParameterSet::gridLength; i++) {
        for(int j=0; j<ParameterSet::gridWidth;j++){
            if(agents[i][j].getActionType() != 1){
                int relativePeriod = period_t % 365;
                for (int k = 0; k < agents[i][j].behaviorPatterns.size(); k++) {
                    if (agents[i][j].behaviorPatterns[k]->actionPlannedOnDay(relativePeriod)) {
                        cout<< "i : " << i << "j : " << j << agents[i][j].behaviorPatterns[k]->getName();
                        agents[i][j].behaviorPatterns[k]->executeAction(&agents[i][j]);
                    
                    }
                    
                }
            }
        }
    }
    //DEPRECATED
    /*
    for (int i = 0; i < ParameterSet::gridLength; i++) {
        for (int j = 0; j < ParameterSet::gridWidth; j++) {
            int relativePeriod = period_t % 365;
            for (int k = 0; k < agents[i][j].behaviorPatterns.size(); k++) {
                plan_func action = agents[i][j].behaviorPatterns[k]->getActionAtDay(relativePeriod);
                 if (action != NULL) {
                    action(&agents[i][j]);
                 }
            }
        }
    }*/
}

/*********************************************************
 * Weibull Survival
 * Returns the survival probability given 
 * vector of characteristsics
 * *******************************************************/
double weibullSurvival(int t, string strategy, double efficacy, string groveID, double alpha) {
    //cout << "Weibull(" << t << ", " << strategy << ", " << efficacy << ", " << groveID << ", " << alpha << "): ";
    double intercept = 6.345658;
    double indv_coef = 0.025527;
    double group_coef = 0.021299;
    double e75_coef = -0.000461;
    double e85_coef = 0.004604;
    double alpha1_coef = 0.396100;
    double g01_coef = -0.135702;
    double g02_coef = 0.001482;
    double g10_coef = -0.288761;
    double g11_coef = -0.510001;
    double g12_coef = -0.270179;
    double g20_coef = -0.082489;
    double g21_coef = -0.242598;
    double g22_coef = -0.067161;
    double threeterm_coef = 0.183670;
    double e75alpha1_coef = 0.279668;
    double e85alpha1_coef = 1.003784;
    double scale = 0.549;

    double lp_alpha0 = intercept;
    bool indv = (strategy == "Individual Action");
    bool group = (strategy == "Group Action");
    bool e75 = (efficacy == 0.75);
    bool e85 = (efficacy == 0.85);
    lp_alpha0 += (indv * indv_coef) + (group * group_coef) + (e75 * e75_coef) + (e85 * e85_coef);
    
    if (groveID == "g01") {
        lp_alpha0 += g01_coef;
    }
    else if (groveID == "g02") {
        lp_alpha0 += g02_coef;
    }
    else if (groveID == "g10") {
        lp_alpha0 += g10_coef;
    }
    else if (groveID == "g11") {
        lp_alpha0 += g11_coef;
    }
    else if (groveID == "g12") {
        lp_alpha0 += g12_coef;
    }
    else if (groveID == "g20") {
        lp_alpha0 += g20_coef;
    }
    else if (groveID == "g21") {
        lp_alpha0 += g21_coef;
    }
    else if (groveID == "g22") {
        lp_alpha0 += g22_coef;
    }
    
    double lp_alpha1 = lp_alpha0;
    bool threeTerm = (e85 && (indv || group));
    lp_alpha1 += alpha1_coef + (e75 * e75alpha1_coef) + (e85 * e85alpha1_coef) + (threeTerm * threeterm_coef);

    double survival_alpha0 = exp(-1 * pow(t / exp(lp_alpha0), 1/scale));
    double survival_alpha1 = exp(-1 * pow(t / exp(lp_alpha1), 1/scale));
    double scaled = survival_alpha0 + alpha*(survival_alpha1 - survival_alpha0);
    //cout << scaled << endl;
    return scaled;
}
/**********************************************************
* Get Expected Risk
* Calculates the growers expected risk of infection based on
* their expectation, the extension agents expectation, and
* their trust in the extension agent
* TEMPORARILY DISABLED FOR MANUSCRIPT 2
***********************************************************/
double getExpectedRisk(Grove*g, int growerI, int growerJ) {
    /*string growerID = string("g") + to_string(growerI) + to_string(growerJ);
    string behaviorName;
    if (g->getBehavior() == behaviorPatterns[0]) {
        behaviorName = "No Action";
    }
    else if (g->getBehavior() == behaviorPatterns[1]) {
        behaviorName = "Individual Action";
    }
    else {
        behaviorName = "Group Action";
    }
    double extensionExpectation = 1 - weibullSurvival(bioABM::getModelDay(), behaviorName, 
                                                  g->getSprayEfficacy(), growerID,
                                                  g->getAlpha());
    double adjustedExpectation = g->getLambda() * (extensionExpectation - g->lastAdjustedRisk) + g->lastAdjustedRisk;
    g->lastGrowerRisk = adjustedExpectation;
    g->lastExtensionRisk = extensionExpectation;
    g->lastAdjustedRisk = adjustedExpectation;
    return adjustedExpectation;*/
    return 0;
    
}


double getMeanHLB(Grove g) {
    int* ibounds = g.getIBounds();
    int* jbounds = g.getJBounds();
    double totalCells = 0.0;
    double totalHLB = 0.0;
    for (int i = ibounds[0]; i < ibounds[1]; i++) {
        for (int j = jbounds[0]; j < jbounds[1]; j++) {
            if (bioABM::isTreeAlive(i,j)) {
                totalCells += 1.0;
                totalHLB += bioABM::getSeverityAt(i, j);
            }
        }
    }
    if (totalCells > 0) {
        return totalHLB / totalCells;
    } else {
        return 0;
    }
}

int getDeadTrees(Grove g) {
    int* ibounds = g.getIBounds();
    int* jbounds = g.getJBounds();
    int deadTrees = 0;
    for (int i = ibounds[0]; i < ibounds[1]; i++) {
        for (int j = jbounds[0]; j < jbounds[1]; j++) {
            if (!bioABM::isTreeAlive(i,j)) {
                deadTrees++;
            }
        }
    }
    return deadTrees;

}
/********************************************************
* Phase 3
* Behavior Determination
* TEMPORARILY DISABLED DURING MANUSCRIPT 2 PLANNING
********************************************************/
void Phase3() {
    /*
    //Only determine during a planning period
    int period_t = bioABM::getModelDay();
    if ( (period_t % ParameterSet::planningLength) != 0 ) { return; }

    

    //Determine behavior for agents with agency
    for (int i = 0; i < ParameterSet::gridLength; i++) {
        for (int j = 0; j < ParameterSet::gridWidth; j++) {
            if (!agents[i][j].hasAgency()) { continue; }
            //Gather info for assessing risk
            //Assess risk
            double risk = getExpectedRisk(&agents[i][j], i, j);
            double meanHLB = getMeanHLB(agents[i][j]);
            bool findsHLB = econ_gen(econ_rng) <= meanHLB;
            if (findsHLB || agents[i][j].foundHLB) {
                risk = 1.0;
                if (!agents[i][j].foundHLB) {
                    agents[i][j].foundHLB = true;
                    agents[i][j].foundHLB_day = bioABM::getModelDay();
                }
            }
            agents[i][j].lastAdjustedRisk = risk;
            double maxExpectedValue = 0;
            int maxEVIndex = -1;
            int totalProjectingDays = ParameterSet::projectionLength + bioABM::getModelDuration();
            for (int k = 0; k < ParameterSet::numBehaviorPatterns; k++) {
                //additional costs are per planning period, given 6 sprays
                double additionalCosts = agents[i][j].behaviorCosts[k] * 6 / 4;
                double EV = behaviorPatterns[k]->getExpectedValue(agents[i][j], risk, totalProjectingDays, 
                                                                period_t, ParameterSet::planningLength,
                                                                agents[i][j].getSprayEfficacy(), agents[i][j].getAlpha(),
                                                                additionalCosts);
                if ( (EV > maxExpectedValue) || k == 0 ) {
                    maxExpectedValue = EV;
                    maxEVIndex = k;
                }
                switch (k) {
                    case 0:
                        agents[i][j].lastNAEV = EV;
                        break;
                    case 1:
                        agents[i][j].lastISEV = EV;
                        break;
                    case 2:
                        agents[i][j].lastGSEV = EV;
                        break;
                    default:
                        break;
                }
            }
            if (agents[i][j].getBehavior() != behaviorPatterns[maxEVIndex]) {
                planBehavior(&agents[i][j], behaviorPatterns[maxEVIndex]);
                agents[i][j].setBehavior(behaviorPatterns[maxEVIndex]);
       
            }
        }
    }*/
    return;

}

/********************************************************
* Phase 4
* Plan Actions
*********************************************************/
void Phase4() {
    //only plan manually during the new year
    int period_t = bioABM::getModelDay();
    if ( (period_t % 365) != 0 ) { return; }
    for (int i = 0; i < ParameterSet::gridLength; i++) {
        for(int j=0; j<ParameterSet::gridWidth;j++){
            if(agents[i][j].getActionType() != 1){
                for (int k = 0; k < agents[i][j].behaviorPatterns.size(); k++) {
                    agents[i][j].behaviorPatterns[k]->PlanActions();
                }
            }
        }
    }
    //DEPRECATED
    /*
    for (int i = 0; i < ParameterSet::gridLength; i++) {
        for (int j = 0; j < ParameterSet::gridWidth; j++) {
            for (int k = 0; k < agents[i][j].behaviorPatterns.size(); k++) {
                agents[i][j].behaviorPatterns[k]->PlanActions(&agents[i][j]);
            }
        }
    }*/
}



/***********************************************
* Phase 5: Accounting
* Economic accounting for crops and mitigation
***********************************************/
void Phase5() {
    int t = bioABM::getModelDay();
    int rel_t = t % 365;
    double probability = 0;
    double prob_delta_pos = 0; 
    double prob_delta_neg = 0; 
    std::vector<double> current_satisfaction;
    double sum = 0;
    double mean = 0;
    double severity = 0;
    double returns = 0;
    double adjustedReturns = 0;
    double adjustedcost = 0;
    
    for (int i = 0; i < ParameterSet::gridLength; i++) {
        for(int j=0; j<ParameterSet::gridWidth;j++){
            if(agents[i][j].getActionType() != 1){
                int* ibounds = agents[i][j].getIBounds();
                int* jbounds = agents[i][j].getJBounds();
                int numCrops = (ibounds[1] - ibounds[0]) * (jbounds[1] - jbounds[0]);

                if (rel_t == 0) {
                    //VC
                    agents[i][j].costs += numCrops * agents[i][j].getCrop()->getVariableCost(); 
                    agents[i][j].costs += agents[i][j].getCrop()->costs;
                    if(!agents[i][j].behaviorPatterns.empty() && 
                        agents[i][j].behaviorPatterns.size() > 1 &&
                        agents[i][j].behaviorPatterns[1] != nullptr &&
                        agents[i][j].behaviorPatterns[1]->getName().find("OTC") != std::string::npos)
                    {
                        agents[i][j].costs += 1.94 * numCrops;
                        adjustedcost =  1.94 * numCrops;
                    }
                    //cout << i<<j<<numCrops * agents[i][j].getCrop()->getVariableCost() << "--" <<  agents[i][j].getCrop()->costs << "--"  << adjustedcost << endl;

                    //agents[i][j].costs += (agents[i][j].getCrop()->costs) * (10833);
                    //FC
                    //agents[i][j].costs += agents[i][j].getFixedCosts();
                    //agents[i][j].costs += (agents[i][j].getCrop()->costs) * (10833);
                }
                //Harvest
                if (agents[i][j].getCrop()->isHarvestPeriod(rel_t)) {
                    //Yield of crop at projected age
                    adjustedReturns =0;
                    for (int k = ibounds[0]; k < ibounds[1]; k++) {
                        for (int l = jbounds[0]; l < jbounds[1]; l++) {
                            //No yield if tree is dead
                            if (!bioABM::isTreeAlive(k,l)) {
                                continue;
                            }
                            //Projected severity based on days since initial infection
                             severity = bioABM::getSeverityAt(k, l);
                            if(!agents[i][j].behaviorPatterns.empty() &&
                                agents[i][j].behaviorPatterns.size() > 1 &&
                                agents[i][j].behaviorPatterns[1] != nullptr &&
                                agents[i][j].behaviorPatterns[1]->getName().find("OTC") != std::string::npos)
                            {                                
                                //Infected yield: Units yielded times projected decay
                                adjustedReturns +=  getInfectedYield(severity) * (1+ agents[i][j].behaviorPatterns[0]->getotcefficacy());
                            }
                            else
                            {   
                     
                                //Infected yield: Units yielded times projected decay
                                adjustedReturns +=  getInfectedYield(severity);
                                //cout<<"adjustedReturns : " << adjustedReturns << " ";   
                            }

                        }
                    }
                    returns = numCrops * agents[i][j].getCrop()->getReturns();
                    agents[i][j].returns += returns + adjustedReturns;
                    returns =0;
                }
                
                // Storing the information in the grovers bank and calculating the memory length

                double meanSeverity = getMeanHLB(agents[i][j]);

        /*
        This is used when we have multiple behaviour type
        
        stringstream strategyNames;
        for (int k = 0; k < agents[i][j].behaviorPatterns.size(); k++) {
                    strategyNames << agents[i][j].behaviorPatterns[k]->getName();
                    if (k != agents[i][j].behaviorPatterns.size() - 1) {
                        strategyNames << "-";
                    }
        }
        agentsinfo[i][j].setgroverbankparameters((agents[i][j].returns - agents[i][j].costs),meanSeverity,strategyNames.str());*/

           
            if(rel_t == 0)
            {
                int hlbyearcount = 0;
                if(meanSeverity != 0)
                {
                    hlbyearcount = agentsinfo[i][j].getgroversbankwithhlbseverityyearcount() + 1;
                    agentsinfo[i][j].setgroversbankwithhlbseverityyearcount(hlbyearcount);
                }
                else
                {
                    hlbyearcount = agentsinfo[i][j].getgroversbanknohlbseverityyearcount() + 1;
                    agentsinfo[i][j].setgroversbanknohlbseverityyearcount(hlbyearcount);
                }
                //cout<<"ishlbseverityexists:"<<rel_t<<"~~"<<meanSeverity<<"~~"<<agentsinfo[i][j].getgroversbankwithhlbseverityyearcount()<<"~~"<<agentsinfo[i][j].getgroversbanknohlbseverityyearcount()<<endl;
            }
            
             
             if(agents[i][j].behaviorPatterns.empty())
             {
                //cout<<"ifblock"<<endl;
                agentsinfo[i][j].setgroverbankparameters((agents[i][j].returns - agents[i][j].costs),meanSeverity, "NoAction");
             }
             else
             {
                //cout<<"elseblock"<<(agents[i][j].returns - agents[i][j].costs)<<meanSeverity<<endl;
                agentsinfo[i][j].setgroverbankparameters((agents[i][j].returns - agents[i][j].costs),meanSeverity, agents[i][j].behaviorPatterns[0]->getName());
             }
            //cout<<"done with the condition"<<endl;
       

       
        //DEPRECATED: Responsibility for strategy accounting moved to strategy objects
        //Mitigation costs 
        /*
        if (rel_t % ParameterSet::planningLength == 0) {
            for (int k = 0; k < agents[i].behaviorPatterns.size(); k++) {
                agents[i].costs += agents[i].behaviorPatterns[k]->getVariableCosts();
            }
        }*/
        }
        else
        {
            agents[i][j].costs = 0;
            agents[i][j].returns = 0;
            agentsinfo[i][j].setgroverbankparameters(0,0,"optout");

        }
     }
    
    }
     //cout<<"end of the phase5"<<endl;
    // DEPRECATED
    
   /* for (int i = 0; i < ParameterSet::gridLength; i++) {
        for (int j = 0; j < ParameterSet::gridWidth; j++) {
            int* ibounds = agents[i][j].getIBounds();
            int* jbounds = agents[i][j].getJBounds();
            int numCrops = (ibounds[1] - ibounds[0]) * (jbounds[1] - jbounds[0]);

            if (rel_t == 0) {
                //VC
                //agents[i][j].costs += numCrops * agents[i][j].getCrop()->getVariableCost();
                //FC
                //agents[i][j].costs += agents[i][j].getFixedCosts();
                agents[i][j].costs += agents[i][j].getCrop()->costs;
            }

            //Harvest
            if (agents[i][j].getCrop()->isHarvestPeriod(rel_t)) {
                for (int k = ibounds[0]; k < ibounds[1]; k++) {
                    for (int l = jbounds[0]; l < jbounds[1]; l++) {
                        //No yield if tree is dead
                        if (!bioABM::isTreeAlive(k,l)) {
                            continue;
                        }
                        //Projected severity based on days since initial infection
                        double severity = bioABM::getSeverityAt(k, l);
                        //Yield of crop at projected age
                        double returns = agents[i][j].getCrop()->getReturns();
                        //Infected yield: Units yielded times projected decay
                        double adjustedReturns = returns * getInfectedYield(severity);
                        agents[i][j].returns += adjustedReturns;

                    }
                }
            }

            //Mitigation costs
            if (rel_t % ParameterSet::planningLength == 0) {
                for (int k = 0; k < agents[i][j].behaviorPatterns.size(); k++) {
                    agents[i][j].costs += agents[i][j].behaviorPatterns[k]->getVariableCosts();
                }
            }
        }
    }*/
    return;
}


/*************************************************************
* Calculate Mean 
***************************************************************/
double calculatemean(std::vector<double> memprf)
{
    double summemory = std::accumulate(memprf.begin(), memprf.end(), 0.0);
    double mean = summemory / memprf.size();

    return mean;
}

/*************************************************************
* Calculate standard deviation and variance
***************************************************************/

double calculatestd(std::vector<double> memprf) {
    double variance = 0;
    double std =0;		
    double meanvalue = calculatemean(memprf);
    for(unsigned int k = 0; k < memorylength; k++)
    {    
         variance = variance + pow(memprf[k] - meanvalue, 2);
    }
    variance = variance / memprf.size();
    std = sqrt(variance);
    return std;
}

/*************************************************************
* Calculate Value Function
***************************************************************/

double calculateSatisfaction(double income, double probability, double probdeltapos, double probdeltaneg) {

        double satisfaction = 0;	
		double value = 0;                 // value function
		double probWeighting = 0;         // probability weighting function
		double probWeightingdelta = 0;
		
		double phi = 0;
        //cout <<referenceincome <<"---"<< double(income) <<"--" << probability << "--" <<probdeltapos <<"--"<<probdeltaneg << endl;
		if (income >= referenceincome) {
            //cout<<"income more"<<endl;
			value = pow(income, alphaplus);
			probability = 1-probability;
			probWeighting = ( pow(probability, phiplus) ) / pow( (pow(probability, phiplus) + pow((1 - probability), phiplus)), (1/phiplus) );
			
			probability = 1 - probdeltapos;
			probWeightingdelta = ( pow(probability, phiplus) ) / pow( (pow(probability, phiplus) + pow((1 - probability), phiplus)), (1/phiplus) );
			
			phi = probWeighting - probWeightingdelta;
            
		} 
		else {
            //cout<<"Income less"<<endl;
			value = (-1)*lambda* pow(income, alphaminus);
			
			probWeighting = ( pow(probability, phiminus) ) / pow( (pow(probability, phiminus) + pow((1 - probability), phiminus)), (1/phiminus) );
			
			probability = probdeltaneg;
			probWeightingdelta = ( pow(probability, phiminus) ) / pow( (pow(probability, phiminus) + pow((1 - probability), phiminus)), (1/phiminus) );
			
			phi = probWeighting - probWeightingdelta;
		}

		satisfaction = value*phi;

        //cout<<income<<"~"<<referenceincome<<"~"<<value<<"~"<<phi<<"~"<<satisfaction<<endl;

		return satisfaction; //satisfaction
 
}

/*************************************************************
* Calculate Satisfaction Function
***************************************************************/
void calcuatesatisfaction(int i,int j)
{
    double probability = 0;
    double prob_delta_pos = 0; 
    double prob_delta_neg = 0; 
    std::vector<double> current_satisfaction;

    double mean = calculatemean(agents[i][j].memorylengthprofit);
    double std = calculatestd(agents[i][j].memorylengthprofit);

    boost::math::normal norm(mean, std);

    for(int k = 0; k< memorylength ; k++)
    {

        probability = cdf(norm,agents[i][j].memorylengthprofit[k]); 
        prob_delta_pos = cdf(norm,agents[i][j].memorylengthprofit[k] + 1); 
        prob_delta_neg =  cdf(norm,agents[i][j].memorylengthprofit[k] - 1); 
        current_satisfaction.push_back(calculateSatisfaction(agents[i][j].memorylengthprofit[k],probability,prob_delta_pos,prob_delta_neg ));
        
    } 


    double meanSatisfaction = std::accumulate(current_satisfaction.begin(), current_satisfaction.end(), 0.0);
    cout<<"meanSatisfaction :" << meanSatisfaction<<endl;
    if(meanSatisfaction > 0)
      agents[i][j].setSatisfaction(1); //1
    else
      agents[i][j].setSatisfaction(0); 

}

/*************************************************************
* Build the memory length
***************************************************************/
void buildingmemorylength(int i, int j,int year)
{
        int count = 0;
     // Building the grovers memory length over the profit value
    if(year == 1)
    {
        while(memorylength != count)
        {
            agents[i][j].memorylengthprofit.push_back(memoryvalue);
            count++;
        }
        double profit = agents[i][j].returns - agents[i][j].costs;
        agents[i][j].memorylengthprofit.erase(agents[i][j].memorylengthprofit.begin());
        agents[i][j].memorylengthprofit.push_back(profit);
        agentsinfo[i][j].setgroverbankprofit(profit);
        agentsinfo[i][j].setgroverbankhlbseverity(getMeanHLB(agents[i][j]));
        agentsinfo[i][j].setgroversbankpreviouyearprofit(profit);

    }
    else
    {
        agents[i][j].memorylengthprofit.erase(agents[i][j].memorylengthprofit.begin());
        double profit = agents[i][j].returns - agents[i][j].costs;

        double current_profit = profit - agentsinfo[i][j].getgroversbankpreviouyearprofit();

        agentsinfo[i][j].setgroversbankpreviouyearprofit(profit);

        agents[i][j].memorylengthprofit.push_back(current_profit);
        agentsinfo[i][j].setgroverbankprofit(current_profit);
        agentsinfo[i][j].setgroverbankhlbseverity(getMeanHLB(agents[i][j]));

       //cout<< year << "~~" <<  current_profit << endl;
    }
    //End of building the memory length
            

}

/*************************************************************
* Calculate Average Population Income Change Rate

*This function updates the income growth rates based on new income for farms.
***************************************************************/
double calculateAveragePopulationIncomeChangeRate(int currentrow,int currentcol)
{
    std::vector<double> differenceIncomeYears ;
	std::vector<double> populationYearlyMeanIncome ;
	string retrunvalue = "No";

    //double currentYearAverageIncome = mean(thisYearIncome);
   // populationYearlyMeanIncome.push_back(calculatemean(thisYearIncome));
				
        std::vector<double> incomeFarmYear ;
        for(int k = 0; k < memorylength; k++) {
            for(int i=0; i< ParameterSet::gridLength;i++){
                for(int j=0; j<ParameterSet::gridWidth;j++){
                        if(agents[i][j].getActionType() != 1){
                            if(typeofnetwork == 1)
                                retrunvalue = gnw.checkbondexists(grovernetwrokinfo,1,10,currentrow,currentcol,i,j);
                            else if(typeofnetwork == 3)
                                retrunvalue = gnw.checkbondexiststypetree(grovernetwrokinfo,9,10,currentrow,currentcol,i,j);
                            else
                                retrunvalue = gnw.checkbondexists(grovernetwrokinfo,9,10,currentrow,currentcol,i,j);

                            if(retrunvalue == "Yes")
                                incomeFarmYear.push_back(agents[i][j].memorylengthprofit[k]); 
                        }
                    }
            }
            //cout<<calculatemean(incomeFarmYear)<<endl;
            populationYearlyMeanIncome.push_back(calculatemean(incomeFarmYear));
        }
        
		
		//int meanYearCount = 2;                                                 // number of years to count in the past for the mean calculation
		
		for(int i = memorylength-1; i > 0; i-- ) {
			double diff = (populationYearlyMeanIncome[i-1] -  populationYearlyMeanIncome[i]) /  populationYearlyMeanIncome[i];
           
			differenceIncomeYears.push_back( diff );   
		}

		double changeRate = calculatemean(differenceIncomeYears);

		return changeRate;
}

/*************************************************************
* Calculate Average Personal Income Change Rate

* Update personal income change rate over the previous time periods and include the current year. 
* Income is an array of [year1_income, year2_income, year3_income, ... yearN_income]
* change rate is calculated for all years as (Year_n-1 - Year_n) / year_n
* We then set the mean of all the rates of change including the most recent year. 

***************************************************************/
double calculateAveragePersonalIncomeChangeRate(int i,int j)
{
    std::vector<double> differenceIncomeYears;
		double historicalIncomeChangeRate = 0;								   
		
		
		for(int k = memorylength-1; k > 0; k-- ) {
			double diff = (agents[i][j].memorylengthprofit[k-1] -  agents[i][j].memorylengthprofit[k]) /  agents[i][j].memorylengthprofit[k];
			differenceIncomeYears.push_back( diff );   
		}
		
		historicalIncomeChangeRate = calculatemean(differenceIncomeYears);
		
		return historicalIncomeChangeRate;
}

/*************************************************************
* Calculate Income Dissimilarity
***************************************************************/
void calculateIncomeDissimilarity(int i, int j,int year)
{
    double AveragePopulationIncomeChangeRate = calculateAveragePopulationIncomeChangeRate(i,j);
    double diff = AveragePopulationIncomeChangeRate - calculateAveragePersonalIncomeChangeRate(i,j);
    //cout<< AveragePopulationIncomeChangeRate<< "__"<<calculateAveragePersonalIncomeChangeRate(i,j) << endl;
    agents[i][j].setIncomeDissimilarity(diff);
}

/*************************************************************
* Optimization
***************************************************************/
void Optimization(int i, int j,int year,double currentprofit,double meanhlbseverity)
{
    double maxamt = 0;
    int highestindex = -1;
    int k = 0;
    for (k = 0; k < totalcountofpreviousdata; k++) {
       /* if (meanhlbseverity == 0) {
            if (pdata[k].getPreviousyeartime() == agentsinfo[i][j].getgroversbanknohlbseverityyearcount() &&
                pdata[k].getPreviousyearcummulative5yearhlbzeroprofit() > currentprofit) {
                if (pdata[k].getPreviousyearcummulative5yearhlbzeroprofit() > maxamt &&
                    pdata[k].getPreviousyearstratergyname().find("OTC") == std::string::npos) { // Select non-OTC values
                    maxamt = pdata[k].getPreviousyearcummulative5yearhlbzeroprofit();
                    highestindex = k;
                }
            }
            if (pdata[k].getPreviousyeartime() > agentsinfo[i][j].getgroversbanknohlbseverityyearcount())
                break;
        } else {*/
            if (pdata[k].getPreviousyeartime() == agentsinfo[i][j].getgroversbankwithhlbseverityyearcount() &&
                pdata[k].getPreviousyearcummulative5yearprofit() > currentprofit) {
                if (pdata[k].getPreviousyearcummulative5yearprofit() > maxamt) {
                    if (meanhlbseverity > 0.003) { // Check if meanhlbseverity > 0.3%
                        if (pdata[k].getPreviousyearstratergyname().find("OTC") != std::string::npos) { // Select OTC values
                            maxamt = pdata[k].getPreviousyearcummulative5yearprofit();
                            highestindex = k;
                        }
                    } else { // Otherwise, pick non-OTC values

                        maxamt = pdata[k].getPreviousyearcummulative5yearprofit();
                        highestindex = k;
                        
                    }
                }
            }
            if (pdata[k].getPreviousyeartime() > agentsinfo[i][j].getgroversbankwithhlbseverityyearcount())
                break;
       // }
    }
     if(highestindex != -1)
     {
        // Stage 1.1 : Change of the grovers behaviour type
        //double hlbmean = getMeanHLB(agents[i][j]);
        vector<string> ssparam = split(pdata[highestindex].getPreviousyearstratergyparameter(), ",");
        vector<string> sparamstringVector;
        for (string s: ssparam) {
            //cout<<s<<":"<<endl;
            sparamstringVector.push_back(s);
        }
        double anualprofit = pdata[highestindex].getPreviousyearcummulative5yearprofit();
        //cout<<"Highest Index"<< i<<"--"<<j<< "~~"<<meanhlbseverity<<"~~" << pdata[highestindex].getPreviousyearstratergyparameter()  << "~~" << pdata[highestindex].getPreviousyearstratergyname() <<"~~"<<anualprofit<<"~~"<<currentprofit << endl;
        ChangeGroverBehaviorType(cropvalue,i,j,sparamstringVector,pdata[highestindex].getPreviousyearstratergyname(),currentprofit,meanhlbseverity); 
        
        
     }
        
}


/*************************************************************
* Imitation
***************************************************************/
void Imitation(int i, int j,int year,double currentprofit,double meanhlbseverity )
{
    double maxamt = -9999999999;
    int highestindex1 = -1,highestindex2 = -1;
    int k1,k2 = 0;
    string retrunvalue = "No";
    for ( k1 = 0; k1 < ParameterSet::gridLength; k1++) {
        for( k2=0; k2<ParameterSet::gridWidth;k2++){

            if(typeofnetwork == 1)
                retrunvalue = gnw.checkbondexists(grovernetwrokinfo,1,10,i,j,k1,k2);
            else if(typeofnetwork == 3)
                retrunvalue = gnw.checkbondexiststypetree(grovernetwrokinfo,9,10,i,j,k1,k2);
            else
                retrunvalue = gnw.checkbondexists(grovernetwrokinfo,9,10,i,j,k1,k2);
            
            //cout<<"Return value~~"<<i<<"-"<<j<<"-"<<k1<<"-"<<k2<<"-"<<retrunvalue<<endl ;
            //cout<<agentsinfo[k1][k2].getgroversbankprofit() << "~~~~" << currentprofit << endl;
            if(retrunvalue == "Yes")
            {
                if(agentsinfo[k1][k2].getgroversbankprofit() > currentprofit )
                {
                    if(agentsinfo[k1][k2].getgroversbankprofit() > maxamt )
                    {
                        maxamt = agentsinfo[k1][k2].getgroversbankprofit();
                        highestindex1 =k1;
                        highestindex2 =k2;
                    
                    }
                    // cout<<pdata[k].getpreviousprofitdata()<<"----"<<pdata[k].getPreviousyearannualprofit()<<"---"<<currentprofit<<endl;
                }
            }
       }   
    }
     if(highestindex1 != -1 && highestindex2 != -1 )
     {
        // Stage 1.1 : Change of the grovers behaviour type
        //double hlbmean = getMeanHLB(agents[highestindex1][highestindex2]);

        stringstream strategyNames;
        stringstream strategyParams;
        if (agents[highestindex1][highestindex2].behaviorPatterns.empty()) {
            strategyNames << "NoAction";
            strategyParams << "NA";
        }

        for (int k = 0; k < agents[highestindex1][highestindex2].behaviorPatterns.size(); k++) {
            strategyNames << agents[highestindex1][highestindex2].behaviorPatterns[k]->getName();
            strategyParams << agents[highestindex1][highestindex2].behaviorPatterns[k]->getParams();
            if (k != agents[highestindex1][highestindex2].behaviorPatterns.size() - 1) {
                strategyNames << "-";
                strategyParams << "-";
            }
        }

        vector<string> ssparam = split(strategyParams.str(), ",");
        vector<string> sparamstringVector;
        for (string s: ssparam) {
            sparamstringVector.push_back(s);
        }
        //double anualprofit = pdata[highestindex].getPreviousyearannualprofit();
        //cout<<"Highest Index"<< highestindex1<<"--"<<highestindex2<< "~~" << strategyParams.str() << "~~" << strategyNames.str() <<"~~"<<maxamt<<"~~"<<currentprofit << endl;
        ChangeGroverBehaviorType(cropvalue,i,j,sparamstringVector,strategyNames.str(),currentprofit,meanhlbseverity);  
     }
        
}


/***********************************************
* Phase 6: Planning for the future
* Grovers Prediciton of Surrouding
***********************************************/
void Phase6() {

      double AveragePopulationIncomeChangeRate ;
       int year = bioABM::getModelDay()/365;
       int i,j;

        for ( i = 0; i < ParameterSet::gridLength; i++) {
            for( j=0; j<ParameterSet::gridWidth;j++){
                if(agents[i][j].getActionType() != 1){
                    buildingmemorylength(i,j,year);
                }
            }
        }

       //cout<<"after memory build"<<endl;

        for ( i = 0; i < ParameterSet::gridLength; i++) {
            for( j=0; j<ParameterSet::gridWidth;j++){
                if(agents[i][j].getActionType() != 1){
                calcuatesatisfaction(i,j);
                calculateIncomeDissimilarity(i,j,year);
                if(agents[i][j].getIncomeDissimilarity() > wh)
                {
                    if(agents[i][j].getSatisfaction() == 0)
                    {
                      cout<<"Firstblock-OptOut";
                       agents[i][j].setActionType(1); // Setting for the opt-out
                    }
                    else
                    {
                      cout<<"Firstblock-Repetition";
                       agents[i][j].setActionType(2); //Repetition
                    }
                }
                else if(wl<= agents[i][j].getIncomeDissimilarity() && agents[i][j].getIncomeDissimilarity()<= wh)
                {
                    
                    if(agents[i][j].getSatisfaction() == 0) 
                    {
                       cout<<"SecondBlock-optimization";
                       agents[i][j].setActionType(3); // optimization
                    }
                    else
                    {
                      cout<<"SecondBlock-Repetition";
                       agents[i][j].setActionType(2); // Repetition
                    }
                }
                else if(agents[i][j].getIncomeDissimilarity() < wl)
                {
                    
                    if(agents[i][j].getSatisfaction() == 0)
                    {
                       cout<<"Lastblock-optout";
                       agents[i][j].setActionType(1); // Setting for the opt-out
                    }
                    else
                    {
                        cout<<"Lastblock-Imitate";
                       agents[i][j].setActionType(4); //Imitate
                    }
                }

            }
            stringstream strategyNames;
            stringstream strategyParams;
            if (agents[i][j].behaviorPatterns.empty()) {
                strategyNames << "NoAction";
                strategyParams << "NA";
            }

            for (int k = 0; k < agents[i][j].behaviorPatterns.size(); k++) {
                strategyNames << agents[i][j].behaviorPatterns[k]->getName();
                strategyParams << agents[i][j].behaviorPatterns[k]->getParams();
                if (k != agents[i][j].behaviorPatterns.size() - 1) {
                    strategyNames << "-";
                    strategyParams << "-";
                }
            }

            cout<< i << "~~" << j <<"~~"<<wl<<"~~"<<wh<<"~~" << agents[i][j].getIncomeDissimilarity()<<"~~"<<agents[i][j].getSatisfaction()<<"~~"<<strategyNames.str()<<"~~"<< strategyParams.str()<<"~~"<<agentsinfo[i][j].getgroversbankprofit()<<endl;
           }

        }

     

}

void Phase7()
{
    int i,j;
    int year = bioABM::getModelDay()/365;
    
     for ( i = 0; i < ParameterSet::gridLength; i++) {
            for( j=0; j<ParameterSet::gridWidth;j++){
                
                
                if(agents[i][j].getActionType() == 3){
                    Optimization(i,j,year,agentsinfo[i][j].getgroversbankprofit(),agentsinfo[i][j].getgroversbankhlbseverity());
                }
                else if(agents[i][j].getActionType() == 4){
                    Imitation(i,j,year,agentsinfo[i][j].getgroversbankprofit(),agentsinfo[i][j].getgroversbankhlbseverity());
                }
            }
     }
     
}


/*************************************************************
* Write CSV Line
***************************************************************/
void writeCSVLine() {
    for (int i = 0; i < ParameterSet::gridLength ; i++) {
        for(int j=0; j<ParameterSet::gridWidth;j++){
            stringstream pcummulative5yearprofitdata;
            stringstream pcummulativehlbzero5yearprofitdata;
            double meanSeverity = getMeanHLB(agents[i][j]);
            if(bioABM::getModelDay()%365 == 0)
            {
                int k;
                for(k = 0 ;k < totalcountofpreviousdata;k++ )
                {

                    if(meanSeverity == 0)
                    {
                        if(pdata[k].getPreviousyeartime() == agentsinfo[i][j].getgroversbanknohlbseverityyearcount())
                        {
                            pcummulative5yearprofitdata <<pdata[k].getPreviousyearcummulative5yearprofit()<<"~"<<pdata[k].getPreviousyearstratergyname()<<"~"<<pdata[k].getPreviousyearstratergyparameter() <<"##";
                            pcummulativehlbzero5yearprofitdata <<pdata[k].getPreviousyearcummulative5yearhlbzeroprofit()<<"~"<<pdata[k].getPreviousyearstratergyname()<<"~"<<pdata[k].getPreviousyearstratergyparameter() <<"##";
                        }
                    }else
                    {
                        if(pdata[k].getPreviousyeartime() == agentsinfo[i][j].getgroversbankwithhlbseverityyearcount())
                        {
                            pcummulative5yearprofitdata <<pdata[k].getPreviousyearcummulative5yearprofit()<<"~"<<pdata[k].getPreviousyearstratergyname()<<"~"<<pdata[k].getPreviousyearstratergyparameter() <<"##";
                            pcummulativehlbzero5yearprofitdata <<pdata[k].getPreviousyearcummulative5yearhlbzeroprofit()<<"~"<<pdata[k].getPreviousyearstratergyname()<<"~"<<pdata[k].getPreviousyearstratergyparameter() <<"##";
                        }
                    }
                }
            }
            else
            {
                pcummulative5yearprofitdata << "0";
                pcummulativehlbzero5yearprofitdata << "0";
            }
             if(agents[i][j].getActionType() != 1){
                   
                    stringstream strategyNames;
                    stringstream strategyParams;
                    if (agents[i][j].behaviorPatterns.empty()) {
                        strategyNames << "NoAction";
                        strategyParams << "0;0;0;0;0;0;0";
                    }
                    else {
                        for (int k = 0; k < agents[i][j].behaviorPatterns.size(); k++) {
                            strategyNames << agents[i][j].behaviorPatterns[k]->getName();
                            strategyParams << agents[i][j].behaviorPatterns[k]->getParams();
                            if (k != agents[i][j].behaviorPatterns.size() - 1) {
                                strategyNames << "-";
                                strategyParams << "-";
                            }
                        }
                    }
                    outputFile << bioABM::getModelDay() << ","; 
                    outputFile << "g" << i << j << ",";
                    outputFile << agents[i][j].costs << ",";
                    outputFile << agents[i][j].returns << ",";
                    outputFile << (agents[i][j].returns - agents[i][j].costs) << ",";
                    outputFile << meanSeverity << ",";
                    outputFile << strategyNames.str() << ",";
                    outputFile << strategyParams.str() << ",";
                    outputFile << getDeadTrees(agents[i][j])<< ",";
                    outputFile << agentsinfo[i][j].getgroversbankprofit()<< ",";
                    if(agents[i][j].getActionType() == 1)
                        outputFile <<"Opt-Out"<< ","; 
                    else if(agents[i][j].getActionType() == 2)
                        outputFile <<"Repetition"<< ",";
                    else if(agents[i][j].getActionType() == 3)
                        outputFile <<"optimization"<< ",";
                    else if(agents[i][j].getActionType() == 4)
                        outputFile <<"Imitate"<< ",";
                    else
                        outputFile <<"Repetition"<< ",";

                    outputFile << agents[i][j].getSatisfaction() << ",";
                    outputFile << agents[i][j].getIncomeDissimilarity() << ",";
                    if(typeofnetwork == 1)
                        outputFile<<"Hub&Spoke-"<< grovernetwrokinfo[0][0]<< ",";
                    else if(typeofnetwork == 2)
                        outputFile<<"EveryBody"<< ",";
                    else if(typeofnetwork == 3)
                        outputFile<<"Neighbours"<< ",";
                    else
                        outputFile<<"EveryBody"<< ",";
                    outputFile << wl << ","; 
                    outputFile << wh << ",";
                    outputFile << pcummulative5yearprofitdata.str() << ",";
                    outputFile << pcummulativehlbzero5yearprofitdata.str() << ",";
                    outputFile << agentsinfo[0][0].getgroversbankinformation()<< ",";
                    outputFile << agentsinfo[0][1].getgroversbankinformation()<< ","; 
                    outputFile << agentsinfo[0][2].getgroversbankinformation()<< ","; 
                    outputFile << agentsinfo[1][0].getgroversbankinformation()<< ","; 
                    outputFile << agentsinfo[1][1].getgroversbankinformation()<< ","; 
                    outputFile << agentsinfo[1][2].getgroversbankinformation()<< ","; 
                    outputFile << agentsinfo[2][0].getgroversbankinformation()<< ","; 
                    outputFile << agentsinfo[2][1].getgroversbankinformation()<< ","; 
                    outputFile << agentsinfo[2][2].getgroversbankinformation()<< endl;  
                }
                else
                {
                        stringstream strategyNames;
                        stringstream strategyParams;
                        if (agents[i][j].behaviorPatterns.empty()) {
                            strategyNames << "NoAction";
                            strategyParams << "NA";
                        }
                        else {
                            for (int k = 0; k < agents[i][j].behaviorPatterns.size(); k++) {
                                strategyNames << agents[i][j].behaviorPatterns[k]->getName();
                                strategyParams << agents[i][j].behaviorPatterns[k]->getParams();
                                if (k != agents[i][j].behaviorPatterns.size() - 1) {
                                    strategyNames << "-";
                                    strategyParams << "-";
                                }
                            }
                        }
                        outputFile << bioABM::getModelDay() << ","; 
                        outputFile << "g" << i << j << ",";
                        outputFile << 0 << ",";
                        outputFile << 0 << ",";
                        outputFile << 0 << ",";
                        outputFile << 0 << ",";
                        outputFile << strategyNames.str() << ",";
                        outputFile << strategyParams.str() << ",";
                        outputFile << getDeadTrees(agents[i][j]) << ",";
                        outputFile << 0<< ",";
                        if(agents[i][j].getActionType() == 1)
                            outputFile <<"Opt-Out"<< ","; 
                        else if(agents[i][j].getActionType() == 2)
                            outputFile <<"Repetition"<< ",";
                        else if(agents[i][j].getActionType() == 3)
                            outputFile <<"optimization"<< ",";
                        else if(agents[i][j].getActionType() == 4)
                            outputFile <<"Imitate"<< ",";
                        else
                            outputFile <<"Repetition"<< ",";

                        outputFile << agents[i][j].getSatisfaction() << ",";
                        outputFile << agents[i][j].getIncomeDissimilarity() << ",";
                        if(typeofnetwork == 1)
                            outputFile<<"Hub&Spoke-"<< grovernetwrokinfo[0][0]<< ",";
                        else if(typeofnetwork == 2)
                            outputFile<<"EveryBody"<< ",";
                        else if(typeofnetwork == 3)
                            outputFile<<"Neighbours"<< ",";
                        else
                            outputFile<<"EveryBody"<< ",";
                        outputFile << wl << ","; 
                        outputFile << wh << ",";
                        outputFile << pcummulative5yearprofitdata.str() << ",";
                        outputFile << pcummulativehlbzero5yearprofitdata.str() << ","; 
                        outputFile << agentsinfo[0][0].getgroversbankinformation()<< ",";
                        outputFile << agentsinfo[0][1].getgroversbankinformation()<< ","; 
                        outputFile << agentsinfo[0][2].getgroversbankinformation()<< ","; 
                        outputFile << agentsinfo[1][0].getgroversbankinformation()<< ","; 
                        outputFile << agentsinfo[1][1].getgroversbankinformation()<< ","; 
                        outputFile << agentsinfo[1][2].getgroversbankinformation()<< ","; 
                        outputFile << agentsinfo[2][0].getgroversbankinformation()<< ","; 
                        outputFile << agentsinfo[2][1].getgroversbankinformation()<< ","; 
                        outputFile << agentsinfo[2][2].getgroversbankinformation()<< endl;  
                }
            }
    }
}





/*************************************************************
* Run Model
* Runs through the phases of the model until the simulation is
* complete
**************************************************************/
void runModel() {
    
    while (bioABM::getModelDay() <= bioABM::getModelDuration()) {
        // Stage 1: Psyllid Growth and Movement
        //cout<<"Before 1" << bioABM::getModelDay() <<endl;
        Phase1();
        
        //cout << "Period " << bioABM::getModelDay() << endl;
        if (!ParameterSet::biologicalRun) {
            if (bioABM::getModelDay() % ParameterSet::planningLength == 0) {
              //cout << "Planning period!\n";
            }

            //cout<<"Before 2" << bioABM::getModelDay() <<endl;
            // Stage 2: Execution of Planned Actions
            Phase2();
            
             
            // Stage 3: Behavior determination
            //Phase3();
            //cout<<"Before 4" << bioABM::getModelDay() <<endl;
            //Stage 4: Planning
            Phase4();
           

            //cout<<"Before 5" << bioABM::getModelDay() <<endl;
            //Stage 5: Accounting
            Phase5();
            
            if((bioABM::getModelDay()%365) == 0)
            {
                cout<<"Before 6" << bioABM::getModelDay() <<endl;
                //Stage 6: Grovers Prediciton of Surrouding
                Phase6();

                cout<<"Before 7" << bioABM::getModelDay() <<endl;
                //Stage 7: Work on the Optimization,Imitation
                Phase7();

                //cout<< bioABM::getModelDay()<<endl;
                
            }

            //cout<<"Before writeline" << bioABM::getModelDay() <<endl;
            writeCSVLine();

            //cout<<"Hello" << bioABM::getModelDay() <<endl;
            //std::cin.get();
        }
        
    }
  
    outputFile.close();
    
    bioABM::finishRun();
}

/******************************************************
* Parse Parameter File
* Parses a parameter file 
* ****************************************************/
void parseParameterFile(string fileName) {
    try {
        ifstream is(fileName);
        cereal::JSONInputArchive archive(is);
        archive(ParameterSet::planningLength, ParameterSet::freshYield, ParameterSet::juiceYield,
            ParameterSet::freshPrice, ParameterSet::juicePrice, ParameterSet::costs, ParameterSet::biologicalRun,
            ParameterSet::projectionLength, outputFilename, harvestDays, strategyParameters, strategyFlags, agencyFlags, experimentID,agentsbankinfo,
           memoryvalue,alphaplus,alphaminus,phiplus,phiminus,lambda,referenceincome,wl,wh,typeofnetwork);
    }
    catch (exception e) {
        cout << "ERROR WITH ECON JSON:" << e.what() << endl;
        exit(-1);
    }
}   

Commodity getCommodity() {
    //get harvests
    vector<string> sHarvests = split(harvestDays, ",");
    vector<int> harvests;
    for (string s: sHarvests) {
        harvests.push_back(stoi(s));
    }
    Commodity retVal(ParameterSet::freshYield, ParameterSet::juiceYield, ParameterSet::freshPrice, ParameterSet::juicePrice, ParameterSet::costs, harvests);
    return retVal;
}


/*****************************************************
* Main
*****************************************************/
int main(int argc, char ** argv) {
    string econConfigFile;
    string bioConfigFile;
    if (argc == 3) {
        econConfigFile = argv[1];
        bioConfigFile = argv[2];
        
    } 
    else {
        cout << "Using default filenames\n";
        econConfigFile = "configs/econConfig.json";
        bioConfigFile = "configs/bioConfig.json";
    }
    parseParameterFile(econConfigFile);
    bioABM::parseParameterFile(bioConfigFile);


    // Set up the behavior patterns, parameters arent used
    //IndividualAction individualSpray = IndividualAction(60, ParameterSet::sprayCost, ParameterSet::planningLength);
    
    //DEPRECATED: No pre-set behavior list for the time being

    /*
    GroupAction groupSpray = GroupAction(21, ParameterSet::sprayCost + 0, ParameterSet::planningLength,
                                        bioABM::getFallStart(), bioABM::getSpringStart(), bioABM::getSummerStart());
    NoAction noaction = NoAction();

    RogueTrees rougeing = RogueTrees(0, ParameterSet::rogueFrequency, ParameterSet::rogueRadius);

    //List of possible behaviors
    behaviorPatterns[0] = &noaction;
    behaviorPatterns[1] = &groupSpray;
    behaviorPatterns[2] = &rougeing;
    */
    bioABM::setExperimentID(experimentID);
    outputFile.open(outputFilename);
    outputFile
        << fixed << "t,id,costs,returns,profit,hlb_severity,strategy_names,strategy_params,RoguedTreecount,AnnualProfit,ActionType,Satisfaction,IncomeDisparity,TypeOfNetwork,Wl,Wh,P5yearcummulativeprofit,P5yearcummulativeHlbzeroprofit,G00Profit,G01Profit,G02Profit,G10Profit,G11Profit,G12Profit,G20Profit,G21Profit,G22Profit" << endl;

    prevdata.ReadPreviousData(pdata,0);
    cropvalue =getCommodity();
    InitialiseCHMA(cropvalue);


    //cout<<"Display the network connection" << endl;
    if(typeofnetwork == 1)
    {
        grovernetwrokinfo.resize(1, std::vector<string>(10));
        gnw.initializeBonds(grovernetwrokinfo,1,10);
        gnw.updateBondstypeone(grovernetwrokinfo, 1 ,10 );
        //gnw.displayBonds(grovernetwrokinfo,1,10);
    }
    else if(typeofnetwork == 2) 
    {
        grovernetwrokinfo.resize(9, std::vector<string>(10));
        gnw.initializeBonds(grovernetwrokinfo,9,10);
        gnw.updateBondstypetwo(grovernetwrokinfo,9,10);
        //gnw.displayBonds(grovernetwrokinfo,9,10);
    }
    else if(typeofnetwork == 3)
    {
        grovernetwrokinfo.resize(9, std::vector<string>(10));
        gnw.initializeBonds(grovernetwrokinfo,9,10);
        gnw.updateBondstypethree(grovernetwrokinfo,9,10);
        //gnw.displayBonds(grovernetwrokinfo,9,10);
    }

    runModel();
    

    

    return 0;
}
