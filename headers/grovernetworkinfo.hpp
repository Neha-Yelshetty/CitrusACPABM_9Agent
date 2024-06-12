#ifndef GROVERNETWORKINFO_HPP
#define GROVERNETWORKINFO_HPP
#include <vector>
#include <sstream>
using namespace std;

class Grovernetworkinfo;

class Grovernetworkinfo {

    public:
        Grovernetworkinfo();
        void initializeBonds(std::vector<std::vector<string>>& bonds,int,int);
        void updateBondstypetwo(std::vector<std::vector<string>>& bonds,int ,int );
        void displayBonds(std::vector<std::vector<string>>& bonds,int ,int );
        void updateBondstypeone(std::vector<std::vector<string>>& bonds,int ,int );
        void updateBondstypethree(std::vector<std::vector<string>>& bonds,int ,int );
        string checkbondexists(std::vector<std::vector<string>>& bonds,int,int,int ,int,int,int);
        string checkbondexiststypetree(std::vector<std::vector<string>>& bonds,int,int,int ,int,int,int);

};

#endif