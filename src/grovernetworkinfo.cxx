#include "../headers/grovernetworkinfo.hpp"
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include <iomanip>
#include <math.h>
#include <algorithm>
#include <random>
#include <boost/random.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>
#include <boost/assign.hpp>

Grovernetworkinfo::Grovernetworkinfo() {
   
}

void Grovernetworkinfo:: initializeBonds(std::vector<std::vector<string>>& bonds,int numrows,int columns)
{
        for (int i = 0; i < numrows; ++i) {
                for (int j = 0; j < columns; ++j) {
                // Initialize all bonds to false (no bond)
                        bonds[i][j] = "";
                }
         }
}
int intRand(const int& min, const int& max) {
    static thread_local std::mt19937 generator(std::random_device{}());
    std::uniform_int_distribution<int> distribution(min, max);
    return distribution(generator);
}

void Grovernetworkinfo:: updateBondstypetwo(std::vector<std::vector<string>>& bonds,int numrows,int columns) {

    std::vector<std::vector<string>> tbonds(numrows, std::vector<string>(columns-1));
    int value1;
    int value2;
    for(int j=0 ;j <numrows;++j){
        // Calculate row and column indices from element index
        int row = j / 3;
        int col = j % 3;
        //bonds[j][0] = std::to_string(j/3) + std::to_string(j%3);
        // Loop through all elements to update bonds
        for (int i = 0; i < columns-1;++i) {
            // Calculate row and column indices for the connected element
            int connectedRow = i / 3;
            int connectedCol = i % 3;

            // Check if the element is not the same as the current one
            if (row != connectedRow || col != connectedCol) {
                // Update bond connections for both directions
                tbonds[j][i] = "true";
                tbonds[i][j] = "true";
            }
        }
    }
    for (int i = 0; i < numrows; ++i) {
        //bonds[i][0] = std::to_string(i/3) + std::to_string(i%3);
        for (int j = 0; j < columns-1;++j) {
            if (tbonds[i][j] == "true") {
                value1 = j / 3;
                value2 = j % 3;
                tbonds[i][j] = std::to_string(value1) + std::to_string(value2);
            }
        }
    }

    for (int i = 0; i < numrows; ++i) {
        value1 = i / 3;
        value2 = i % 3;
        bonds[i][0] = std::to_string(value1) + std::to_string(value2);
        for (int j = 0; j < columns-1; ++j) {
             
             bonds[i][j+1] = tbonds[i][j];
                //cout << tbonds[i][j] << " ";
             
        }
        //cout<<endl;
    }
}


void Grovernetworkinfo:: updateBondstypeone(std::vector<std::vector<string>>& bonds,int numrows,int columns) {
    int randomNumber = intRand(0,8);
    int value1 = randomNumber / 3;
    int value2 = randomNumber % 3;
    bonds[0][0] = std::to_string(value1) + std::to_string(value2);
    for (int j = 0; j < columns-1; ++j) {
        value1 = j / 3;
        value2 = j % 3;
        if (bonds[0][0] != (std::to_string(value1) + std::to_string(value2)))
            bonds[0][j+1] = std::to_string(value1) + std::to_string(value2);
    }
}

void Grovernetworkinfo:: updateBondstypethree(std::vector<std::vector<string>>& bonds,int numrows,int columns) {
    
    /*bonds[0][0] = {"01","10","11"};
    bonds[0][1] = {"00","02","10","11","12"};
    bonds[0][2] = {"01","11","12"};
    bonds[1][0] = {"00","01","11","21","20"};
    bonds[1][1] = {"00","01","02","10","12","20","21","22"};
    bonds[1][2] = {"02","01","11","21","22"};
    bonds[2][0] = {"10","11","21"};
    bonds[2][1] = {"20","10","11","12","22"};
    bonds[2][2] = {"12","11","21"};*/

    bonds.resize(numrows);
    bonds[0] = {"00","01", "10", "11"};
    bonds[1] = {"01","00", "02", "10", "11", "12"};
    bonds[2] = {"02","01", "11", "12"};
    bonds[3] = {"10","00","01","11","21","20"};
    bonds[4] = {"11","00","01","02","10","12","20","21","22"};
    bonds[5] = {"12","02","01","11","21","22"};
    bonds[6] = {"20","10","11","21"};
    bonds[7] = {"21","20","10","11","12","22"};
    bonds[8] = {"22","12","11","21"};

}



void Grovernetworkinfo:: displayBonds(std::vector<std::vector<string>>& bonds,int numrows,int columns) {

    for (int i = 0; i < numrows; ++i) {
        for (int j = 0; j < columns; ++j) {
             
                cout << bonds[i][j] << " ";
             
        }
        cout<<endl;
    }
}

string Grovernetworkinfo:: checkbondexiststypetree(std::vector<std::vector<string>>& bonds,int numrows,int columns,int ownerrow,int ownercol,int friendrow,int friendcol) {

    for (int i = 0; i < numrows; ++i) {
        if(static_cast<int>(bonds[i][0][0]-'0') == ownerrow && static_cast<int>(bonds[i][0][1]-'0') == ownercol)
        {
            for (int j = 0; j < columns; ++j) {
                if(static_cast<int>(bonds[i][j][0]-'0') == friendrow && static_cast<int>(bonds[i][j][1]-'0') == friendcol)
                {
                    return "Yes";
                }
            }
        }
    }
    return "No";
}


string Grovernetworkinfo:: checkbondexists(std::vector<std::vector<string>>& bonds,int numrows,int columns,int ownerrow,int ownercol,int friendrow,int friendcol) {

    for (int i = 0; i < numrows; ++i) {
        //cout<<bonds[i][0][0] << "~~" << bonds[i][0][1] << "~~" << ownerrow <<"~~"<< ownercol << "~~"<< friendrow << "~~"<< friendcol<< "///";
        if(static_cast<int>(bonds[i][0][0]-'0') == friendrow && static_cast<int>(bonds[i][0][1]-'0') == friendcol)
        {
            return "Yes";
        }
        if(static_cast<int>(bonds[i][0][0]-'0') == ownerrow && static_cast<int>(bonds[i][0][1]-'0') == ownercol)
        {
            for (int j = 0; j < columns; ++j) {
                //cout << i << "-"<<j <<"-"<<bonds[i][j]<<"-"<<static_cast<int>(bonds[i][j][0]-'0')<<"-"<<static_cast<int>(bonds[i][j][1]-'0')<<":";
                if(static_cast<int>(bonds[i][j][0]-'0') == friendrow && static_cast<int>(bonds[i][j][1]-'0') == friendcol)
                {
                    return "Yes";
                }
                
            }
        }
    }
    return "No";
}

