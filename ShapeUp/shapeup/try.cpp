#include <fstream>
#include <iostream>
using namespace std;
int main(){
    ifstream inFile;
    int x;
    inFile.open("../data/vertices_index.txt");
    for(int leg = 0; leg < 8; leg++){
        for(int i = 0; i < 4; i++){
            inFile >> x;
            cout << x << " ";
        }
        cout << endl;
    }
}