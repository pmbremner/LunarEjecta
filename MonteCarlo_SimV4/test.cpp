#include <string>
#include <iostream>
#include <fstream>


using namespace std;




int main(int argc, char const *argv[])
{
	
	string filename = argv[1];

	cout << filename << endl;

	ofstream myfile(filename);

	myfile << 1 << endl;
	myfile.close(); 


	return 0;
}