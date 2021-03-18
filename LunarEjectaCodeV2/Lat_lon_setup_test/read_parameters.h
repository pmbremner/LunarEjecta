#ifndef READ_PARAMETERS_H
#define READ_PARAMETERS_H

#include <fstream>
#include <iostream>
#include <sstream>

using namespace std;

// https://www.fluentcpp.com/2019/07/23/how-to-define-a-global-constant-in-cpp/
inline bool initError = 0;

// https://www.cplusplus.com/doc/oldtutorial/templates/
template <class paramType>
void getParam (string param_fn, string paramLabel, paramType& param, bool defaultExists) {
	ifstream input_file;

	input_file.open("./param_files/" + param_fn);
	char C_Line[64];
	paramType tparam;
	bool isError = 1;

	int lineNumber = 1;

	while (input_file.getline(C_Line, 64, '#'))
	{
		string paramLabel_file(C_Line);
		
		if(paramLabel_file == paramLabel)
		{
			if(input_file.getline(C_Line, 64, '#'));
				isError = 0;
			string tstr(C_Line);
			stringstream SS_param(tstr);
			SS_param >> param;
			//input_file.ignore(256, '\n'); // temp
			if(defaultExists)
				cout << "  Overriding Default value...\n";
			cout << "Line " << lineNumber << " | " << paramLabel << " = " << param << endl;
			//return;
		}
		else
		{
			input_file.ignore(256, '\n');
		}

		lineNumber++;
	}
	// if there's an error, set the initError flag
	if(!defaultExists && isError){
		initError = 1;
		cout << "ERROR: Parameter '" << paramLabel << "' not found...\n";
	}

	input_file.close();
}


#endif 