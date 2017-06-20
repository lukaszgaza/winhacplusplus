#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <vector>
#include <stdlib.h>
#include <cmath>
#include "../Functions.h"


using namespace std;


void Tokenize(const string& str,
                      vector<string>& tokens,
                      const string& delimiters = " ");

int main(int argc,char **argv){
	if(argc<3) cout<<"Usage : Combine <output_path> [job 1 path] [job 2 path] ..."<<endl;

	cout<<"Combining jobs outputs"<<endl;
	int precision = 13;
	string path(argv[1]);

	map<string, vector<double> > xSections;
	map<string, vector<double> > xSectionsErrors;
	vector<double> numberOfEvents;
	vector<double> crudes;
	vector<double> crudesErrors;

	for(unsigned i = 2; i < argc ; ++i){
		string jobpath(argv[i]);
		fstream xsecfile((jobpath+"/XSections.dat").c_str(),fstream::in);
		fstream numberfile((jobpath+"/numberOfEvents.dat").c_str(),fstream::in);
		fstream crudefile((jobpath+"/XSectionCrude.dat").c_str(),fstream::in);

		cout<<" Reading folder: "<<jobpath<<endl;
		int eventsNo;
		numberfile>>eventsNo;
		numberfile.close();
		numberOfEvents.push_back(eventsNo);

		string tmpcrd;
		crudefile>>tmpcrd;
		vector<string> crdtokens;
		Tokenize(tmpcrd,crdtokens,":");
		if(crdtokens.size()==2){
			crudes.push_back(atof(crdtokens[0].c_str()));
			crudesErrors.push_back(atof(crdtokens[1].c_str()));
		}
		crudefile.close();

		while(!xsecfile.eof()){
			string tmp;
			xsecfile >> tmp;
			vector<string> line;
			Tokenize(tmp,line,":");
			if(line.size()!=3) continue;
			string name = line[0];
			double xsec = atof(line[1].c_str());
			double error = atof(line[2].c_str());
			cout<<"   "<<name<<" ";
			cout.precision(precision);
			cout<<xsec<<" ";
			cout.precision(precision);
			cout<<error<<endl;
			xSections[name].push_back(xsec);
			xSectionsErrors[name].push_back(error);
		}

		xsecfile.close();
	}


	double nTot = 0;
	map<string,double> sumW;
	map<string,double> sumW2;

	for(map<string, vector<double> >::iterator it = xSections.begin(); it!=xSections.end(); ++it){
		sumW[(*it).first] = 0;
		sumW2[(*it).first] = 0;
	}

	for(unsigned i = 0 ; i<numberOfEvents.size(); ++i){
		nTot += numberOfEvents[i];

		for(map<string, vector<double> >::iterator it = xSections.begin(); it!=xSections.end(); ++it){
			double w = numberOfEvents[i]*((*it).second[i]);
			double w2 = pow(w,2)*(pow(xSectionsErrors[(*it).first][i]/((*it).second[i]),2)+1.0/numberOfEvents[i] );
			sumW[(*it).first]+=w;
			sumW2[(*it).first]+=w2;
		}

	}

	string outpath(argv[1]);
	fstream totxsecfile((outpath+"/XSections.dat").c_str(),fstream::out);
	fstream totnumberfile((outpath+"/numberOfEvents.dat").c_str(),fstream::out);
	fstream totcrudefile((outpath+"/XSectionCrude.dat").c_str(),fstream::out);

	totnumberfile<<nTot;
	totnumberfile.close();

	for(map<string,double>::iterator it = sumW.begin(); it!=sumW.end(); ++it){
		totxsecfile.precision(precision);
		totxsecfile<<(*it).first<<":";
		totxsecfile.precision(precision);
		totxsecfile<<(*it).second/nTot<<":";
		totxsecfile.precision(precision);
		totxsecfile<<((*it).second/nTot)*sqrt(sumW2[(*it).first]/pow((*it).second,2)-1.0/nTot)<<endl;
	}

	totxsecfile.close();

	double totcrud=0;
	double totcruderr=0;

	for(unsigned i= 0 ; i<crudes.size(); ++i){
		totcrud += crudes[i];
		totcruderr += pow(crudesErrors[i],2);
	}
	totcrud = totcrud/static_cast<double>(crudes.size());
	totcruderr = sqrt(totcruderr)/static_cast<double>(crudes.size());

	totcrudefile<<totcrud<<":"<<totcruderr;
	totcrudefile.close();

	return 0;
}


void Tokenize(const string& str,
                      vector<string>& tokens,
                      const string& delimiters)
{
    // Skip delimiters at beginning.
    string::size_type lastPos = str.find_first_not_of(delimiters, 0);
    // Find first "non-delimiter".
    string::size_type pos     = str.find_first_of(delimiters, lastPos);

    while (string::npos != pos || string::npos != lastPos)
    {
        // Found a token, add it to the vector.
        tokens.push_back(str.substr(lastPos, pos - lastPos));
        // Skip delimiters.  Note the "not_of"
        lastPos = str.find_first_not_of(delimiters, pos);
        // Find next "non-delimiter"
        pos = str.find_first_of(delimiters, lastPos);
    }
}
