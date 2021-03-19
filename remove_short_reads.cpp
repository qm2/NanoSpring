#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>

int main(int argc, char **argv)
{
	std::string infile = std::string(argv[1]);
    int min_readlen = atoi(argv[2]);
	std::string s[4];
	std::string outfile = infile + ".filtered";	
	std::ifstream f_in(infile);
	std::ofstream f_out(outfile);
	while(std::getline(f_in,s[0]))
	{
		std::getline(f_in,s[1]);
		std::getline(f_in,s[2]);
		std::getline(f_in,s[3]);
		if(s[1].length() >= min_readlen)
			f_out << s[0] <<"\n" <<  s[1] <<"\n" <<  s[2] <<"\n" <<  s[3] <<"\n"; 
	}
	f_in.close();
	f_out.close();
	return 0;	
}
