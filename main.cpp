#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>

#include "defs.h"
#include "vorticity.h"

using namespace std;

#define PARAM_COUNT 3

#define FL_RAW_BORDER 0

void print_eng_usage()
{
	cout << "This Software is for calculating vorticity between two points\n";
	cout << "USAGE:\tVorticity.exe <vp_out_file> <boundary_points_list> <results_file>\n\n";

	cout << "Example: ""Vorticity.exe out_2006-05-04_0730_n27799.m.pro_2006-05-04_1300_n70056.m.pro.txt boundaries.txt result.txt""\n\n"

		 << "\t<vp_out_file>\t\tVecPlotter output file\n"
		 	<< "\t    string format:\n"
				<< "\t\tgeo longitude start\n"
				<< "\t\tgeo latitude start\n"
				<< "\t\tgeo longitude end\n"
				<< "\t\tgeo latitude end\n"
				<< "\t\tpixel longitude start\n"
				<< "\t\tpixel latitude start\n"
				<< "\t\tpixel longitude end\n"
				<< "\t\tpixel latitude end\n"				
				<< "\t\tcorrelation\n"
				<< "\t\tvelocity of movement\n"
				<< "\t\ta priori error\n"

		 << "\t<boundary_points_list>\tList of cut definitions\n"
			<< "\t    string format:\n"
				<< "\t\tgeo longitude start\n"
				<< "\t\tgeo latitude start\n"
				<< "\t\tgeo longitude end\n"
				<< "\t\tgeo latitude end\n"
				<< "\t\tcut width (radius), [degrees]\tor -1 by default\n"
				<< "\t\tinterpolation interval (diameter), [degrees]\tor -1 by default\n"
				<< "\t\tweight coefficient\tor -1 by default\n"

		 << "\t<results_file>\t\tfile for result output\n"
			<< "\t    string format:\n"
				<< "\t\tgeo longitude start\n"
				<< "\t\tgeo latitude start\n"
				<< "\t\tgeo longitude end\n"
				<< "\t\tgeo latitude end\n"
				<< "\t\tcut width (radius), [degrees]\n"
				<< "\t\tinterpolation interval (diameter), [degrees]\n"
				<< "\t\tweight coefficient\n"
				<< "\t\tvector count\n"
				<< "\t\ta priori error\n"
				<< "\t\tinterpolation accuracy\n"
				<< "\t\tinterpolation mean-square deviation\n"
				<< "\t\tangle\n"
				<< "\t\tvorticity component -dV/dn [1/second]\n"
				<< "\t\talpha (regression line coefficient)\n"
				<< "\t\tleft boundary of alpha confidence interval\n"
				<< "\t\talpha right boundary\n"
				<< "\t\tcut length, [meters]\n"
				<< "\t\tinterpolation step count\n"
				<< "\t\tintegration step size [meters]\n"
		 << "\n";

}

char *movements_field_filename;
char *flow_line_filename;
char *flow_line_result_filename;
char *fl_piece_filename;
char *cut_filename;

char *result_filename;

bool file_exists(const char *fname)
{
	return ifstream(fname) != NULL;
}

bool analize_options(int argc, char** argv)
{
	if (argc != PARAM_COUNT + 1) return false;
	// flow_line_filename = argv[1];
	// flow_line_result_filename = argv[2];
	movements_field_filename = argv[1];
	cut_filename = argv[2];
	result_filename = argv[3];
	return /*file_exists(flow_line_filename) && */file_exists(movements_field_filename) && 
		   file_exists(cut_filename);
}

void read_flow_line(char *file_name, vector <flow_line_point> &fl)
{
	fstream ffl;
	ffl.open(flow_line_filename);
	int line, pixel, row;
	double longitude, latitude, geo;

	ofstream fgr;
	fgr.open("FL.vec");

	string s;
	getline(ffl, s);

	int max_row = 0;
	int count[15];
	memset(count, 0, sizeof(count));

	while(ffl >> line >> pixel >> longitude >> latitude >> row >> geo)
	{
		flow_line_point flp(line, pixel, longitude, latitude, row, geo);
		fl.push_back(flp);		

		max_row = max(max_row, row);
		++count[row];

		if (row > FL_RAW_BORDER) fgr << flp.toVecFormat() << endl;
	}

	ffl.close();
	fgr.close();
}

void write_nonzero_flow_line(char *file_name, vector <flow_line_point> &fl)
{
	ofstream fflr;
	fflr.open(file_name);

	cout << "fl size is " << fl.size() << endl;

	for (int i = 0; i < fl.size(); ++i)
		if (fl[i].row != 0)
			fflr << fl[i].toString() << endl;

	fflr.close();
}

void read_movements_field(char *filename, vector <movement> &m)
{
	fstream fmvn;
	fmvn.open(filename);

	double gsx, gsy, gex, gey, crl, vlc, err;
	int psx, psy, pex, pey;
	while(fmvn >> gsx >> gsy >> gex >> gey >> psx >> psy >> pex >> pey >> crl >> vlc >> err)
	{
		vec v(point(gsx, gsy), point(gex, gey));
		m.push_back(movement(v, crl, vlc, err));
	}

	fmvn.close();
}

void read_fl_pieces(char *filename, vector <flow_line_piece> &flp)
{
	fstream fflp;
	fflp.open(filename);

	double x1, y1, x2, y2;

	while (fflp >> x1 >> y1 >> x2 >> y2)
		flp.push_back(flow_line_piece(x1, y1, x2, y2));

	fflp.close();
}

void read_cuts(char *file_name, vector <scut> &cut)
{
	fstream fcut;
	fcut.open(file_name);
	double gsx, gsy, gex, gey, cut_width, itp_diameter, weight_coef;

	while(fcut >> gsx >> gsy >> gex >> gey >> cut_width >> itp_diameter >> weight_coef)
	{
		vec v(point(gsx, gsy), point(gex, gey));
		cut.push_back(scut(v, cut_width, itp_diameter, weight_coef));
	}
	fcut.close();
}

int main(int argc, char** argv)
{

	if (!analize_options(argc, argv))
	{
		cout << "Incorrect arguments!\n\n";
		print_eng_usage();
		return 0;
	}

	// vector <flow_line_point> flow_line;
	// read_flow_line(flow_line_filename, flow_line);

	// write_nonzero_flow_line(flow_line_result_filename, flow_line);

	// vector <flow_line_piece> fl_piece;
	// read_fl_pieces(fl_piece_filename, fl_piece);

	vector <scut> cut;
	read_cuts(cut_filename, cut);

	vector <movement> mvn;
	read_movements_field(movements_field_filename, mvn);

	Vorticity vorticity(mvn);

	ofstream fres;
	fres.open(result_filename);

	fres << vorticity_result::header() << endl;

	vector <vorticity_result> vort_res;

	for (int i = 0; i < cut.size(); ++i)
	{
		vort_res.clear();

		vorticity.set_cut(cut[i]);

		try
		{
			vorticity.calculate(vort_res);
		}
		catch(int a)
		{
			cout << "Caught exception number:  " << a << endl;
    		return -1;
		}
		
		for (int i = 0; i < vort_res.size(); ++i)
			fres << vort_res[i].toString() << endl;
	}

/*	Line l1(point(0, 0), point(1, 0));
	Line l2(point(0, 0), point(2, 1));
	cout << "angle " << l1.angle(l2.interval());*/

	// point p1(120.398, 77.1539), p2(129.55, 77.1804);
	// cout << vec(p1, p2).toString("for vec") << "\tlen is " << p1.distance_to(p2) << endl;

	fres.close();

	return 0;
}