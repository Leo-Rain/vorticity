#ifndef DEFS_H
#define DEFS_H

#include "geometry.h"

#define BREADTH_SEARCH_RADIUS 0.1 // degrees

#define NORMAL_COMPONENT_FILENAME_PREFIX "NC"
#define PROJECTION_FILENAME_PREFIX 		 "PNC"
#define INTERPOLATION_FILENAME_PREFIX 	 "ITP"
#define APPROX_LINE_PREFIX 				 "AL"
#define GLANCE_VECTOR_FILE_EXTENSION 	 "vec"

#define BY_DEFAULT -1

#define WEIGHT_COEF_BY_DEFAULT 200.0

#define BLACK	0
#define VINOUS  2
#define VIOLET	5
#define RED		10
#define YELLOW	11
#define GREEN	12
#define BLUE	14

#define			   CUT_COLOR VINOUS
#define		 NORM_COMP_COLOR VIOLET
#define NORM_COMP_PROJ_COLOR VINOUS
#define    APPROX_LINE_COLOR GREEN
#define		 FLOW_LINE_COLOR VIOLET
#define  INTERPOLATION_COLOR YELLOW

#define NO_ANCHOR "NO_ANCHOR"

typedef vector <vec> vecs;
typedef vector <point> points;

ofstream flog;

struct movement
{
	vec v;
	double correlation;
	double velocity;
	double error;		// априорная ошибка

	movement(vec _v, double crl, double vl, double err) : v(_v), correlation(crl), velocity(vl), error(err) {}

	string toString()
	{
		char str[70];
		sprintf(str, "%2f %2f %2f %2f %2f %2f %2f", v.start.x, v.start.y, v.end.x, v.end.y, correlation, velocity, error);
		return string(str);
	}

	static string header() 
	{ 
		char str[150];
		sprintf(str, "longitude1 latitude1 longitude2 latitude2 correlat velocity a_apriori_error");
		return string(str); 
	}
};

struct flow_line_point
{
	int line;
	int pixel;
	double longitude;
	double latitude;
	int row;
	double geo;

	flow_line_point(int ln, int pxl, double lgt, double ltt, int _row, double _geo) :
		line(ln), pixel(pxl), longitude(lgt), latitude(ltt), row(_row), geo(_geo) {}

	string toString()
	{
		char str[50];
		sprintf(str, "%d %d %2f %2f %d %2f", line, pixel, longitude, latitude, row, geo);
		return string(str);
	}

	string toVecFormat()
	{
		char str[80];
		sprintf(str, "TYPE = POINT\tCOLOR = 12\tNATION = RU\tLEGEND=\"\"\tGEO=(%2f,%2f)\tWIDTH = 1", longitude, latitude);
		return string(str);
	}
};

struct flow_line_piece
{
	vec v;

	flow_line_piece(double x1, double y1, double x2, double y2) : 
		v(point(x1, y1), point(x2, y2)) {}

};

struct scut // Разрез
{
	point start;
	point end;
	double width;
	double itp_interval;
	double weight_coef;

	// cut(vec v, double w) : 
		// start(v.start), end(v.end), width(w) {}

	scut() : width(0.0), itp_interval(0.0), weight_coef(0.0) {};
	scut(vec v, double w, double ii = 0.0, double wc = 0.0) : 
		start(v.start), end(v.end), width(w),
		itp_interval(ii), weight_coef(wc) {}

	vec v()	{ return vec(start, end); }

	void set_v(vec _v)
	{
		start = _v.start;
		end = _v.end;
	}

	string toString(string name = string(""))
	{
		char str[30];
		sprintf(str, "%s(%2f, %2f)", name.c_str(), width, weight_coef);
		return string(str);
	}
};


//--------------------------------------BORDER---------------------------------------------------

class Border
{
	double to_start;
	double to_end;
	point start;
	point end;
	vec cut;

public:
	Border(vec v);
	void refresh(point p);
	vec get_border();
};

Border::Border(vec v) : cut(v)
{
	to_start = to_end = v.length();
	start = v.end;
	end = v.start;
}

void Border::refresh(point p)
{
	if (p.distance_to(cut.start) < to_start)
	{
		to_start = p.distance_to(cut.start);
		start = p;
	}
	if (p.distance_to(cut.end) < to_end)
	{
		to_end = p.distance_to(cut.end);
		end = p;
	}
}

vec Border::get_border()
{
	return vec(start, end);
}

#endif // DEFS_H