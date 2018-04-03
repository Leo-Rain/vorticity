#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <fstream>
#include <string>

#define EPS 0.0000001
#define GR 4e7 / 360 	// метров в градусе

#define EARTH_RADIUS 6372795 // in meters

#define sign_(a) (a == 0.0) ? 0.0 : (a < 0.0 ? -1.0 : 1.0)
#define deg2rad(a) a * M_PI / 180
#define rad2deg(a) a * 180 / M_PI
#define loc_sys(a) a 
#define deg2LS(a) a /** GR*/
#define LS2deg(a) a /*/ GR*/

using namespace std; //	delete (escape it)

/*double sign(double a)
{
	if (a == 0.0) return 0.0;
	if (a >= 0.0) return 1.0;
	if (a <= 0.0) return -1.0;
}*/

struct point
{
	double x; // longitude
	double y; // latitude

	point() : x(0), y(0) {};
	point(double a, double b) : x(a), y(b) {};
	// point(const point &p) : x(p.x), y(p.y) {};

	bool equal(const point &p)
	{
		return (x == p.x && y == p.y);
	}

	double distance_to(point p)
	{
		return sqrt((x - p.x) * (x - p.x) + (y - p.y) * (y - p.y));
		double f1 = deg2rad(y), l1 = deg2rad(x), f2 = deg2rad(p.y), l2 = deg2rad(p.x);
		double df2 = (f2 - f1) / 2, dl2 = (l2 - l1) / 2;
		double ds = 2 * asin(sqrt(sin(df2) * sin(df2) + cos(f1) * cos(f2) * sin(dl2) * sin(dl2)));
		// cout << ds << "\t" << point(x, y).toString("p1") << "\t" << p.toString("p2") << "\t";
		// cout << ds * EARTH_RADIUS << endl;
		return loc_sys(ds * EARTH_RADIUS);
	}

	double distance_to_m(point p)
	{
		double f1 = deg2rad(y), l1 = deg2rad(x), f2 = deg2rad(p.y), l2 = deg2rad(p.x);
		double df2 = (f2 - f1) / 2, dl2 = (l2 - l1) / 2;
		double ds = 2 * asin(sqrt(sin(df2) * sin(df2) + cos(f1) * cos(f2) * sin(dl2) * sin(dl2)));
		// cout << ds << "\t" << point(x, y).toString("p1") << "\t" << p.toString("p2") << "\t";
		// cout << ds * EARTH_RADIUS << endl;
		return loc_sys(ds * EARTH_RADIUS);
	}

	string toString(string name = string(""))
	{
		char str[30];
		sprintf(str, "%s(%2f, %2f)", name.c_str(), x, y);
		return string(str);
	}

	string toStr()
	{
		char str[20];
		sprintf(str, "%2f %2f", x, y);
		return string(str);
	}

	static string header()
	{
		return string("longitude  latitude");
	}

};

struct vec
{
	point start;
	point end;

	vec() : start(), end() {};
	vec(point a, point b) : start(a), end(b) {};

	double dx()
	{
		return end.x - start.x;
	}

	double dy()
	{
		return end.y - start.y;
	}

	bool contains(point p)
	{
		return sign_(p.x - start.x) * sign_(p.x - end.x) < 0;
	}

	double length()
	{
		return start.distance_to(end);
	}

	double length_m()
	{
		return start.distance_to_m(end);
	}

/*	double sign_vv(vec ln, vec v)
	{
		double l1 = dist_vp(ln, v.start), l2 = dist_vp(ln, v.end);
		// if ((p2v(ln, v.start) < 0 && p2v(ln, v.end) > 0) || sign(p2v(ln, v.start)) * sign(l2 - l1) >= 0) return 1.0;
		if ((p2v(ln, v.start) < 0 && (p2v(ln, v.end) > 0 || (p2v(ln, v.end) < 0 && l2 - l1 > 0)))) return 1.0;
		if (p2v(ln, v.start) > 0 && l2 - l1 > 0) return 1.0;
		return -1.0;
	}*/

	void shorten(double len/*, double dl = 0.0*/)
	{
		double m = len, 
		/*if (dl == 0.0)*/ n = LS2deg(start.distance_to(end)) - len;
		/*else n = dl;*/
		double x1 = start.x, y1 = start.y, x2 = end.x, y2 = end.y;
		double l = m / n;
		double x = (x1 + l * x2) / (1 + l),
			   y = (y1 + l * y2) / (1 + l);
		end = point(x, y);
		return;
	}

	point middle()
	{
		return point((start.x + end.x) / 2, (start.y + end.y) / 2);
	}

	static string header()
	{
		return string("longitude1 latitude1 longitude2 latitude2");
	}

	string toString(string name = string(""))
	{
		char str[100];
		sprintf(str, "%s[%s -> %s]", name.c_str(), start.toString().c_str(), end.toString().c_str());
		return string(str);
	}

	string toStr()
	{
		char str[50];
		sprintf(str, "%s %s", start.toStr().c_str(), end.toStr().c_str());
		return string(str);
	}

	string toGlanceFormat(int color = 10, int width = 1, string end_cap = string("ARROW_ANCHOR"))
	{
		char str[120];
		sprintf(str, "TYPE = VECTOR\tCOLOR = %d\tWIDTH = %d\tSCALE = 1.00\tEND_CAP = %s\tGEO = (%2f,%2f %2f,%2f)\n",
					 color, width, end_cap.c_str(), start.x, start.y, end.x, end.y);
		return string(str);
	}
};

class Line // line segment
{
	double a, b, c;
	point p1, p2;

	double f(point p);
public:
	Line();
	Line(vec v);
	Line(point _p1, point _p2);
	Line(double _a, double _b, double _c); // p1 and p2 are not initialized

	double get_a();
	double get_b();
	double get_c();

	point start();
	point end();

	string toString();

	double x2y(double _x);

	vec interval();
	void set_interval(point _p1, point _p2);

	void set_bounds(Line &l);	// for Line(_a, _b, _c)

	bool contains(point p);
	bool contains_projection(point p);

	double distance_to(point p);
	double distance_to_m(point p);
	double angle(vec v);

	point intersection(Line l);
	point projection(point p);

	vec perpendicular(point p, double r = 0.0);
	vec normal_component(vec v);
	vec normal_component_projection(vec v);
	Line parallel(point p);
};

Line::Line() : a(0), b(0), c(0), p1(point(0, 0)), p2(point(0, 0)) {}

Line::Line(vec v) : p1(v.start), p2(v.end)
{
	a = p2.y - p1.y;
	b = p1.x - p2.x;
	c = - p1.y * b - p1.x * a;
}

Line::Line(point _p1, point _p2) : p1(_p1), p2(_p2)
{
	a = p2.y - p1.y;
	b = p1.x - p2.x;
	c = - p1.y * b - p1.x * a;
}

Line::Line(double _a, double _b, double _c) : a(_a), b(_b), c(_c), p1(point(0, 0)), p2(point(0, 0))
{
	// p1 and p2 are NOT INITIALIZED
}

double Line::f(point p)
{
	return a * p.x + b * p.y + c;
}

double Line::get_a()
{
	return a;
}

double Line::get_b()
{
	return b;
}

double Line::get_c()
{
	return c;
}

point Line::start()
{
	return p1;
}

point Line::end()
{
	return p2;
}

string Line::toString()
{
	char str[100];
	sprintf(str, "a = %2f, b = %2f, c = %2f", a, b, c);
	return string(str);
}

double Line::x2y(double _x)
{
	if (b == 0.0) return 0.0;
	return -(a * _x + c) / b;
}

vec Line::interval()
{
	return vec(p1, p2);
}

void Line::set_interval(point _p1, point _p2)
{
	// Line(_p1, _p2);
	p1 = _p1; p2 = _p2;
	a = p2.y - p1.y;
	b = p1.x - p2.x;
	c = - p1.y * b - p1.x * a;
}

void Line::set_bounds(Line &l)	// for Line(_a, _b, _c)
{
	Line st_prnd(l.perpendicular(l.start()));
	p1 = intersection(st_prnd);
	Line end_prnd(l.perpendicular(l.end()));
	p2 = intersection(end_prnd);
}

bool is_t(double numerator, double denominator)
{
	return (denominator == 0 && numerator == 0) || 
		   (denominator != 0 && numerator / denominator >= 0 && numerator / denominator <= 1);
}

bool Line::contains(point p)
{
	return fabs(f(p)) < EPS && is_t(p.x - p1.x, p2.x - p1.x) && is_t(p.y - p1.y, p2.y - p1.y);
}

bool Line::contains_projection(point p)
{
	point pt = projection(p);
	return contains(pt);
}

double Line::distance_to(point p)
{
	return p.distance_to(projection(p));
	// return fabs(f(p)) / sqrt(a * a + b * b);
}

double Line::distance_to_m(point p)
{
	return p.distance_to_m(projection(p));
	// return fabs(f(p)) / sqrt(a * a + b * b);
}

double Line::angle(vec v)
{
	Line l(v);
	return rad2deg(atan2(a * l.get_b() - l.get_a() * b, a * l.get_a() + b * l.get_b()));
}

point Line::intersection(Line l)
{
	double y = (a * l.get_c() - l.get_a() * c) / (l.get_a() * b - a * l.get_b());
	double x = (b * l.get_c() - l.get_b() * c) / (a * l.get_b() - l.get_a() * b);
	return point(x, y);
}

point Line::projection(point p)
{
	Line l = perpendicular(p);
	return intersection(l);
}

vec Line::perpendicular(point p, double r)
{
	vec pv1(p, point(p.x + a, p.y + b));
	// vec pv2(p, point(p.x - a, p.y - b));
	if (r > 0)
	{
		pv1.shorten(r);
		// cout << ">";
		// pv2.shorten(r);
	}
	return vec(p, pv1.end);
}

vec Line::normal_component(vec v)
{
	point prj = projection(v.start);
	Line prl = parallel(v.end);
	point norm = prl.projection(v.start);
	return vec(v.start, norm);
}

vec Line::normal_component_projection(vec v)
{
	vec nc = normal_component(v);
	point st_proj = projection(v.start);
	return vec(st_proj, point(st_proj.x + nc.dx(), st_proj.y + nc.dy()));
}

Line Line::parallel(point p)
{
	return Line(get_a(), get_b(), -(get_a() * p.x + get_b() * p.y));
}

#endif // GEOMETRY_H