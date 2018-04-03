#ifndef INTERPOLATION_H
#define INTERPOLATION_H

#include "defs.h"

#include <vector>
#include <algorithm>
#include <iterator>

#define WEIGHT_COEF 300.0

// #define DEFAULT_MS EMS_Degree

/*enum EMS//  measurement system
{
	EMS_Metter,
	EMS_Degree
};
*/

enum E_Output_Mode
{
	EOM_WITH_OUTPUT,
	EOM_WITHOUT_OUTPUT
};

class Interpolation
{
	int file_index;
	int calc_index;

	double R;
	double weight_coef;
	vec interval;
	vecs wv;
	vecs itp;

	E_Output_Mode output_mode;

	string itp_filename();
	string nc_filename();

	double get_norm_comp(int idx);

	double weight_func(double _r);
	double calc_weight_sum(vector <int> &idx, point pt, int omit_idx);
	double get_interpolation_result(vector <int> &idx, point pt, double sum);

public:
	Interpolation(vec itv, vecs &_wv, int _file_index, double _R, double wc);

	void set_file_index(int index);

	void calc_radius();
	void set_radius(double r);
	void set_weight_coef(double coef);

	double get_radius();
	double get_interval();
	double get_weight_coef();

	double take_for(point pt, int esc_ind);

	vecs take(int step_count, double *step);

	void calc_accuracy(double *itp_acc, double *mds, int *point_count, point start, point end);

};

string Interpolation::itp_filename()
{
	char str[20];
	sprintf(str, "%s%d_%d.%s", INTERPOLATION_FILENAME_PREFIX, file_index, calc_index, GLANCE_VECTOR_FILE_EXTENSION);
	return string(str);
}

string Interpolation::nc_filename()
{
	char str[20];
	sprintf(str, "%s%d_%d.%s", NORMAL_COMPONENT_FILENAME_PREFIX, file_index, calc_index, GLANCE_VECTOR_FILE_EXTENSION);
	return string(str);
}

double Interpolation::get_norm_comp(int idx)
{
	Line cut_line(interval);
	double ang = cut_line.angle(wv[idx]);
	return wv[idx].length()/* * - sign(ang)*/;
}

double Interpolation::weight_func(double r)
{
	return exp(- weight_coef * r * r) - exp(- weight_coef * R * R);
	// return 1 / (r * r);
}

double Interpolation::calc_weight_sum(vector <int> &idx, point pt, int omit_idx)
{
	double S = 0.0;
	for (int j = 0; j < wv.size(); ++j)
		if (omit_idx != j && pt.distance_to(wv[j].start) <= R /* 2*/)
		{
			idx.push_back(j);
			double r = pt.distance_to(wv[j].start);

			S += weight_func(r);
		}
	return S;
}

double Interpolation::get_interpolation_result(vector <int> &idx, point pt, double sum)
{
	if (sum == 0.0) return 0.0;
	double val = 0.0;
	for (int i = 0; i < idx.size(); ++i)
	{
		double r = pt.distance_to(wv[idx[i]].start);
		val += get_norm_comp(idx[i]) * weight_func(r) / sum;
	}
	return val;
}

Interpolation::Interpolation(vec itv, vecs &_wv, int _file_index, double _R, double wc) : 
	file_index(_file_index), calc_index(0), R(_R), weight_coef(wc), interval(itv), wv(_wv), output_mode(EOM_WITH_OUTPUT)
{
	if (_R == BY_DEFAULT) R = itv.length();
	else R /= 2;
	if (wc == BY_DEFAULT) weight_coef = WEIGHT_COEF_BY_DEFAULT;
}

void Interpolation::set_file_index(int index)
{
	file_index = index;
}

void Interpolation::calc_radius()
{
	double _min, _max, res;
	vector <double> dist;
	_min = _max = interval.start.distance_to(wv[0].start);
	for (int i = 0 ; i < wv.size(); ++i)
	{
		double d = interval.start.distance_to(wv[i].start);
		dist.push_back(d);
		// _min = min(_min, d);
		// _max = max(_max, d);
	}
	sort(dist.begin(), dist.end());
	// res = max(dist.front(), dist_pp(st.start, st.end) - dist.back());
	res = 0.0;
	for (int i = 1; i < dist.size(); ++i)
	{
		res = max(res, dist[i] - dist[i - 1]);
	}
	R = res * 1.5;
}

void Interpolation::set_radius(double _r)
{
	R = _r;
}

void Interpolation::set_weight_coef(double coef)
{
	weight_coef = coef;
}


double Interpolation::get_radius()
{
	return R;
}

double Interpolation::get_interval()
{
	return R * 2;
}

double Interpolation::get_weight_coef()
{
	return weight_coef;
}

double Interpolation::take_for(point pt, int esc_ind = -1)
{

	vector <int> act_p_ind;

	double S = calc_weight_sum(act_p_ind, pt, esc_ind);

	double val = get_interpolation_result(act_p_ind, pt, S);

	return val;
}

vecs Interpolation::take(int step_count, double *step)
{
	if (step_count == 0) return itp;

	int n = step_count;
	double len = interval.length();
	double h = len / n;
	double dx = (interval.end.x - interval.start.x) / n;
	double dy = (interval.end.y - interval.start.y) / n;

	*step = interval.length_m() / n;

	Line cut_line(interval);

	for (int i = 0; i < n; ++i)
	{
		point gr1(interval.start.x + i * dx, interval.start.y + i * dy);
		point gr2(interval.start.x + (i + 1) * dx, interval.start.y + (i + 1) * dy);
		point gr_avr = vec(gr1, gr2).middle();

		double val = take_for(gr_avr);
		vec nv = cut_line.perpendicular(gr_avr, LS2deg(val));

		itp.push_back(nv);

	}

	return itp;
}


void Interpolation::calc_accuracy(double *itp_acc, double *msd, int *point_count, point start, point end)
{
	++calc_index;

	vector <double> vl;
	vector <double> nc;

	vector <double> err;

	Line itv(start, end);

	double acr_sum = 0.0;
	for (int i = 0; i < wv.size(); ++i)
	{
		if (itv.contains(wv[i].start))
		{
			double val = take_for(wv[i].start, i);		
			err.push_back(val - get_norm_comp(i));
			acr_sum += err.back();

			vl.push_back(val);
			nc.push_back(get_norm_comp(i));	
		}
		
	}
	*point_count = err.size();
	if (*point_count == 0) return;

	*itp_acc = acr_sum / *point_count;

	double msd_sum = 0.0;
	for (int i = 0; i < err.size(); ++i)
		msd_sum += (*itp_acc - err[i]) * (*itp_acc - err[i]);
	*msd = sqrt(msd_sum / *point_count);
/*
	char str[15];
	sprintf(str, "err%d.txt", err.size());
	ofstream f;
	f.open(str);
	copy(err.begin(), err.end(), ostream_iterator<double>(f, " "));
	f.close();
*/
}

#endif // INTERPOLATION_H