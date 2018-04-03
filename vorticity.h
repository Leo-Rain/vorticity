#ifndef VORTICITY_H
#define VORTICITY_H

#include <algorithm>

#include "interpolation.h"
#include "lsm.h"
#include "defs.h"
#include "geometry.h"

struct line_index
{
	int begin;
	int end;

	line_index() : begin(0), end(0) {}

	line_index(int a, int b) : begin(a), end(b) {}
};

struct vorticity_result
{
	scut cut;
	int vec_count;
	double a_priori_error;
	double itp_accuracy;
	double msd;
	double angle;
	double _dV_dn;
	double approx_line_alpha;
	double approx_line_alpha_left;
	double approx_line_alpha_right;
	double approx_line_accuracy;
	double approx_line_average;
	double approx_line_msd;
	int itp_step_count;
	double itp_step_size;

	vorticity_result() : vec_count(0), a_priori_error(1.0), itp_accuracy(0.0), msd(0.0), angle(0.0), _dV_dn(0.0), approx_line_alpha(0.0), approx_line_alpha_left(0.0), 
						 approx_line_alpha_right(0.0), approx_line_accuracy(0.0), approx_line_average(0.0), approx_line_msd(0.0),
						 itp_step_count(0), itp_step_size(0.0) {}

	vorticity_result(vec v, double width) : cut(v, width, 0.0, 0.0), vec_count(0), a_priori_error(1.0), itp_accuracy(0.0), msd(0.0), angle(0.0), _dV_dn(0.0), 
						 approx_line_alpha(0.0), approx_line_alpha_left(0.0), 
						 approx_line_alpha_right(0.0), approx_line_accuracy(0.0), approx_line_average(0.0), approx_line_msd(0.0),
						 itp_step_count(0), itp_step_size(0.0) {}

	vorticity_result(vec v, double width, double itp_itv, double weight_coef, int _vec_count) : cut(v, width, itp_itv, weight_coef), vec_count(_vec_count), a_priori_error(1.0), 
						 itp_accuracy(0.0), msd(0.0), angle(0.0), _dV_dn(0.0), 
						 approx_line_alpha(0.0), approx_line_alpha_left(0.0), 
						 approx_line_alpha_right(0.0), approx_line_accuracy(0.0), approx_line_average(0.0), approx_line_msd(0.0),
						 itp_step_count(0), itp_step_size(0.0) {}

	vorticity_result(vec v, double width, double itp_itv, double weight_coef, int _vec_count, double apr_err, double itp_acc, double _msd, double _angle, double dVn,
							double apl_alpha, double apl_alphal, double apl_alphar, double apl_acc, double apl_avr, double apl_msd, int isc, double iss) : 
		cut(v, width, itp_itv, weight_coef), vec_count(_vec_count), a_priori_error(apr_err), itp_accuracy(itp_acc), msd(_msd), angle(_angle), _dV_dn(dVn),
		approx_line_alpha(apl_alpha), approx_line_alpha_left(apl_alphal), approx_line_alpha_right(apl_alphar), 
		approx_line_accuracy(apl_acc), approx_line_average(apl_avr), approx_line_msd(apl_msd), itp_step_count(isc), itp_step_size(iss) {}

	static string header()
	{
		char str[200];
		sprintf(str, "longitude1 latitude1 longitude2 latitude2 cut_width itp_itv weight_K vecs apr_err itp_acc  itp_MSD   angle    -dV/dn    alpha     left     right\t cut_len   steps step_size");
		return string(str);
	}

	string toString()
	{
		char str[220];
		memset(str, 0, sizeof(str));
		sprintf(str, "%2f %2f %2f %2f %2f %2f %2f %d %2f %2f %2f %2f %2f %2f %2f %2f\t%2f %d %2f", cut.start.x, cut.start.y, cut.end.x, cut.end.y, 
			cut.width, cut.itp_interval, cut.weight_coef, vec_count, a_priori_error, itp_accuracy, msd, angle, _dV_dn, approx_line_alpha, approx_line_alpha_left,
			approx_line_alpha_right, /*approx_line_accuracy, approx_line_average, approx_line_msd,*/ cut.v().length_m(), itp_step_count, itp_step_size);
		return string(str);
	}
};

class Vorticity
{
	scut cut;

	int calc_index;
	int sub_index;
	vector <movement> mvn;

	string filename(string prefix);
	string subfilename(string prefix);

	void find_extremum(vecs &vs, vector <line_index> &idx, vecs &crd);

	void output_graphic_info(vector <int> &idx, vecs &pnc, vecs &itp, point begin, point end, int start_idx, int end_idx, double *apr_err);

public:
	Vorticity(vector <movement> &m);
	~Vorticity();

	void set_cut(scut _cut);
	void set_fl_segment(vec v);
	// void calc_prnd();

	string get_result_header();

	void calculate(vector <vorticity_result> &result);
};

string Vorticity::filename(string prefix)
{
	char str[30];
	sprintf(str, "%s%d.%s", prefix.c_str(), calc_index, GLANCE_VECTOR_FILE_EXTENSION);
	return string(str);
}

string Vorticity::subfilename(string prefix)
{
	char str[30];
	sprintf(str, "%s%d_%d.%s", prefix.c_str(), calc_index, sub_index, GLANCE_VECTOR_FILE_EXTENSION);
	return string(str);
}

void Vorticity::find_extremum(vecs &vs, vector <line_index> &idx, vecs &crd)
{
	if (vs.size() == 0) return;

	line_index li;
	li.begin = 0;

	if (vs.size() < 2) return;
	vec vextreme;
	int ext_count = 0;
	double prev_level = vs[0].length();
	int prev_grad = sign_(vs[1].length() - vs[0].length());
	vector <double> gradient;
	for (int i = 2; i < vs.size(); ++i)
	{
		int cur_grad = sign_(vs[i].length() - vs[i - 1].length());
		if (cur_grad * prev_grad < 0)
		{
			vextreme = vec(vs[i - 1]);
			// ext_ind.push_back(i - 1);
			li.end = i - 1;
			idx.push_back(li);
			crd.push_back(vec(vs[li.begin].start, vs[li.end].start));
			li.begin = i - 1;
			gradient.push_back(vs[i - 1].length() - prev_level);

			prev_level = vs[i - 1].length();
			prev_grad = cur_grad;
		}
	}
	li.end = vs.size() - 1;
	idx.push_back(li);
	crd.push_back(vec(vs[li.begin].start, vs[li.end].start));
}

void Vorticity::output_graphic_info(vector <int> &idx, vecs &pnc, vecs &itp, point begin, point end, int start_idx, int end_idx, double *apr_err)
{
	ofstream fnc, fpnc, fitp;
	fnc.open(subfilename(NORMAL_COMPONENT_FILENAME_PREFIX).c_str());
	fpnc.open(subfilename(PROJECTION_FILENAME_PREFIX).c_str());
	fitp.open(subfilename(INTERPOLATION_FILENAME_PREFIX).c_str());

	Line ct(begin, end);
	*apr_err = 0.0;
	int count = 0;

	for (int i = 0; i < pnc.size(); ++i)
		if (ct.contains(pnc[i].start))
		{
			fnc << ct.normal_component(mvn[idx[i]].v).toGlanceFormat(NORM_COMP_COLOR);
			fpnc << pnc[i].toGlanceFormat(NORM_COMP_PROJ_COLOR);
			*apr_err += mvn[idx[i]].error;
			++count;
		}

	*apr_err /= count;

	for (int i = start_idx; i <= end_idx; ++i)
		fitp << itp[i].toGlanceFormat(INTERPOLATION_COLOR);

	fnc << ct.interval().toGlanceFormat(CUT_COLOR);
	fpnc << ct.interval().toGlanceFormat(CUT_COLOR);
	fitp << ct.interval().toGlanceFormat(CUT_COLOR);

	fnc.close();
	fpnc.close();
	fitp.close();
}

Vorticity::Vorticity(vector <movement> &m) : calc_index(0), sub_index(0), mvn(m)
{
	flog.open("log.txt");
}

Vorticity::~Vorticity()
{
	flog.close();
}

void Vorticity::set_cut(scut _cut)
{
	cut = _cut;
	if (cut.width == -1) cut.width = BREADTH_SEARCH_RADIUS;
}

/*void Vorticity::calc_prnd()
{
	point m = fl_segment.middle();
	Line fl_line(fl_segment);
	fl_prnd = fl_line.perpendicular(m, BREADTH_SEARCH_RADIUS);
}*/

void Vorticity::calculate(vector <vorticity_result> &result)
{
	++calc_index;
	sub_index = 0;

	Line flp_line(cut.v());

	vecs vs;
	vector <int> idx;

	Border border(cut.v());

	flog << "-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-   Cut #" << calc_index << "   -|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-\n";

	ofstream fnc;
	ofstream fpnc;
	fnc.open(filename(NORMAL_COMPONENT_FILENAME_PREFIX).c_str());
	fpnc.open(filename(PROJECTION_FILENAME_PREFIX).c_str());

	double sum_S = 0.0, sum_V = 0.0;

	vector <string> slog;

	for (int i = 0; i < mvn.size(); ++i)
	{
		if (/*fls_line.contains_projection(mvn[i].v.start) &&*/ flp_line.contains_projection(mvn[i].v.start) && 
			(flp_line.distance_to(mvn[i].v.start) < deg2LS(cut.width)) && mvn[i].velocity > 0.00001)
		{
			vec nc_proj = flp_line.normal_component_projection(mvn[i].v);
			vs.push_back(nc_proj);
			idx.push_back(i);

			fnc << flp_line.normal_component(mvn[i].v).toGlanceFormat(NORM_COMP_COLOR);
			fpnc << nc_proj.toGlanceFormat(NORM_COMP_PROJ_COLOR);

			border.refresh(nc_proj.start);

			slog.push_back(mvn[i].toString());

			sum_S += mvn[i].v.length_m();
			sum_V += mvn[i].velocity;
		}
	}

	fnc << border.get_border().toGlanceFormat(CUT_COLOR);
	fpnc << border.get_border().toGlanceFormat(CUT_COLOR);
	fnc.close();
	fpnc.close();
	

	if (vs.size() == 0)
	{
		flog << "Favorable vectors are not found!\n\n";
		vorticity_result vort_res(flp_line.interval(), cut.width);
		return result.push_back(vort_res);	
	} 

	Line cut_line(border.get_border());

	flog << "Cut coordinates  are\t" << vec::header() << endl;
	flog << "\t\t\t" << cut_line.interval().toStr() << endl;
	flog << "    coefficients are\t" << cut_line.toString() << endl;
	flog << vs.size() << " vectors are found:\n";
	flog << movement::header() << endl;
	copy(slog.begin(), slog.end(), ostream_iterator<string>(flog, "\n"));

	double *itp_step_size = new double(0.0);
	Interpolation itp(border.get_border(), vs, calc_index, cut.itp_interval, cut.weight_coef);
	vecs ivs = itp.take(vs.size() * 5, itp_step_size);

	flog << "As a result of the interpolation " << ivs.size() << " vectors are received.\n";

	ofstream fitp;
	fitp.open(filename(INTERPOLATION_FILENAME_PREFIX).c_str());
	for (int i = 0; i < ivs.size(); ++i)
		fitp << ivs[i].toGlanceFormat(INTERPOLATION_COLOR);
	fitp << border.get_border().toGlanceFormat(CUT_COLOR);
	fitp.close();


	vector <line_index> ext_ind;
	vecs lc;
	find_extremum(ivs, ext_ind, lc);
	if (lc.size() > 0)
	{
		lc[0].start = border.get_border().start;
		lc[lc.size() - 1].end = border.get_border().end;
		if (lc.size() > 1)
		{
			flog << lc.size() - 1 << " extremum " << ((lc.size() > 2) ? "are" : "is") << " found:\n";
			flog << point::header() << endl;
			for (int i = 1; i < lc.size(); ++i)
				flog << lc[i].start.toStr() << endl;

		}
		
	}

	ofstream fal;
	fal.open(filename(APPROX_LINE_PREFIX).c_str());
	fal << border.get_border().toGlanceFormat(CUT_COLOR);

	// double K = mvn[idx[0]].v.length_m() / mvn[idx[0]].velocity;
	double T = sum_S / sum_V;

	flog << "T = " << T << " seconds\n\n";

	Least_Squares_Method lsm_orig(vs);
	Least_Squares_Method lsm_itp(ivs);

	for (int i = 0; i < ext_ind.size(); ++i)
	{
		++sub_index;

		flog << "------------------------Piece #" << sub_index << "---------------------------------------\n";

		double *lsm_err = new double(0.0);
		double *lsm_avr = new double(0.0);
		double *lsm_msd = new double(0.0);
		double *alpha = new double(0.0);
		double *alpha_left = new double(0.0);
		double *alpha_right = new double (0.0);
		flog << "Regression line:\talpha\t\tleft\tright\n";
		Line l = lsm_itp.calculate(lc[i], lsm_err, lsm_avr, lsm_msd, alpha, alpha_left, alpha_right);
		flog << "\tinterpolation\t" << *alpha << " " << *alpha_left << " " << *alpha_right << endl;
		// Line l = lsm_orig.calculate(lc[i], lsm_err, lsm_avr, lsm_msd, alpha, alpha_left, alpha_right);
		// flog << "\toriginal\t" << *alpha << " " << *alpha_left << " " << *alpha_right << endl;

		flog << "Regression line coefficients are " << l.toString() << endl;

		Line gl(vec(ivs[ext_ind[i].begin].start, ivs[ext_ind[i].end].start));
		l.set_bounds(gl);
		// li.set_bounds(gl);

		flog << "After setting bounds regression line\n" <<
				"  coordinates  are\t" << vec::header() << endl <<
				"\t\t\t" << l.interval().toStr() << endl /*<<
				"  coefficients are\t" << l.toString() << endl << endl*/;

		point p1 = gl.projection(l.start());
		double d1 = l.start().distance_to_m(p1);
		point p2 = gl.projection(l.end());
		double d2 = l.end().distance_to_m(p2);

		double V1 = d1 / T;
		double V2 = d2 / T;

		flog << "Vorticity calculate\n";
		flog << "V1 = " << V1 << " V2 = " << V2 << endl;		

		double dn = p1.distance_to_m(p2);
		flog << "dn = " << dn << endl;

		double Vn = - fabs(V1 - V2) / dn;
		flog << "-dV/dn = " << Vn << endl << endl;

		ofstream fsubal;
		fsubal.open(subfilename(APPROX_LINE_PREFIX).c_str());
		fsubal << l.interval().toGlanceFormat(APPROX_LINE_COLOR, 3, NO_ANCHOR);
		// fal << li.interval().toGlanceFormat(INTERPOLATION_COLOR, 3, "NO_ANCHOR");
		fsubal.close();

		fal << l.interval().toGlanceFormat(APPROX_LINE_COLOR, 3, NO_ANCHOR);

		double angle = l.angle(cut.v());
		// double diff = tan(angle);
		double diff = angle;
		
		double *itp_acc = new double(0.0);
		double *msd = new double(0.0);
		int *point_count = new int(0);
		itp.calc_accuracy(itp_acc, msd, point_count, lc[i].start, lc[i].end);
		if (*point_count < 3)
		{
			result.push_back(vorticity_result(vec(lc[i].start, lc[i].end), cut.width, itp.get_interval(), itp.get_weight_coef(), *point_count));
			continue;
		}

		double *apr_err = new double(0.0);
		output_graphic_info(idx, vs, ivs, lc[i].start, lc[i].end, ext_ind[i].begin, ext_ind[i].end, apr_err);

		int itp_step_count = ext_ind[i].end - ext_ind[i].begin;

		vorticity_result vort_res(vec(lc[i].start, lc[i].end), cut.width, itp.get_interval(), itp.get_weight_coef(), *point_count, *apr_err, *itp_acc, *msd, diff, 
									Vn, *alpha, *alpha_left, *alpha_right, *lsm_err, *lsm_avr, *lsm_msd, itp_step_count, *itp_step_size);
		
		result.push_back(vort_res);
	}

	fal.close();

}

#endif // VORTICITY_H