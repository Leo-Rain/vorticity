#ifndef LSM_H
#define LSM_H

#include <stdio.h>
#include <iostream>
#include <fstream>

#include "geometry.h"
#include "ap.h"
#include "specialfunctions.h"


#define STUDENT_P 0.95

class Least_Squares_Method
{
	vecs vs;

public:
	Least_Squares_Method(vecs &_vs);

	Line calculate(vec itv, double *d, double *avr_dist, double *msd_dist, double *alpha, double *alpha_left, double *alpha_right);
};

Least_Squares_Method::Least_Squares_Method(vecs &_vs) : vs(_vs)
{

}

Line Least_Squares_Method::calculate(vec itv, double *d, double *avr_dist, double *msd_dist, double *alpha, double *alpha_left, double *alpha_right)
{
	vector <int> idx;
	double avr_x = 0.0, avr_y = 0.0;
	Line cut(itv);
	for (int i = 0; i <= vs.size(); ++i)
		if (cut.contains_projection(vs[i].start))
		{
			avr_x += vs[i].end.x;
			avr_y += vs[i].end.y;	
			idx.push_back(i);
		}

	int count = idx.size();

	avr_x /= count;
	avr_y /= count;

	double alpha_num = 0.0, alpha_den = 0.0, Syz = 0.0;
	for (int i = 0; i < count; ++i)
	{
		alpha_num += (vs[idx[i]].end.y - avr_y) * (vs[idx[i]].end.x - avr_x);
		alpha_den += (vs[idx[i]].end.x - avr_x) * (vs[idx[i]].end.x - avr_x);
		Syz += vs[idx[i]].end.y * (vs[idx[i]].end.x - avr_x);
	}
	*alpha = alpha_num / alpha_den;

	Line l(*alpha, -1, avr_y - *alpha * avr_x);

	double S2 = 0.0, SR2 = 0.0, sum_dist = 0.0;
	for (int i = 0; i < count; ++i)
	{
		S2 += (vs[idx[i]].end.y - avr_y) * (vs[idx[i]].end.y - avr_y);
		SR2 += (vs[idx[i]].end.y - l.x2y(vs[idx[i]].end.x)) * (vs[idx[i]].end.y - l.x2y(vs[idx[i]].end.x));
		sum_dist += l.distance_to(vs[idx[i]].end);
	}
	*d = 1 - (SR2 / S2);
	*avr_dist = sum_dist / count;
	double sum_dist2 = 0.0;
	for (int i = 0; i < count; ++i)
		sum_dist2 += (l.distance_to(vs[idx[i]].end) - *avr_dist) * (l.distance_to(vs[idx[i]].end) - *avr_dist);
	*msd_dist = sqrt(sum_dist2 / count);

	double t;
	if (count > 2)
		t =  alglib::invstudenttdistribution(count - 2, STUDENT_P);
	else
	{
		flog << "regression: number of points < 2\n";
		t = 3.5;
	}
	flog << "regression: student param \"t\" = " << t << endl;



	double msd_y = sqrt(SR2 / (count - 2));
	double Syz_z = Syz / alpha_den;
	double ts_Sz = t * msd_y / sqrt(alpha_den);

	*alpha_left = Syz_z - ts_Sz;
	*alpha_right = Syz_z + ts_Sz;

	// cout << "left " << Syz_z - ts_Sz << " right " << Syz_z + ts_Sz << endl;

	return l;
}

#endif // LSM_H