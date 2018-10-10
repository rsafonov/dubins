#include <dubins/dubinscar.h>

const char* PathTypeName[11] = {"LSL", "LSR", "RSL", "RSR", "RLR", "LRL", "TSL", "TSR", "LST", "RST", "TST"};
const double M_2PI = 2*M_PI;

inline double mod2pi(double theta)
{
	return theta - M_2PI*floor(theta/M_2PI);    
}

DubinsCar::DubinsCar()
{
	cache_created = false;
	debug = false;
	nn = 0;
	len_allocated = 0;
	h = 0.0;
	len = 0;
	pathlen = 0.0;
	rho = 0.0; 
	theta = 0.0;
	alpha = 0.0;
	beta = 0.0;
	sa = 0.0;
	sb = 0.0;
	ca = 0.0;
	cb = 0.0;
	cab = 0.0;
	dist = 0.0;
	dist2= 0.0;
	required_type = -1;

	traj_segment_length[0] = 0;
	traj_segment_length[1] = 0;
	traj_segment_length[2] = 0;

	traj = NULL;
	sp = NULL;
	cp = NULL;
}

DubinsCar::~DubinsCar()
{
	if (traj != NULL) delete [] traj;
	if (sp != NULL) delete [] sp;
	if (cp != NULL) delete [] cp;
}

void DubinsCar::CreateCache()
{
	nn = rho*dconfig.number_of_angles;
	h = M_2PI/nn; 
	if(debug) printf("Dubins: rho=%f number_of_angles=%d nn=%d h=%f\n", rho, dconfig.number_of_angles, nn, h);

	cp = new double[nn];
	sp = new double[nn];

	double tt;
	//Pre-calculate sin/cos arrays to speed calculations
	for (int j=0; j<nn; j++)
	{
		tt = j*h;
		sp[j] = sin(tt);
		cp[j] = cos(tt);
		//SBPL_PRINTF("i %d t %f sin %f cos %f\n", i, t, sp[i], cp[i]);
	}
	cache_created = true;
}

bool DubinsCar::Init(DubinsConfig _dconfig)
{
	dconfig.vel_tr = _dconfig.vel_tr;
	dconfig.vel_ang = _dconfig.vel_ang;
	dconfig.vel_turn = _dconfig.vel_turn;
	dconfig.number_of_angles = _dconfig.number_of_angles;
	dconfig.max_path_angle = _dconfig.max_path_angle;

	rho = _dconfig.vel_tr/_dconfig.vel_ang;
	dconfig.pathfile = _dconfig.pathfile;
	dconfig.pointfile = _dconfig.pointfile;

	if (debug)
	{
		printf("\nDubinsCar::Init start.\n");
		printf("vel_tr %f vel_ang %f rho = %f number_of_angles=%d max_path_angle=%f\n", dconfig.vel_tr, dconfig.vel_ang, rho, dconfig.number_of_angles, dconfig.max_path_angle);
	}
	return true;
}

/*
 Find the Dubins path given initial coordinates (start.x, start.y) and direction (start.yaw),  
 final coordinates (end.x, end.y) and direction (end.yaw) 
 and minimal turning radius rho = dconfig.vel_tr/dconfig.vel_ang.
 GetBestPath returns dpath structure populated with dubins path parameters
 */

int DubinsCar::GetBestPath(Pose2D start, Pose2D end)
{
	double dx, dy, tmp1, tmp2, s1, s2, s3, cost, w1, w2, w3, w4, w5, pangle; 
	
	if (debug) printf("\nGetBestPath start.\n");

	//Cleanup dpath structure
	dpath.type = -1; dpath.cost = INFINITY;
	dpath.s1 = 0.0; dpath.s2 = 0.0; dpath.s3 = 0.0;
	dpath.gamma1 = 0.0; dpath.gamma3 = 0.0;
	dpath.sgn1 = 1; dpath.sgn3 = 1;
	
	start_ = start;
	end_ = end;

	//Real distance between start and end
	dx = end.x - start.x;
	dy = end.y - start.y;
	dist = sqrt(dx*dx + dy*dy)/rho;

	if(debug) printf("start: %f %f end %f %f\n", start.x, start.y, end.x, end.y);
	if (debug) printf("dx=%f dy=%f dist=%f\n", dx, dy, dist);

	//Real angle between x-axis 
	theta = mod2pi(atan2(dy, dx));
	//Normalized directions at the start and end points 
	alpha = mod2pi(mod2pi(start.yaw) - theta);
	beta  = mod2pi(mod2pi(end.yaw) - theta);
	if (debug) printf("theta=%f alpha=%f beta=%f\n", theta, alpha, beta);

	//Some values used more then once in calculations
	sa = sin(alpha);            
	sb = sin(beta);             
	ca = cos(alpha);           
	cb = cos(beta);             
	cab = 2*cos(alpha - beta);
	dist2 = dist*dist;
	w1 = 2*dist*(sa - sb);
	w2 = 2*dist*(sa + sb);
	w3 = 6.0 - dist2 + cab;
	w4 = dist + sa - sb;
	w5 = dist - sa + sb;

	//Determine what path type is requested
	if (required_type >= 0)
	{
		switch (required_type)
		{
			case LSL:
				goto LSL_;
			case LSR:
				goto LSR_;
			case RSL:
				goto RSL_;
			case RSR:
				goto RSR_;
			case RLR:
				goto RLR_;
			case LRL:
				goto LRL_;
			default:
				return dpath.type;
		}
	}

LSL_: //0
	tmp2 = 2 + dist2 - cab + w1;
	if( tmp2 >= 0 && w4 >= 0) 
	{
		tmp1 = atan((cb - ca)/w4);  
		s1 = mod2pi(-alpha + tmp1);
		s2 = sqrt(tmp2);
		s3 = mod2pi(beta - tmp1);
		cost = mincost(s1, s2, s3, LSL, dconfig.vel_tr, alpha, beta, dpath);
		pangle = path_angle(LSL, s1, s2, s3);
		if (debug) printf("LSL (%d): t=%f p=%f q=%f cost=%f pangle=%f\n", LSL, s1, s2, s3, cost, pangle);
	}
	if (required_type >= 0) return LSL;


LSR_: //1
	tmp1 = -2 + dist2 + cab + w2;
	if( tmp1 >= 0 ) 
	{
		s2 = sqrt(tmp1);
		tmp2 = atan((-ca - cb)/(dist + sa + sb)) - atan(-2.0/s2);
		s1 = mod2pi(-alpha + tmp2);
		s3 = mod2pi(-mod2pi(beta) + tmp2);
		cost = mincost(s1, s2, s3, LSR, dconfig.vel_tr, alpha, beta, dpath);
		pangle = path_angle(LSR, s1, s2, s3);
		if (debug) printf("LSR (%d): t=%f p=%f q=%f cost=%f pangle=%f\n", LSR, s1, s2, s3, cost, pangle);	
	}
	if (required_type >= 0) return LSR;

RSL_: //2
	tmp1 = dist2 - 2 + cab - w2;
	if( tmp1 > 0 ) 
	{
		s2 = sqrt(tmp1);
		tmp2 = atan((ca + cb)/(dist - sa - sb)) - atan(2.0/s2);
		s1 = mod2pi(alpha - tmp2);
		s3 = mod2pi(beta - tmp2);
		cost = mincost(s1, s2, s3, RSL, dconfig.vel_tr, alpha, beta, dpath);
		pangle = path_angle(RSL, s1, s2, s3);
		if (debug) printf("RSL (%d): t=%f p=%f q=%f cost=%f pangle=%f\n", RSL, s1, s2, s3, cost, pangle);
	}
	if (required_type >= 0) return RSL;
    
RSR_: //3
	tmp2 = 2 + dist2 - cab - w1;
	if( tmp2 >= 0 && w5 >= 0) 
	{
		tmp1 = atan((ca - cb)/w5);
		s1 = mod2pi(alpha - tmp1);
		s2 = sqrt(tmp2);
		s3 = mod2pi(-beta + tmp1);
		cost = mincost(s1, s2, s3, RSR, dconfig.vel_tr, alpha, beta, dpath);
		pangle = path_angle(RSR, s1, s2, s3);
		if (debug) printf("RSR (%d): t=%f p=%f q=%f cost=%f pangle=%f\n", RSR, s1, s2, s3, cost, pangle);
	}
	if (required_type >= 0) return RSR;
RLR_: //4
	tmp1 = (w3 + w1)/8.0;
	if( fabs(tmp1) < 1) 
	{
		s2 = acos(tmp1);
		s1 = mod2pi(alpha - atan2(ca - cb, w5) + s2/2.0);
		s3 = mod2pi(alpha - beta - s1 + s2);
		cost = mincost(s1, s2, s3, RLR, dconfig.vel_tr, alpha, beta, dpath);
		pangle = path_angle(RLR, s1, s2, s3);
		if (debug) printf("RLR (%d): t=%f p=%f q=%f cost=%f pangle=%f\n", RLR, s1, s2, s3, cost, pangle);
	}
	if (required_type >= 0) return RLR;

LRL_: //5
	tmp1 = (w3 - w1)/8.0;
	if( fabs(tmp1) < 1) 
	{
		s2 = acos(tmp1);
		s1 = mod2pi(-alpha - atan2(ca - cb, w4) + s2/2.0);
		s3 = mod2pi(beta - alpha - s1 + s2);
		cost = mincost(s1, s2, s3, LRL, dconfig.vel_tr, alpha, beta, dpath);
		pangle = path_angle(LRL, s1, s2, s3);
		if (debug) printf("LRL (%d): t=%f p=%f q=%f cost=%f pangle=%f\n", LRL, s1, s2, s3, cost, pangle);
	}
	if (required_type >= 0) return LRL;

	if (debug) printf("dpath.type = %d\n", dpath.type);
	return dpath.type;
}

int DubinsCar::GetBestPathEx(Pose2D start, Pose2D end)
{
	int sgn1, sgn3;  
	double gamma1, gamma3, tmp1, s1, s2, s3, dsb, dsa, cost; 
	
	if (debug) printf("\nGetBestPathEx start.\n");

	int best_type = GetBestPath(start, end);

	dsb = 2*dist*sb;
	dsa = 2*dist*sa;

	if (required_type >= 0)
	{
		switch (required_type)
		{
			case TSL:
				goto TSL_;
			case TSR:
				goto TSR_;
			case LST:
				goto LST_;
			case RST:
				goto RST_;
			case TST:
				goto TST_;
			default:
				return best_type;
		}
	}

TSL_: //6
	tmp1 = dist2 - dsb;
	if( tmp1 > 0) 
	{
		s2 = sqrt(tmp1);
		gamma1 = mod2pi(asin(-(dist - sb - s2*cb)/(s2*s2 + 1)));
		sgn1 = turn_in_place(alpha, gamma1, s1);
		s3 = M_2PI - gamma1 + beta;
		cost = mincost_ex(s1, s2, s3, TSL, gamma1, 0.0, sgn1, 1, dconfig.vel_tr, dconfig.vel_turn, alpha, beta, dpath);
		if (debug) printf("TSL (%d): t=%f p=%f q=%f cost=%f\n", TSL, s1, s2, s3, cost);
	}
	if (required_type >= 0) return TSL;
    
TSR_: //7
	tmp1 = dist2 + dsb;
	if( tmp1 > 0) 
	{
		s2 = sqrt(tmp1);
		gamma1 = mod2pi(asin((dist + sb - s2*cb)/(s2*s2 + 1)));
		sgn1 = turn_in_place(alpha, gamma1, s1);
		s3 = M_2PI + gamma1 - beta;
		cost = mincost_ex(s1, s2, s3, TSR, gamma1, 0.0, sgn1, 1, dconfig.vel_tr, dconfig.vel_turn, alpha, beta, dpath);
		if (debug) printf("TSR (%d): t=%f p=%f q=%f cost=%f\n", TSR, s1, s2, s3, cost);
	}
	if (required_type >= 0) return TSR;

LST_: //8
	tmp1 = dist2 + dsa;
	if( tmp1 > 0) 
	{
		s2 = sqrt(tmp1);
		gamma3 = mod2pi(asin((dist + sa - s2*ca)/(s2*s2 + 1)));
		sgn3 = turn_in_place(gamma3, beta, s3);
		s1 = M_2PI + gamma3 - alpha;
		cost = mincost_ex(s1, s2, s3, LST, 0.0, gamma3, 1, sgn3, dconfig.vel_tr, dconfig.vel_turn, alpha, beta, dpath);
		if (debug) printf("LST (%d): t=%f p=%f q=%f cost=%f\n", LST, s1, s2, s3, cost);
	}
	if (required_type >= 0) return LST;
    
RST_: //9
	tmp1 = dist2 - dsa;
	if( tmp1 > 0) 
	{
		s2 = sqrt(tmp1);
		gamma3 = mod2pi(asin(-(dist - sa - s2*ca)/(s2*s2 + 1)));
		sgn3 = turn_in_place(gamma3, beta, s3);
		s1 = M_2PI - gamma3 + alpha;
		cost = mincost_ex(s1, s2, s3, RST, 0.0, gamma3, 1, sgn3, dconfig.vel_tr, dconfig.vel_turn, alpha, beta, dpath);
		if (debug) printf("RST (%d): t=%f p=%f q=%f cost=%f\n", RST, s1, s2, s3, cost);
	}
	if (required_type >= 0) return RST;

TST_: //10

	s2 = dist;		
	sgn1 = turn_in_place(alpha, 0.0, s1);
	sgn3 = turn_in_place(0.0, beta, s3);
	cost = mincost_ex(s1, s2, s3, TST, 0.0, 0.0, sgn1, sgn3, dconfig.vel_tr, dconfig.vel_turn, alpha, beta, dpath);
	if (debug) printf("TST (%d): t=%f p=%f q=%f cost=%f\n", TST, s1, s2, s3, cost);
	if (required_type >= 0) return TST;

	if (debug) printf("dpath.type = %d\n", dpath.type);
	return dpath.type;
}

bool DubinsCar::BuildTrajectory()
{
	int i = 0, j, dir_type = 0, k = 0, k0, k1, k2, flag = 0, sgn = 0;
	double t = 0.0, sn0, cs0, sn1, cs1, sn2, cs2, d, pp, delta; 
	double ss = 0.0, cc = 0.0, ss0 = 0.0, cc0 = 0.0;
	double q10x, q10y, q10t;		//Normalized coordinates of the last point of the first segment 
	double q20x, q20y, q20t;		//Normalized coordinates of the last point of the second segment 
	double qtx = 0.0, qty = 0.0, qtt = 0.0;
	double qqx, qqy, qqt;
	double qt0x = 0.0, qt0y = 0.0, qt0t = 0.0;
	double qix = 0.0, qiy = 0.0, qit = alpha; //Normalized coordinates and orientation of the start point
	const int* types = DIRDATA[dpath.type];

	double sin_theta = sin(theta);
	double cos_theta = cos(theta);

	if (!cache_created) CreateCache();
	
	pathlen = dpath.s1 + dpath.s2 + dpath.s3;	
	len= pathlen/h + 3;
	if (len > len_allocated)
	{
		if (debug) printf("\nReallocating memory for traj (len=%d len_allocated=%d).\n", len, len_allocated);
		if (traj != NULL) delete [] traj;		
		len_allocated = 2*len;
		traj = new Pose2D[len_allocated];
	}
	else
	{
		for (j=0; j<len_allocated; j++)
		{
			traj[j].x = 0.0; traj[j].y = 0.0; traj[i].yaw = 0.0;
		}
	}
	
	if (debug)
	{
		printf("\nBuildTrajectory start.\n");
		printf("s1=%f s2=%f s3=%f\n", dpath.s1, dpath.s2, dpath.s3);
		printf("types: %d %d %d\n", types[0], types[1], types[2]);
 		printf("pathlen=%f number_of_angles=%d nn=%d h=%f len=%d\n", pathlen, dconfig.number_of_angles, nn, h, len);
		fflush(stdout);
	}

	sn0 = sin(qit);
	cs0 = cos(qit);
	k0 = qit/h;
	if (debug) printf("qit=%f sn0=%f cs0=%f k0=%d\n", qit, sn0, cs0, k0);
	
	//Calculate normalized coordinates of the last point of the first segment
	if (types[0] == L_SEG) 
	{
		ss = sin(qit + dpath.s1);
		cc = cos(qit + dpath.s1);
		q10x = qix + ss - sn0;
		q10y = qiy - cc + cs0;
		q10t = mod2pi(qit + dpath.s1);
	}
	else if (types[0] == R_SEG) 
	{
		ss = sin(qit - dpath.s1);
		cc = cos(qit - dpath.s1);
		q10x = qix - ss + sn0;
		q10y = qiy + cc - cs0;
		q10t = mod2pi(qit - dpath.s1);
	}
	else if (types[0] == T_SEG) 
	{
		q10x = qix;
		q10y = qiy;
		q10t = mod2pi(qit + dpath.s1*dpath.sgn1); 
	}
	else
	{
		printf("First segment is invalid: types[0] = %d\n", types[0]);
		fflush(stdout);
		return false;
	}
	if (debug) printf("q10=%f %f %f\n", q10x, q10y, q10t);

	//Calculate real coordinates of the last point of the first segment	
	q1.x = (q10x*cos_theta - q10y*sin_theta)*rho + start_.x;
	q1.y = (q10x*sin_theta + q10y*cos_theta)*rho + start_.y;
	q1.yaw = mod2pi(q10t + theta);
	if (debug) printf("q1=%f %f %f\n", q1.x, q1.y, q1.yaw);
    
	sn1 = sin(q10t);
	cs1 = cos(q10t);
	k1 = q10t/h;
	if (debug) printf("q10t=%f sn1=%f cs1=%f k1=%d\n", q10t, sn1, cs1, k1);

	//Calculate normalized coordinates of the last point of the second segment
	if (types[1] == L_SEG) 
	{
		ss = sin(q10t + dpath.s2);
		cc = cos(q10t + dpath.s2);
		q20x = q10x + ss - sn1;
		q20y = q10y - cc + cs1;
		q20t = mod2pi(q10t + dpath.s2);
	}
	else if (types[1] == R_SEG) 
	{
		ss = sin(q10t - dpath.s2);
		cc = cos(q10t - dpath.s2);
		q20x = q10x - ss + sn1;
		q20y = q10y + cc - cs1;
		q20t = mod2pi(q10t - dpath.s2);
	}
	else if (types[1] == S_SEG) 
	{
		q20x = q10x + cs1*dpath.s2;
		q20y = q10y + sn1*dpath.s2;
		q20t = mod2pi(q10t);
	}
	else
	{
		printf("Second segment is invalid: types[1] = %d\n", types[1]);
		fflush(stdout);
		return false;
	}
	if (debug) printf("q20=%f %f %f\n", q20x, q20y, q20t);

	//Calculate real coordinates of the last point of the second segment
	q2.x = (q20x*cos_theta - q20y*sin_theta)*rho + start_.x;
	q2.y = (q20x*sin_theta + q20y*cos_theta)*rho + start_.y;
	q2.yaw = mod2pi(q20t + theta);
	if (debug) printf("q2=%f %f %f\n", q2.x, q2.y, q2.yaw);

	sn2 = sin(q20t);
	cs2 = cos(q20t);
	k2 = q20t/h;
	if (debug) 
	{
		printf("q20t=%f sn2=%f cs2=%f k2=%d\n", q20t, sn2, cs2, k2);
		printf("k0=%d k1=%d k2=%d\n", k0, k1, k2);
		printf("types[0]=%d\n", types[0]);
		fflush(stdout);
	}
    
	pp = dpath.s1 + dpath.s2;

	traj[i].x = start_.x;
	traj[i].y = start_.y;
	traj[i].yaw = start_.yaw;
	i++;
	traj_segment_length[0]++;
 
	while (t < pathlen)
	{
		if (t < dpath.s1)
		{
			ss0 = sn0; 
			cc0 = cs0;
			qqx = qix;
			qqy = qiy;
			qqt = qit;
			dir_type = types[0];
	    		
			if (flag == 0)
			{
				if (dir_type == L_SEG || dir_type == R_SEG) 
				{
					delta = abs(k0*h - qqt);			
					t = delta;
					k = k0;
				}
				else if (dir_type == T_SEG) 
				{
					delta = 0.0;
					t = 0.0;
					k = k0;
					sgn = dpath.sgn1;
				}
				flag = 1;
				if (debug) printf("\ni=%d dir_type=%d k0=%d k=%d qqt=%f t=%f k0*h=%f delta = %f\n", i, dir_type, k0, k, qqt, t, k0*h, delta);
			}
			d = t;
			traj_segment_length[0]++;
		}
		else if (t < pp)
		{
			ss0 = sn1; 
			cc0 = cs1;
			qqx = q10x;
			qqy = q10y;
			qqt = q10t;
			dir_type = types[1];

			if (flag == 1)
			{
				if (dir_type == L_SEG || dir_type == R_SEG) 
				{
					delta = abs(k1*h - qqt);			
					t = dpath.s1 + delta;
					k = k1;
				}
				else if (dir_type == T_SEG || dir_type == S_SEG) 
				{
					delta = 0.0;
					t = dpath.s1;
					k = k1;
				}
				flag = 0;
				if (debug) printf("\ni=%d dir_type=%d k1=%d k=%d qqt=%f t=%f k1*h=%f delta = %f\n", i, dir_type, k1, k, qqt, t, k1*h, delta);
			}
			d = t - dpath.s1;
			traj_segment_length[1]++;
		}
		else
		{
			ss0 = sn2; 
			cc0 = cs2;
			qqx = q20x;
			qqy = q20y;
			qqt = q20t;
			dir_type = types[2];

			if (flag == 0)
			{
				if (dir_type == L_SEG  || dir_type == R_SEG) 
				{
					delta = abs((k2+1)*h - qqt);			
					t = pp + delta;
					k = k2;
				}
				else if (dir_type == T_SEG) 
				{
					delta = 0.0;
					t = pp;
					k = k2;
					sgn = dpath.sgn3;
				}
				flag = 1;
				if (debug) printf("\ni=%d dir_type=%d k2=%d k=%d qqt=%f t=%f k1*h=%f delta = %f pp=%f\n", i, dir_type, k2, k, qqt, t, k1*h, delta, pp);
			}
			d = t - pp;
			traj_segment_length[2]++;
		}
        
		switch (dir_type)
		{
			case L_SEG:
				k = k+1;
				if (k > nn-1) k = 0;
                     
				qt0x = qqx + sp[k] - ss0;
				qt0y = qqy - cp[k] + cc0;
				qt0t = mod2pi(qqt + d);
				break;

			case R_SEG:
				qt0x = qqx - sp[k] + ss0;
				qt0y = qqy + cp[k] - cc0;
				qt0t = mod2pi(qqt - d);

				k = k-1;
				if (k < 0) k = nn-1;            
				break;

			case S_SEG:
				qt0x = qqx + cc0 * d;
				qt0y = qqy + ss0 * d;
				qt0t = mod2pi(qqt);
				break;

			case T_SEG:
				qt0x = qqx;
				qt0y = qqy;
				qt0t = mod2pi(qqt + d*sgn); 
				break;
		} 

		qtx = (qt0x*cos_theta - qt0y*sin_theta)*rho + start_.x;
		qty = (qt0x*sin_theta + qt0y*cos_theta)*rho + start_.y;
		qtt = mod2pi(qt0t + theta);
		if (debug)
		{
			if (i/10*10 == i) printf("i=%d k=%d t=%f qt: %f %f %f dir_type=%d ss0=%f sp[k]=%f\n", i, k, t, qtx, qty, qtt, dir_type, ss0, sp[k]);
			////printf("k=%d t=%f qt: %f %f %f dir_type=%d\n", k, t, qtx, qty, qtt, dir_type);
			//fflush(stdout);
		}
       
		traj[i].x = qtx;
		traj[i].y = qty;
		traj[i].yaw = qtt;
	
		i = i+1;
		t = t+h; 
	}

	traj[i].x = end_.x;
	traj[i].y = end_.y;
	traj[i].yaw = end_.yaw;
	len = i+1;
	traj_segment_length[2]++;

	if (debug) 
	{
		printf("start: %f %f %f end: %f %f %f\n", start_.x, start_.y, start_.yaw, end_.x, end_.y, end_.yaw);
		printf("len=%d traj_segment_length: %d %d %d\n", len, traj_segment_length[0], traj_segment_length[1], traj_segment_length[2]);		
		printf("BuildTrajectory end.\n");
		fflush(stdout);
	}
	return true;
}

bool DubinsCar::ExportTrajectory(string suffix)
{
	int i;
	double px[4], py[4], pyaw[4];
	string file_name = "";

	//Write generated dubins path trajectory to a file
	file_name = GetExportFileName(dconfig.pathfile, suffix);
	FILE* ptr1 = fopen(file_name.c_str(),"w");
	if (ptr1 == NULL)
	{
		printf("Could not open file %s\n", file_name.c_str());
		fflush(stdout);
		return false;
	}
	//printf("Opened file %s\n", dconfig.pathfile.c_str());
	//fflush(stdout);

	fprintf(ptr1, "len=%d\n", len);

	for (i=0; i<len; i++)
	{
		fprintf(ptr1, "%f %f %f\n", traj[i].x, traj[i].y, traj[i].yaw);
	}
	fflush(stdout);
	fclose(ptr1); 

	//Write start/end segment points to a file
	px[0] = start_.x;
	px[1] = q1.x;
	px[2] = q2.x;
	px[3] = end_.x;

	py[0] = start_.y;
	py[1] = q1.y;
	py[2] = q2.y;
	py[3] = end_.y;

	pyaw[0] = start_.yaw;
	pyaw[1] = q1.yaw;
	pyaw[2] = q2.yaw;
	pyaw[3] = end_.yaw;

	file_name = GetExportFileName(dconfig.pointfile, suffix);
	FILE* ptr = fopen(file_name.c_str(),"w");
	if (ptr == NULL)
	{
		printf("Could not open file %s\n", file_name.c_str());
		fflush(stdout);
		return false;
	}

	//printf("Opened file %s\n", dconfig.pointfile.c_str());
	//fflush(stdout);	

	fprintf(ptr, "type=%d\ngamma1=%f\ngamma3=%f\nsgn1=%d\nsgn3=%d\ns1=%f\ns2=%f\ns3=%f\nH=%f\n", dpath.type, dpath.gamma1, dpath.gamma3, dpath.sgn1, dpath.sgn3, dpath.s1, dpath.s2, dpath.s3, h);

	for (i=0; i<4; i++)
	{
		fprintf(ptr, "%f %f %f\n", px[i], py[i], pyaw[i]);
	}
	fflush(stdout);
	fclose(ptr);     

	return true;
}

string DubinsCar::GetExportFileName(string fname, string suffix)
{
	string a = "";
	int len1 = fname.length();
	if (len1 == 0) return a;
	int len2 = suffix.length();
	if (len2 == 0) return fname;
	size_t pos = fname.find_last_of('.');
	if ((int)pos == -1) pos = len1;
	a = fname.substr(0, pos) + "_" + suffix + fname.substr(pos);
	if (debug) printf("a = %s\n", a.c_str());
	return a;
}
