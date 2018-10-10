#ifndef DUBINSCAR_H_
#define DUBINSCAR_H_

#include <cmath>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <dirent.h>
#include <sys/stat.h>
#include <errno.h>

using namespace std; 

// L stands for  left arc turn,  R - right arc turn, S - straight line, T - turn in place

#define LSL (0)
#define LSR (1)
#define RSL (2)
#define RSR (3)
#define RLR (4)
#define LRL (5)
#define TSL (6)
#define TSR (7)
#define LST (8)
#define RST (9)
#define TST (10)

#define SZ (3)

#define L_SEG (0)
#define S_SEG (1)
#define R_SEG (2)
#define T_SEG (3)

//const char* PathTypeName[11] = {"LSL", "LSR", "RSL", "RSR","RLR", "LRL", "TSL", "TSR", "LST", "RST", "TST"};
extern const char* PathTypeName[11];

//const double M_2PI = 2*M_PI;

const int DIRDATA[][3] = {
	{ L_SEG, S_SEG, L_SEG },
	{ L_SEG, S_SEG, R_SEG },
	{ R_SEG, S_SEG, L_SEG },
	{ R_SEG, S_SEG, R_SEG },
	{ R_SEG, L_SEG, R_SEG },
	{ L_SEG, R_SEG, L_SEG },
	{ T_SEG, S_SEG, L_SEG },
	{ T_SEG, S_SEG, R_SEG },
	{ L_SEG, S_SEG, T_SEG },
	{ R_SEG, S_SEG, T_SEG },
	{ T_SEG, S_SEG, T_SEG }
};

struct Pose2D
{	
	Pose2D(double _x, double _y, double _yaw): x(_x), y(_y), yaw(_yaw){}

	Pose2D()
	{
		x = 0.0;
		y = 0.0;
		yaw = 0.0;
	}
	
	Pose2D &operator=(const Pose2D &rhs)
	{
		if (this == &rhs)
			return *this;
		else
		{
			memcpy(this, &rhs, sizeof(Pose2D));
			return *this;
		}
	}

	double x; 
	double y;
	double yaw;
};

struct DubinsPath
{
	DubinsPath()
	{
		type = -1;
		s1 = 0.0;
		s2 = 0.0;
		s3 = 0.0;
		cost = INFINITY;
		gamma1 = 0.0;
		gamma3 = 0.0;
		sgn1 = 1;
		sgn3 = 1;
		path_angle = 0.0;
	}

	int type;		//type of the dubins path
	double s1;		//length of the first segment
	double s2;		//length of the second segment
	double s3;		//length of the third segment
	double cost;		//cost of the dubins path
	double gamma1;		//turn in place angle for the first segment
	double gamma3;		//turn in place angle for the last segment
	int sgn1;		//direction of the turn in place for the first segment
	int sgn3;		//direction of the turn in place for the last segment
	double path_angle;	//overall turning angle of the dubins path;
};

struct DubinsConfig
{
	DubinsConfig()
	{
		vel_tr = 0.0;
		vel_ang = 0.0;
		vel_turn = 0.0;
		number_of_angles = 100;
		pathfile = "";
		pointfile = "";
		outdir = "";
		max_path_angle = 0.0;
	}

	DubinsConfig(double _vel_tr, double _vel_ang, double _vel_turn, int _number_of_angles, string _pathfile, string _pointfile): 
			vel_tr(_vel_tr), vel_ang(_vel_ang), number_of_angles(_number_of_angles), pathfile(_pathfile), pointfile(_pointfile) {}

	double vel_tr;		//translational velocity
	double vel_ang;		//angular velocity along arc
	double vel_turn;	//angular velocity when turning in place
	int number_of_angles;	//number of angles - number of discrete angles between the start heading and the end heading when building a trajectory.
	string pathfile;	//file where dubins trajectory gets exported.
	string pointfile;	//file where start, end and intermideate poses of the trajectory get exported.
	string outdir;		//directory where files with results get exported;
	double max_path_angle;	//Maximum acceptable overall turning angle of the entire dubins path.
};

/*
inline double mod2pi(double theta)
{
	return theta - M_2PI*floor(theta/M_2PI);    
}
*/

inline int turn_in_place(double alpha, double gamma, double& turn_angle)
{
	int turn_sgn = 1;

	turn_angle = abs(gamma - alpha);
	if (turn_angle > M_PI) turn_angle = abs(2*M_PI - turn_angle); //M_2PI
                
	if (0 <= alpha && alpha < M_PI)
	{
		if (gamma <= alpha || gamma > alpha + M_PI) turn_sgn=-1;
	}
	else
	{
		if (alpha - M_PI <= gamma && gamma < alpha) turn_sgn=-1;
	}
	return turn_sgn;
}

inline double segment_time(double s, int dirtype, double vel_tr, double vel_turn)
{
	double t;
	if (dirtype == T_SEG)
		t = s/vel_turn;
	else
		t = s/vel_tr;
	return t;
}

inline double path_angle(int type, double s1, double s2, double s3)
{
	double path_angle = 0.0;
	
	const int* types = DIRDATA[type];	
	if (types[0] == L_SEG || types[0] == R_SEG || types[0] == T_SEG) 
		path_angle = path_angle + s1;
	if (types[1] == L_SEG || types[1] == R_SEG || types[1] == T_SEG) 
		path_angle = path_angle + s2;
	if (types[2] == L_SEG || types[2] == R_SEG || types[2] == T_SEG) 
		path_angle = path_angle + s3;
	//printf("path_angle=%f\n", path_angle);
	return path_angle;
}

inline double mincost(double s1, double s2, double s3, int type, double vel_tr, double alpha, double beta, DubinsPath& dpath)
{
	double cost = 	(s1 + s2 + s3)/vel_tr;
	if (cost < dpath.cost)
	{
		dpath.type = type;
		dpath.cost = cost;
		dpath.s1 = s1;
		dpath.s2 = s2;
		dpath.s3 = s3;
		dpath.path_angle = path_angle(type, s1, s2, s3);
	}
	return cost;
}

inline double mincost_ex(double s1, double s2, double s3, int type, double gamma1, double gamma3, int sgn1, int sgn3, 
			 double vel_tr, double vel_turn, double alpha, double beta, DubinsPath& dpath)
{
	double cost = 	segment_time(s1, DIRDATA[type][1], vel_tr, vel_turn) + 
			segment_time(s2, DIRDATA[type][2], vel_tr, vel_turn) + 
			segment_time(s3, DIRDATA[type][3], vel_tr, vel_turn);
	if (cost < dpath.cost)
	{
		dpath.type = type;
		dpath.cost = cost;
		dpath.s1 = s1;
		dpath.s2 = s2;
		dpath.s3 = s3;
		dpath.gamma1 = gamma1;
		dpath.sgn1 = sgn1;
		dpath.gamma3 = gamma3;
		dpath.sgn3 = sgn3;
		dpath.path_angle = path_angle(type, s1, s2, s3);
	}
	return cost;
}

class DubinsCar
{
public:
	// Constructor
	DubinsCar();

	// Destructor
	~DubinsCar();
	
	int GetBestPath(Pose2D start, Pose2D end);
	int GetBestPathEx(Pose2D start, Pose2D end);
	bool BuildTrajectory();
	bool ExportTrajectory(string suffix);
	bool Init(DubinsConfig _dconfig);

	DubinsConfig dconfig;
	DubinsPath dpath;
	Pose2D start_;
	Pose2D end_;

	Pose2D q1;
	Pose2D q2;

	Pose2D *traj;
	int traj_segment_length[3];

	bool debug;
	int required_type;
	int len;  //number of points in dubins the trajecory
	double rho; // forward velocity / angular velocity 

	double h; //parameter (time?) increment betwen sequential points on the dubins trajectory.
	double pathlen;	//total length of h dubins path

protected:
	string GetExportFileName(string fname, string suffix);

	void CreateCache();
	bool cache_created;
	int nn;
	int len_allocated;
	double alpha;		//Normalized start orientation 
	double beta;		//Normalized end orientation
	double theta;		//Original angle between the line connecting start and end and the x-axis
	double sa, sb, ca, cb, cab, dist, dist2;

	double *sp;
	double *cp;
};

#endif // DUBINSCAR_H
