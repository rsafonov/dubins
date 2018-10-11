/*
 * Copyright (c) 2015, Maxim Likhachev
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the Carnegie Mellon University nor the names of its
 *       contributors may be used to endorse or promote products derived from
 *       this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */

#include <tinyxml.h>
#include <dubins/dubinscar.h>

extern bool ReadConfig(string ConfigXmlFile, Pose2D& start, Pose2D& end, DubinsConfig& dconfig, long& number_of_repeats, int& path_type, bool& debug);
extern string IntToStr(int i);
int DeleteFiles(string dir);

int main()
{
	double dt = 0.0, v = 0.0; 
	bool bres = false, debug;
	Pose2D start, end;
	DubinsConfig dconfig;
	long number_of_repeats;
	int path_type, countDeleted;
	time_t t0, t1, t2;
	string suffix = "";

	string config = "config/dubins_params.xml";
	
	if (ReadConfig(config, start, end, dconfig, number_of_repeats, path_type, debug))
	{ 
		if (number_of_repeats > 1)
		{
			t0 = time(0);
			printf("t0=%ld\n", t0);
		}

		DubinsCar dcar;
		dcar.debug = debug;
		dcar.required_type= path_type;
		dcar.Init(dconfig);
   
		//Calculate parameters of the optimal dubins path
		//int best_type = dcar.GetBestPathEx(start, end);
		int best_type = dcar.GetBestPath(start, end);
		printf("dubins_path returned: type=%s (%d) type=%d s1=%f s2=%f s3=%f cost=%f gamma1=%f sgn1=%d gamma3=%f sgn3=%d rho=%f path_angle=%f\n", 
			PathTypeName[best_type], best_type, dcar.dpath.type, dcar.dpath.s1, dcar.dpath.s2, dcar.dpath.s3, dcar.dpath.cost, 
			dcar.dpath.gamma1, dcar.dpath.sgn1, dcar.dpath.gamma3, dcar.dpath.sgn3, dcar.rho, dcar.dpath.path_angle);
		fflush(stdout);
		if (dcar.dpath.type == -1) return 0;
		if (dcar.dpath.path_angle > dconfig.max_path_angle)
		{
			printf("Dubins path angle exceeded the limit!\n");
			fflush(stdout);
			//return 0;
		}

		printf("outdir = %s\n", dconfig.outdir.c_str());
		if (dconfig.outdir == "")
		{
			size_t pos = dconfig.pathfile.find_last_of("/");
			if (pos != -1) dconfig.outdir = dconfig.pathfile.substr(0, pos+1);
		}
		printf("outdir = %s\n", dconfig.outdir.c_str());
		if (dconfig.outdir != "") countDeleted = DeleteFiles(dconfig.outdir);
	
		//Generate the dubins path trajectory    
		if (dconfig.number_of_angles > 0)
		{
			dcar.BuildTrajectory();
			printf("len = %d\n", dcar.len);
			//for (int i=0; i<dcar.len; i++)
			//{
			//	printf("%d %f %f %f\n", i, dcar.traj[i].x, dcar.traj[i].y, dcar.traj[i].yaw); 
			//}
			
			printf("Calling ExportTrajectory.\n");
			dcar.ExportTrajectory("");
		}

		//if number_of_repeats > 1, repeat dubins path calculation/generatrion multiple times to fund the average speed of calculations
		if (number_of_repeats > 1)
		{
			t1 = time(0);
			dt = difftime(t1, t0);
			printf("t1=%ld dt = %15.10f\n", t1, dt);

			for (int i=0; i<number_of_repeats; i++)
			{
				printf("*******i %d start **************\n", i);
				start.x = start.x + 10;
				start.y = start.y + 10;
				end.x = end.x - 10;
				end.y = end.y - 10;

				best_type = dcar.GetBestPath(start, end);
				//dcar.GetBestPathEx(start, end);
				printf("dubins_path returned: type=%s (%d) type=%d s1=%f s2=%f s3=%f cost=%f gamma1=%f sgn1=%d gamma3=%f sgn3=%d rho=%f\n", 
					PathTypeName[best_type], best_type, dcar.dpath.type, dcar.dpath.s1, dcar.dpath.s2, dcar.dpath.s3, dcar.dpath.cost, 
					dcar.dpath.gamma1, dcar.dpath.sgn1, dcar.dpath.gamma3, dcar.dpath.sgn3, dcar.rho);

				if (dconfig.number_of_angles > 0)
				{
					dcar.BuildTrajectory();
					suffix = IntToStr(i);
					dcar.ExportTrajectory(suffix);
				}
				printf("*******i %d end **************\n", i);
			}

			t2 = time(0);
			dt = difftime(t2, t1);
			v = dt/number_of_repeats;
			printf("t2=%ld dt = %15.10f v = %15.10f\n", t2, dt, v);
			fflush(stdout);
		}
	}
	return 0;    
}

bool ReadConfig(string ConfigXmlFile, Pose2D& start, Pose2D& end, DubinsConfig& dconfig, long& number_of_repeats, int& path_type, bool& debug)
{
	bool bload, bres = false;
	string stmp="", path_file = "path.txt", point_file = "points.txt", out_dir = "";
	const char *pText;
	double tmp = -1000000.0;
	int ntmp = -100000, error_count = 0;

	TiXmlElement *root, *element, *test_element;
	TiXmlDocument doc(ConfigXmlFile);
	TiXmlHandle hDoc(&doc), hRoot(0), hTestElement(0), hElement(0);

	printf("ReadConfig start for %s\n", ConfigXmlFile.c_str());
	fflush(stdout);

	number_of_repeats = 1;
	path_type = -1;
	debug = false;
	dconfig.max_path_angle = 0.0;
    
	bload = doc.LoadFile();
	if (bload)
	{
		//printf("Xml file loaded.\n");
		//fflush(stdout);

		root = hDoc.FirstChildElement("config").Element();
		if (root)
		{
			//printf("config element found.\n");
			//fflush(stdout);
			hRoot=TiXmlHandle(root); 
			element = hRoot.FirstChildElement("start").Element();
			if (element)
			{
				//printf("start element found.\n");
				//fflush(stdout);
				hElement = TiXmlHandle(element); 
				if (element->QueryDoubleAttribute("x", &start.x) != 0) error_count++;
				if (element->QueryDoubleAttribute("y", &start.y) != 0) error_count++;
				if (element->QueryDoubleAttribute("yaw", &start.yaw) != 0) error_count++;
				printf("start: x=%f y=%f yaw=%f\n", start.x, start.y, start.yaw);
			}
			else
			{
				printf("start element not found.\n");
				fflush(stdout);
				error_count++;
			}

			element = hRoot.FirstChildElement("end").Element();
			if (element)
			{
				hElement = TiXmlHandle(element); 
				if (element->QueryDoubleAttribute("x", &end.x) != 0) error_count++;
				if (element->QueryDoubleAttribute("y", &end.y) != 0) error_count++;
				if (element->QueryDoubleAttribute("yaw", &end.yaw) != 0) error_count++;
				printf("end: x=%f y=%f yaw=%f\n", end.x, end.y, end.yaw);
			}
			else
			{
				printf("end element not found.\n");
				fflush(stdout);
				error_count++;
			}

			element = hRoot.FirstChildElement("options").Element();
			if (element)
			{
				hElement = TiXmlHandle(element); 
				if (element->QueryDoubleAttribute("trans_velocity", &dconfig.vel_tr) != 0) error_count++;
				if (element->QueryDoubleAttribute("angular_velocity", &dconfig.vel_ang) != 0) error_count++;
				if (element->QueryDoubleAttribute("turn_in_place_velocity", &dconfig.vel_turn) != 0) error_count++;
				element->QueryIntAttribute("number_of_angles", &dconfig.number_of_angles);
				element->QueryStringAttribute("number_of_repeats", &stmp);
				if (stmp != "") number_of_repeats = atol(stmp.c_str());
				element->QueryIntAttribute("path_type", &path_type);
				element->QueryBoolAttribute("debug", &debug);
				if (element->QueryDoubleAttribute("max_path_angle", &dconfig.max_path_angle) != 0) error_count++;

				element = hElement.FirstChildElement("out_dir").Element();
				if (!element) error_count++;
				pText=element->GetText(); out_dir = pText;
				element = hElement.FirstChildElement("path_file").Element();
				if (!element) error_count++; 
				pText=element->GetText(); path_file = pText;
				element = hElement.FirstChildElement("point_file").Element();
				if (!element) error_count++; 
				pText=element->GetText(); point_file = pText;

				dconfig.pathfile = out_dir + path_file;
				dconfig.pointfile = out_dir + point_file;
				dconfig.outdir = out_dir;

				printf("vel_tr=%f vel_ang=%f vel_turn=%f number_of_angles=%d number_of_repeats=%ld path_type=%d debug=%d max_path_angle=%f\n", 
					dconfig.vel_tr, dconfig.vel_ang, dconfig.vel_turn, dconfig.number_of_angles, number_of_repeats, path_type, debug, 
					dconfig.max_path_angle);
				printf("pathfile = %s\n", dconfig.pathfile.c_str());
				printf("pointfile = %s\n", dconfig.pointfile.c_str());
				printf("outdir = %s\n", dconfig.outdir.c_str());
			}
			else
			{
				printf("options element not found.\n");
				fflush(stdout);
				error_count++;
			}
              
			if (error_count > 0 || dconfig.vel_ang <= 0)
			{
				printf("ReadConfig failed. Some configuraiton options are invalid.\n");
			}
			else
			{
				bres = true;
				printf("ReadConfig succeeded.\n");
			}
		}
		else
		{
			printf("ReadConfig failed. XML file is invalid: %s.\n", ConfigXmlFile.c_str());
		}
		doc.Clear();
	}
	else
	{
		printf("ReadConfig failed. XML load failed for %s.\n", ConfigXmlFile.c_str());
	}
	fflush(stdout);
	return bres;
}

string IntToStr(int i)
{	
	string s = "";
	char buf[sizeof(int)*8+1];
	sprintf(buf, "%d", i);	
	s = buf;
	return s;
}

int DeleteFiles(string dir)
{
	struct dirent *dirp; 
	DIR *dp;
	struct stat buf;
	int count = 0;
	int lendir = dir.length();

	if (dir[lendir-1] != '/') dir = dir + "/";

	if((dp  = opendir(dir.c_str())) == NULL) 
	{
	        cout << "Error(" << errno << ") opening " << dir << endl;
       		return errno;
	}
		
	while ((dirp = readdir(dp)) != NULL) 
	{
		//printf("d_name = %s\n", dirp->d_name);
		
		int len = strlen(dirp->d_name);
		if (dirp->d_name[len-1] == '~' || strcmp(dirp->d_name, ".") == 0 || strcmp(dirp->d_name, "..") == 0) continue;

	        string fname = dir + dirp->d_name;
		//printf("%d %s\n", count, fname.c_str());
		
		if (stat(fname.c_str(), &buf) < 0) continue;

	        bool is_dir = S_ISDIR(buf.st_mode);
		if (is_dir) 
		{ 
			printf("%s is a directory.\n", fname.c_str());
			fflush(stdout);
		}
		else
		{
			count++;
			printf("%s is a file.\n", fname.c_str());
			fflush(stdout);
			//remove(fname.c_str());
		}
	}
	printf("count = %d\n", count);
	return count;
}








		
