#include <iostream>
#include <vector>
#include <cstdlib>
#include "stdio.h"
#include <fstream>
#include "assert.h"
#include <algorithm>
#include <cstring>
#include <cmath>

bool comparator( std::pair< double, double > l, std::pair< double, double > r){ return l.first < r.first; }
bool comparator2( std::pair< std::pair<double,double>, int > l, std::pair< std::pair<double,double>, int > r){ return l.first.first < r.first.first; }

#define NUM 1
#define t_step_max 100
#define t_step_interval 1
#define xcen 1
#define ycen 1
#define zcen 1
#define print_interval 0.0001 
#define num_1D 3000
#define radius 0.5
#define PMIN 5
#define PMAX 15

using namespace std;

const char* right_flush(int	n, int	ndigits) {
	static	char	s[20];
	int		i;

	if (n == 0)
		ndigits--;
	for (i = n; i > 0; i /= 10) ndigits--;

	s[0] = '\0';
	for (;ndigits > 0; ndigits--)
		(void) sprintf(s,"%s0",s);
	(void) sprintf(s,"%s%d",s,n);
	return s;
}

int main() {
	vector< pair<double,double> > p_norm, rho_norm;
	for(int t_step=0; t_step<=t_step_max; t_step+=t_step_interval) {
		cout<<"time="<<t_step*print_interval<<endl;
		
		ifstream f[NUM];
		char filename[NUM][100];
		const char* outputname = "out";
		int num_of_par[NUM];
		vector<double> xp, yp, zp, up, vp, wp, p, rho, e;

		for(int n=0; n<NUM; n++) {
			sprintf(filename[n],"out%s.vtk",right_flush(t_step,7));
			cout<<"filename = "<<filename[n]<<endl;
			f[n].open(filename[n]);
			
			double x, y, z, u, v, w, pressure, density, energy;
			char c1[50], c2[50];
			string s;
		
			// Skip the first 4 lines
			for(int line=1; line<=4; line++) {
				getline(f[n],s);
				cout<<s<<endl;
			}
		
			// Read the number of particles
			f[n] >> c1 >> num_of_par[n] >> c2;
//			cout<<c1<<" "<<num_of_par[n]<<" "<<c2<<endl;
		
			// Read location
			for(int i=0; i<num_of_par[n]; i++) {
				f[n] >> x >> y >> z;
				xp.push_back(x);
				yp.push_back(y);
				zp.push_back(z);
	//			cout<<x<<" "<<y<<" "<<z<<endl;
			}
			cout<<"xp.size()="<<xp.size()<<endl;
		
			// Skip 2 additional lines (Note the first will be a blank line so we do 3 times)
			for(int line=1; line<=3; line++) {
				getline(f[n],s);
//				cout<<s<<endl;
			}
		
			// Read velocity
			for(int i=0; i<num_of_par[n]; i++) {
				f[n] >> u >> v >> w;
				up.push_back(u);
				vp.push_back(v);
				wp.push_back(w);
	//			cout<<u<<" "<<v<<" "<<w<<endl;
			}
//			cout<<"up.size()="<<up.size()<<endl;
		
			// Skip 2 additional lines (Note the first will be a blank line so we do 3 times)
			for(int line=1; line<=3; line++) {
				getline(f[n],s);
//				cout<<s<<endl;
			}
		
			// Read pressure
			for(int i=0; i<num_of_par[n]; i++) {
				f[n] >> pressure;
				p.push_back(pressure);
	//			cout<<pressure<<endl;
			}
//			cout<<"p.size()="<<p.size()<<endl;
		
			// Skip 2 additional lines (Note the first will be a blank line so we do 3 times)
			for(int line=1; line<=3; line++) {
				getline(f[n],s);
//				cout<<s<<endl;
			}
		
			// Read density
			for(int i=0; i<num_of_par[n]; i++) {
				f[n] >> density;
				rho.push_back(1./density);
	//			cout<<density<<endl;
			}
//			cout<<"rho.size()="<<rho.size()<<endl;
			

		}
	
		vector<double> dis;
		vector<double> vel;
		for(int i=0; i<xp.size(); i++) {
			dis.push_back(sqrt((xp[i]-xcen)*(xp[i]-xcen) + (yp[i]-ycen)*(yp[i]-ycen) + (zp[i]-zcen)*(zp[i]-zcen)));
			vel.push_back(sqrt(up[i]*up[i] + vp[i]*vp[i] + wp[i]*wp[i]));
		}
//		cout<<"dis.size()="<<dis.size()<<endl;
//		cout<<"vel.size()="<<vel.size()<<endl;
	
		vector< pair<double,double> > dis_vel;
		vector< pair<double,double> > dis_p;
		vector< pair<double,double> > dis_rho;
		vector< pair<double,double> > dis_e;
	
		for(int i=0; i<xp.size(); i++) {
			double signed_dis = dis[i];
			if(xp[i]<xcen) signed_dis*=-1.0;
//			if(dis[i]<=radius) {
				dis_vel.push_back(make_pair(signed_dis,vel[i]));
				dis_p.push_back(make_pair(signed_dis,p[i]));
				dis_rho.push_back(make_pair(signed_dis,rho[i]));
				dis_e.push_back(make_pair(signed_dis,e[i]));
//			}
		}
	
//		cout<<"dis_vel.size()="<<dis_vel.size()<<endl;
//		cout<<"dis_p.size()="<<dis_p.size()<<endl;
//		cout<<"dis_rho.size()="<<dis_rho.size()<<endl;
//		cout<<"dis_e.size()="<<dis_e.size()<<endl;
	
		sort(dis_vel.begin(),dis_vel.end(),comparator);
		sort(dis_p.begin(),dis_p.end(),comparator);
		sort(dis_rho.begin(),dis_rho.end(),comparator);
		sort(dis_e.begin(),dis_e.end(),comparator);
		
		// Open files to write to
		FILE *f_vel, *f_p, *f_rho, *f_e;
		char fn1[50], fn2[50], fn3[50], fn4[50];
		sprintf(fn1,"z002_vel_%.4g.txt",t_step*print_interval);
		sprintf(fn2,"z002_p_bd%.4g.txt",t_step*print_interval);
		sprintf(fn3,"z002_rho_bd%.4g.txt",t_step*print_interval);
		sprintf(fn4,"z002_e_%.4g.txt",t_step*print_interval);
//		f_vel = fopen(fn1,"w");
		f_p = fopen(fn2,"w");
		f_rho = fopen(fn3,"w");
//		f_e = fopen(fn4,"w");
	
		for(int i=0; i<dis_p.size(); i++) {
	//		cout<<"dis="<<dis_p[i].first<<" p="<<dis_p[i].second<<endl;
//			fprintf(f_vel,"%.16g %.16g\n",dis_vel[i].first,dis_vel[i].second);
			fprintf(f_p,"%.16g %.16g\n",dis_p[i].first,dis_p[i].second);
			fprintf(f_rho,"%.16g %.16g\n",dis_rho[i].first,dis_rho[i].second);
//			fprintf(f_e,"%.16g %.16g\n",dis_e[i].first,dis_e[i].second);
		}
	
//		fclose(f_vel);
		fclose(f_p);
		fclose(f_rho);
//		fclose(f_e);
		
//		vector< pair< pair<double,double>,int > > dis_p_tag;
//		vector< pair< pair<double,double>,int > > dis_rho_tag;
//		for(int i=0; i<dis_p_1D.size(); i++) {
//			dis_p_tag.push_back(make_pair(make_pair(dis_p_1D[i].first,dis_p_1D[i].second),1));
//			dis_rho_tag.push_back(make_pair(make_pair(dis_rho_1D[i].first,dis_rho_1D[i].second),1));
//		}
//		for(int i=0; i<dis_p.size(); i++) {
//			dis_p_tag.push_back(make_pair(make_pair(dis_p[i].first,dis_p[i].second),0));
//			dis_rho_tag.push_back(make_pair(make_pair(dis_rho[i].first,dis_rho[i].second),0));
//		}
//		cout<<"dis_p_tag.size()="<<dis_p_tag.size()<<endl;
//		cout<<"dis_rho_tag.size()="<<dis_rho_tag.size()<<endl;
//		
//		sort(dis_p_tag.begin(),dis_p_tag.end(),comparator2);
//		sort(dis_rho_tag.begin(),dis_rho_tag.end(),comparator2);
//		
//		vector< pair<double,double> > p_error;
//		vector< pair<double,double> > rho_error;
//		for(int i=0; i<dis_p_tag.size(); i++) {
//			int tag = dis_p_tag[i].second;
////			cout<<"tag="<<tag<<endl;
//			if(tag == 1) {
//				bool no_left = false;
//				bool no_right = false;
//				pair<double,double> left_p, right_p, left_rho, right_rho;
//				int j = i-1;
//				while(j>=0) {
//					if(dis_p_tag[j].second == 0) {
//						left_p = make_pair(dis_p_tag[j].first.first,dis_p_tag[j].first.second);
//						left_rho = make_pair(dis_rho_tag[j].first.first,dis_rho_tag[j].first.second);
//						break;
//					}
//					else
//						j--;
//				}
//				if(j<0) {
//					left_p = make_pair(0,PMAX);
////					left_rho = make_pair(0,PMAX);
//					no_left = true;
//				}
//				
//				j = i+1;
//				while(j<dis_p_tag.size()) {
//					if(dis_p_tag[j].second == 0) {
//						right_p = make_pair(dis_p_tag[j].first.first,dis_p_tag[j].first.second);
//						right_rho = make_pair(dis_rho_tag[j].first.first,dis_rho_tag[j].first.second);
//						break;
//					}
//					else
//						j++;
//				}
//				if(j>=dis_p_tag.size()) {
//					right_p = make_pair(radius,PMIN+(PMAX-PMIN)*exp(-100.*radius*radius));
////					right_rho = make_pair(radius,PMIN+(PMAX-PMIN)*exp(-100.*radius*radius));
//					no_right = true;
//				}	
//				double total_dis = right_p.first - left_p.first;
//				double left_dis = dis_p_tag[i].first.first - left_p.first;
//				double right_dis = right_p.first - dis_p_tag[i].first.first;
//				double p_inter = left_p.second*(right_dis/total_dis) + right_p.second*(left_dis/total_dis);
//				double p_err = dis_p_tag[i].first.second - p_inter;
//				double rho_inter = left_rho.second*(right_dis/total_dis) + right_rho.second*(left_dis/total_dis);
//				double rho_err = dis_rho_tag[i].first.second - rho_inter;
////				cout<<"total_dis="<<total_dis<<endl;
////				cout<<"left_dis="<<left_dis<<endl;
////				cout<<"right_dis="<<right_dis<<endl;
////				cout<<"error="<<p_err<<endl;
//				if(no_left == false && no_right == false) {
////					if(dis_p_tag[i].first.first<=0.25) {
//						p_error.push_back(make_pair(dis_p_tag[i].first.first,p_err));
//						rho_error.push_back(make_pair(dis_rho_tag[i].first.first,rho_err));
////					}
//				}
//			}
//		}
//		cout<<"p_error.size()="<<p_error.size()<<endl;
//		cout<<"rho_error.size()="<<rho_error.size()<<endl;
//		
//		// Open files to write to
//		FILE *f_p_error, *f_rho_error;
//		char fn_p_error[50], fn_rho_error[50];
//		sprintf(fn_p_error,"z002_p_error%.4g.txt",t_step*print_interval);
//		sprintf(fn_rho_error,"z002_rho_error%.4g.txt",t_step*print_interval);
//		f_p_error = fopen(fn_p_error,"w");
//		f_rho_error = fopen(fn_rho_error,"w");
//	
//		for(int i=0; i<p_error.size(); i++) {
//			fprintf(f_p_error,"%.16g %.16g\n",p_error[i].first,p_error[i].second);
//			fprintf(f_rho_error,"%.16g %.16g\n",rho_error[i].first,rho_error[i].second);
//		}
//	
//		fclose(f_p_error);
//		fclose(f_rho_error);
//		
//		double pe = 0;
//		double rhoe = 0;
//		for(int i=0; i<p_error.size(); i++) {
////			if(p_error[i].first>0.05 && p_error[i].first<0.25)
//				pe += p_error[i].second*p_error[i].second;
//				rhoe += rho_error[i].second*rho_error[i].second;
//		}
//		pe /= p_error.size();
//		rhoe /= rho_error.size();
//		p_norm.push_back(make_pair(t_step*print_interval,sqrt(pe)));
//		rho_norm.push_back(make_pair(t_step*print_interval,sqrt(rhoe)));
	}
//	
//	// Open files to write to
//	FILE *f_p_norm, *f_rho_norm;
//	char fn_p_norm[50], fn_rho_norm[50];
//	sprintf(fn_p_norm,"z002_p_norm.txt");
//	sprintf(fn_rho_norm,"z002_rho_norm.txt");
//	f_p_norm = fopen(fn_p_norm,"w");
//	f_rho_norm = fopen(fn_rho_norm,"w");
//	
//	for(int i=0; i<p_norm.size(); i++) {
//		fprintf(f_p_norm,"%.16g %.16g\n",p_norm[i].first,p_norm[i].second);
//		fprintf(f_rho_norm,"%.16g %.16g\n",rho_norm[i].first,rho_norm[i].second);
//	}
//	
//	fclose(f_p_norm);
//	fclose(f_rho_norm);
	
	return 0;

}
