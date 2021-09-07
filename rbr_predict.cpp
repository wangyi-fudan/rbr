#include	<iostream>
#include	<stdint.h>
#include	<fstream>
#include	<sstream>
#include	<cstdlib>
#include	<cfloat>
#include	<string>
#include	<vector>
#include	<cmath>
#include	<omp.h>
using	namespace	std;
const	float	missing=1e15;
const	float	epsilon=1e-4;

int	main(int	ac,	char	**av){
	if(ac<2){	cerr<<"rbr_predict model [threads] <testX >yhat\n";	return	0;	}
	
	size_t	feature,	bsize,	twist;
	vector<uint32_t>	fids;
	vector<float>	x,	mean,	prec,	fwei,	thre,	weig;
	bool	binary;
	float	bias;
	
	ifstream	fi(av[1]);
	if(!fi){	cerr<<"fail to open "<<av[1]<<endl;	return	0;	}
	string	temp;	getline(fi,	temp);
	istringstream	si(temp);
	float	t;	feature=0;
	while(!si.eof()){	t=-FLT_MAX;	si>>t;	if(t>-FLT_MAX)	feature++;	}
	mean.resize(feature);	prec.resize(feature);	x.resize(feature);
	fi>>temp>>temp;
	for(size_t	i=0;	i<feature;	i++)	fi>>mean[i]>>prec[i];
	fi>>temp>>bsize>>temp>>twist>>temp>>temp;
	binary=(temp=="binary");
	fwei.resize(bsize*twist);	fids.resize(bsize*twist);	thre.resize(bsize);	weig.resize(bsize);
	fi>>temp>>temp;
	for(size_t	i=0;	i<twist;	i++)	fi>>temp>>temp;
	for(size_t	i=0;	i<bsize;	i++){
		fi>>weig[i]>>thre[i];
		for(size_t	j=0;	j<twist;	j++)	fi>>fids[i*twist+j]>>fwei[i*twist+j];
	}
	fi>>bias;
	fi.close();
	
	size_t	npro=ac<3?omp_get_num_procs():atoi(av[2]);
	omp_set_num_threads(npro);
	vector<float>	sum(npro);
	while(!cin.eof()){
		getline(cin,	temp);	if(!temp.size())	continue;
		istringstream	si(temp);
		for(size_t	i=0;	i<feature;	i++){
			si>>x[i];
			x[i]=fabs(x[i])<missing?prec[i]*(x[i]-mean[i]):0;
		}
		for(size_t	i=0;	i<npro;	i++)	sum[i]=0;
		
		#pragma omp parallel for
		for(size_t	i=0;	i<bsize;	i++){
			float	s=0;
			float	*fwe=&fwei[i*twist];
			uint32_t	*fid=&fids[i*twist];
			for(size_t	j=0;	j<twist;	j++)	s+=x[fid[j]]*fwe[j];
			if(s+epsilon>=thre[i])	sum[omp_get_thread_num()]+=weig[i];			
		}
		float	y=bias;
		for(size_t	i=0;	i<npro;	i++)	y+=sum[i];
		if(binary)	y=1.0f/(1.0f+expf(-y));
		cout<<y<<'\n';
	}
	return	0;
}

