#define	USE_SSE
#include	<gsl/gsl_randist.h>
#include	<algorithm>
#include	<iostream>
#include	<stdint.h>
#include	<unistd.h>
#include	<cstdlib>
#include	<fstream>
#include	<sstream>
#include	"lbfgs.c"
#include	"sse.cpp"
#include	<cfloat>
#include	<string>
#include	<vector>
#include	<cmath>
#include	<ctime>
#include	<omp.h>
using	namespace	std;
const	float	missing=1e15;

typedef	float	v4sf	__attribute__	((__vector_size__	(16)));
union	XMM{
    uint32_t	i[4];
    v4sf	xmm;
}mask[16]={
	{  0U,  0U,  0U,  0U },	{ ~0U,  0U,  0U,  0U },	{  0U, ~0U,  0U,  0U },	{ ~0U, ~0U,  0U,  0U },
 	{  0U,  0U, ~0U,  0U },	{ ~0U,  0U, ~0U,  0U },	{  0U, ~0U, ~0U,  0U },	{ ~0U, ~0U, ~0U,  0U },
 	{  0U,  0U,  0U, ~0U },	{ ~0U,  0U,  0U, ~0U },	{  0U, ~0U,  0U, ~0U },	{ ~0U, ~0U,  0U, ~0U },
 	{  0U,  0U, ~0U, ~0U },	{ ~0U,  0U, ~0U, ~0U },	{  0U, ~0U, ~0U, ~0U },	{ ~0U, ~0U, ~0U, ~0U },
};

class	RBR_Model{
private:
	vector<uint32_t>	bits,	fids;
	float	*dbeta;
	size_t	sample,	feature,	asize,	fsize;
	vector<float>	mean,	prec,	impo,	fwei,	thre,	weig;
	bool	binary,	pret;

	void	set_one(vector<uint32_t>	&V,	uint64_t	I){	V[I>>5]|=1U<<(I&31);	}
	lbfgsfloatval_t	function_binary(const	lbfgsfloatval_t	*x,	lbfgsfloatval_t	*g);
	static	lbfgsfloatval_t	evaluate_binary(void	*instance,	const	lbfgsfloatval_t	*x,	lbfgsfloatval_t	*g,	const int n,	const lbfgsfloatval_t step){
		return	reinterpret_cast<RBR_Model*>(instance)->function_binary(x, g);
	}
	lbfgsfloatval_t	function_real(const	lbfgsfloatval_t	*x,	lbfgsfloatval_t	*g);
	static	lbfgsfloatval_t	evaluate_real(void	*instance,	const	lbfgsfloatval_t	*x,	lbfgsfloatval_t	*g,	const int n,	const lbfgsfloatval_t step){
		return	reinterpret_cast<RBR_Model*>(instance)->function_real(x, g);
	}
    static	int	progress(void	*instance,	const	lbfgsfloatval_t	*x,	const lbfgsfloatval_t *g,	const	lbfgsfloatval_t	fx,	const	lbfgsfloatval_t	xnorm,	const lbfgsfloatval_t gnorm,	const	lbfgsfloatval_t	step,	int n,	int	k,	int	ls){
		cerr<<k<<'\t'<<fx<<"        \r";
		return	0;
	}
public:
	vector<float>	yvec;
	float	*xmat,	ridge;
	size_t	bsize,	twist,	threads,	correction,	max_iteration;
	string	weight_file;

	void	document(void);
	bool	load_data(const	char	*X,	const	char	*Y);
	bool	load_weight(void);
	void	standardization(void);
	void	features(void);
	void	estimate(void);
	bool	save_model(const	char	*M);
};

bool	RBR_Model::load_data(const	char	*X,	const	char	*Y){
	cout<<"loading data...\n";
	vector<float>	temp;
	ifstream	fx(X);
	if(!fx){	cerr<<"fail to open "<<X<<endl;	return	false;	}
	sample=0;	string	buf;	float	x;
	while(!fx.eof()){
		getline(fx,	buf);	if(!buf.size())	continue;
		for(size_t	i=0;	i<buf.size();	i++)	if(buf[i]==','||buf[i]==';')	buf[i]=' ';
		istringstream	si(buf);
		while(!si.eof()){	x=-FLT_MAX;	si>>x;	if(x>-FLT_MAX)	temp.push_back(x);	}
		if(!(++sample%1000))	cerr<<sample/1000<<"k\r";
	}
	fx.close();
	
	ifstream	fy(Y);
	if(!fy){	cerr<<"fail to open "<<Y<<endl;	return	false;	}
	for(fy>>x;	!fy.eof();	fy>>x)	yvec.push_back(x);
	fy.close();
	
	if(temp.size()%sample){	cerr<<"rows with different number of columns in "<<X<<endl;	return	false;	}
	feature=temp.size()/sample;
	if(yvec.size()!=sample){	cerr<<X<<" and "<<Y<<" has different number of rows\n";	return	false;	}
	cout<<"sample\t"<<sample<<"\nfeature\t"<<feature<<endl;
	
	asize=a_size_f(sample);
	if(posix_memalign((void**)&xmat,	16,	feature*asize*4)){	cerr<<"not enough memory\n";	return	false;	}
	for(size_t	i=0;	i<sample;	i++)	for(size_t	j=0;	j<feature;	j++)	xmat[j*asize+i]=temp[i*feature+j];
	return	true;
}

bool	RBR_Model::load_weight(void){
	impo.assign(feature,	1);
	if(!weight_file.size())	return	false;
	ifstream	fi(weight_file.c_str());
	if(!fi){	cerr<<"fail to open "<<weight_file<<endl;		return	false;	}
	for(size_t	i=0;	i<feature;	i++){	fi>>impo[i];	impo[i]=fmaxf(0,	impo[i]);	}
	fi.close();
	return	true;
}

void	RBR_Model::standardization(void){
	cout<<"standardization...\n";
	mean.resize(feature);	prec.resize(feature);
	#pragma omp parallel for
	for(size_t	i=0;	i<feature;	i++){
		double	sx=0,	sxx=0;	size_t	sn=0;	float	*p=xmat+i*asize;
		for(size_t	j=0;	j<sample;	j++)	if(fabs(p[j])<missing){	sx+=p[j];	sxx+=p[j]*p[j];	sn++;	}
		if(!sn)	sx=sxx=0;
		else{	sx/=sn;	sxx=sxx-sx*sx*sn>0?sqrt(sn/(sxx-sx*sx*sample)):0;	}
		mean[i]=sx;	prec[i]=sxx;
		for(size_t	j=0;	j<sample;	j++)	p[j]=fabs(p[j])<missing?sxx*(p[j]-sx):0;
	}
}

void	RBR_Model::features(void){
	cout<<"feature generation...\n";
	fsize=(bsize+1)%32?((bsize+1)/32+1)*32:bsize+1;
	bits.resize(sample*(fsize>>5));
	fids.resize(bsize*twist);	fwei.resize(bsize*twist);	thre.resize(bsize);
	float	*temp;	if(posix_memalign((void**)&temp,	16,	threads*asize*4))	return;	
	size_t	t0=time(NULL);	vector<gsl_rng*>	rng(threads);
	for(size_t	i=0;	i<threads;	i++){	rng[i]=gsl_rng_alloc(gsl_rng_default);	gsl_rng_set(rng[i],	t0+i);	}

	float	sumw=0;	for(size_t	i=0;	i<feature;	i++)	sumw+=impo[i];
	for(size_t	i=0;	i<sample;	i++)	set_one(bits,	(uint64_t)i*fsize+bsize);
	size_t	block=fsize/32;
	#pragma omp parallel for
	for(size_t	b=0;	b<block;	b++){
		for(size_t	f=b*32;	f<b*32+32;	f++)	if(f<bsize){
			uint32_t	*fid=&fids[f*twist];
			float	*fwe=&fwei[f*twist],	*tem=temp+omp_get_thread_num()*asize;
			gsl_rng	*r=rng[omp_get_thread_num()];
			for(size_t	i=0;	i<twist;	i++){
				fwe[i]=gsl_ran_ugaussian(r);
				float	ran=gsl_rng_uniform(r)*sumw,	sum=0;
				for(size_t	j=0;	j<feature;	j++){	sum+=impo[j];	if(sum>=ran){	fid[i]=j;	break;	}	}
			}
			sort(fid,	fid+twist);
			memset(tem,	0,	sample*4);
			for(size_t	i=0;	i<twist;	i++)	vec_con_fma_f(tem,	xmat+fid[i]*asize,	fwe[i],	sample);
			float	cut=thre[f]=tem[gsl_rng_get(r)%sample];
			for(size_t	i=0;	i<sample;	i++)	if(tem[i]>=cut)	set_one(bits,	(uint64_t)i*fsize+f);
			if(!((f+1)%10000))	cerr<<'=';
		}
	}
	cerr<<endl;
	for(size_t	i=0;	i<threads;	i++)	gsl_rng_free(rng[i]);
	free(temp);	free(xmat);
}

lbfgsfloatval_t	RBR_Model::function_binary(const	lbfgsfloatval_t	*x,	lbfgsfloatval_t	*g){
	float	*weight=dbeta+threads*fsize,	like=0;
	memset(dbeta,	0,	threads*fsize*4);	memcpy(weight,	x,	fsize*4);
	#pragma omp parallel for
	for(size_t	i=0;	i<sample;	i++)	if(!pret||!(i&15)){
		size_t	k=omp_get_thread_num();
		uint32_t	*p=&bits[i*(fsize>>5)];
		float	*q=dbeta+k*fsize;
		v4sf	sum={0,0,0,0};
		for(size_t	j=0;	j<fsize;	j+=32){
			uint32_t	b=p[j>>5];
			sum=__builtin_ia32_addps(sum,	__builtin_ia32_andps(*(v4sf*)(weight+j),	mask[b&0xf].xmm));
			sum=__builtin_ia32_addps(sum,	__builtin_ia32_andps(*(v4sf*)(weight+j+4),	mask[(b>>4)&0xf].xmm));
			sum=__builtin_ia32_addps(sum,	__builtin_ia32_andps(*(v4sf*)(weight+j+8),	mask[(b>>8)&0xf].xmm));
			sum=__builtin_ia32_addps(sum,	__builtin_ia32_andps(*(v4sf*)(weight+j+12),	mask[(b>>12)&0xf].xmm));
			sum=__builtin_ia32_addps(sum,	__builtin_ia32_andps(*(v4sf*)(weight+j+16),	mask[(b>>16)&0xf].xmm));
			sum=__builtin_ia32_addps(sum,	__builtin_ia32_andps(*(v4sf*)(weight+j+20),	mask[(b>>20)&0xf].xmm));
			sum=__builtin_ia32_addps(sum,	__builtin_ia32_andps(*(v4sf*)(weight+j+24),	mask[(b>>24)&0xf].xmm));
			sum=__builtin_ia32_addps(sum,	__builtin_ia32_andps(*(v4sf*)(weight+j+28),	mask[(b>>28)&0xf].xmm));
		}
		float	*r=(float*)&sum,	z=r[0]+r[1]+r[2]+r[3];
		z=1.0f/(1.0f+expf(-z));
		float	dl=yvec[i]?logf(fmaxf(z,	FLT_MIN)):logf(fmaxf(1.0f-z,	FLT_MIN));
		#pragma omp atomic
		like+=dl;
		r[0]=r[1]=r[2]=r[3]=yvec[i]-z;
		for(size_t	j=0;	j<fsize;	j+=32){
			uint32_t	b=p[j>>5];	
			*(v4sf*)(q+j)=__builtin_ia32_addps(*(v4sf*)(q+j),	__builtin_ia32_andps(sum,	mask[b&0xf].xmm));
			*(v4sf*)(q+j+4)=__builtin_ia32_addps(*(v4sf*)(q+j+4),	__builtin_ia32_andps(sum,	mask[(b>>4)&0xf].xmm));
			*(v4sf*)(q+j+8)=__builtin_ia32_addps(*(v4sf*)(q+j+8),	__builtin_ia32_andps(sum,	mask[(b>>8)&0xf].xmm));
			*(v4sf*)(q+j+12)=__builtin_ia32_addps(*(v4sf*)(q+j+12),	__builtin_ia32_andps(sum,	mask[(b>>12)&0xf].xmm));
			*(v4sf*)(q+j+16)=__builtin_ia32_addps(*(v4sf*)(q+j+16),	__builtin_ia32_andps(sum,	mask[(b>>16)&0xf].xmm));
			*(v4sf*)(q+j+20)=__builtin_ia32_addps(*(v4sf*)(q+j+20),	__builtin_ia32_andps(sum,	mask[(b>>20)&0xf].xmm));
			*(v4sf*)(q+j+24)=__builtin_ia32_addps(*(v4sf*)(q+j+24),	__builtin_ia32_andps(sum,	mask[(b>>24)&0xf].xmm));
			*(v4sf*)(q+j+28)=__builtin_ia32_addps(*(v4sf*)(q+j+28),	__builtin_ia32_andps(sum,	mask[(b>>28)&0xf].xmm));
		}
	}
	for(size_t	j=0;	j<fsize;	j++){	g[j]=j?2*ridge*x[j]:0;	like-=j?ridge*x[j]*x[j]:0;	}
	for(size_t	i=0;	i<threads;	i++){
		float	*p=dbeta+i*fsize;
		for(size_t	j=0;	j<fsize;	j++)	g[j]-=p[j];
	}
	return	-like;
}

lbfgsfloatval_t	RBR_Model::function_real(const	lbfgsfloatval_t	*x,	lbfgsfloatval_t	*g){
	float	*weight=dbeta+threads*fsize,	square=0;
	memset(dbeta,	0,	threads*fsize*4);	memcpy(weight,	x,	fsize*4);
	#pragma omp parallel for
	for(size_t	i=0;	i<sample;	i++)	if(!pret||!(i&15)){
		size_t	k=omp_get_thread_num();
		uint32_t	*p=&bits[i*(fsize>>5)];
		float	*q=dbeta+k*fsize;
		v4sf	sum={0,0,0,0};
		for(size_t	j=0;	j<fsize;	j+=32){
			uint32_t	b=p[j>>5];
			sum=__builtin_ia32_addps(sum,	__builtin_ia32_andps(*(v4sf*)(weight+j),	mask[b&0xf].xmm));
			sum=__builtin_ia32_addps(sum,	__builtin_ia32_andps(*(v4sf*)(weight+j+4),	mask[(b>>4)&0xf].xmm));
			sum=__builtin_ia32_addps(sum,	__builtin_ia32_andps(*(v4sf*)(weight+j+8),	mask[(b>>8)&0xf].xmm));
			sum=__builtin_ia32_addps(sum,	__builtin_ia32_andps(*(v4sf*)(weight+j+12),	mask[(b>>12)&0xf].xmm));
			sum=__builtin_ia32_addps(sum,	__builtin_ia32_andps(*(v4sf*)(weight+j+16),	mask[(b>>16)&0xf].xmm));
			sum=__builtin_ia32_addps(sum,	__builtin_ia32_andps(*(v4sf*)(weight+j+20),	mask[(b>>20)&0xf].xmm));
			sum=__builtin_ia32_addps(sum,	__builtin_ia32_andps(*(v4sf*)(weight+j+24),	mask[(b>>24)&0xf].xmm));
			sum=__builtin_ia32_addps(sum,	__builtin_ia32_andps(*(v4sf*)(weight+j+28),	mask[(b>>28)&0xf].xmm));
		}
		float	*r=(float*)&sum,	z=r[0]+r[1]+r[2]+r[3];
		z-=yvec[i];
		float	ds=z*z;
		#pragma omp atomic
		square+=ds;
		r[0]=r[1]=r[2]=r[3]=2*z;
		for(size_t	j=0;	j<fsize;	j+=32){
			uint32_t	b=p[j>>5];	
			*(v4sf*)(q+j)=__builtin_ia32_addps(*(v4sf*)(q+j),	__builtin_ia32_andps(sum,	mask[b&0xf].xmm));
			*(v4sf*)(q+j+4)=__builtin_ia32_addps(*(v4sf*)(q+j+4),	__builtin_ia32_andps(sum,	mask[(b>>4)&0xf].xmm));
			*(v4sf*)(q+j+8)=__builtin_ia32_addps(*(v4sf*)(q+j+8),	__builtin_ia32_andps(sum,	mask[(b>>8)&0xf].xmm));
			*(v4sf*)(q+j+12)=__builtin_ia32_addps(*(v4sf*)(q+j+12),	__builtin_ia32_andps(sum,	mask[(b>>12)&0xf].xmm));
			*(v4sf*)(q+j+16)=__builtin_ia32_addps(*(v4sf*)(q+j+16),	__builtin_ia32_andps(sum,	mask[(b>>16)&0xf].xmm));
			*(v4sf*)(q+j+20)=__builtin_ia32_addps(*(v4sf*)(q+j+20),	__builtin_ia32_andps(sum,	mask[(b>>20)&0xf].xmm));
			*(v4sf*)(q+j+24)=__builtin_ia32_addps(*(v4sf*)(q+j+24),	__builtin_ia32_andps(sum,	mask[(b>>24)&0xf].xmm));
			*(v4sf*)(q+j+28)=__builtin_ia32_addps(*(v4sf*)(q+j+28),	__builtin_ia32_andps(sum,	mask[(b>>28)&0xf].xmm));
		}
	}
	for(size_t	j=0;	j<fsize;	j++){	g[j]=j?2*ridge*x[j]:0;	square+=j?ridge*x[j]*x[j]:0;	}
	for(size_t	i=0;	i<threads;	i++){
		float	*p=dbeta+i*fsize;
		for(size_t	j=0;	j<fsize;	j++)	g[j]+=p[j];
	}
	return	square;
}

void	RBR_Model::estimate(void){
	binary=true;
	double	sx=0,	sxx=0;
	for(size_t	i=0;	i<sample;	i++){
		if(yvec[i]!=0&&yvec[i]!=1)	binary=false;
		sx+=yvec[i];	sxx+=yvec[i]*yvec[i];
	}
	sx/=sample;	sxx=sxx/sample-sx*sx;
	if(!sxx)	sxx=1;	
	ridge=binary?ridge*bsize:ridge*bsize/sxx;
	cout<<"type\t"<<(binary?"binary":"real")<<endl;

	if(posix_memalign((void**)&dbeta,	16,	(threads+1)*fsize*4)){	cerr<<"not enough memory\n";	return;	}	
	lbfgsfloatval_t	score,	*x=lbfgs_malloc(fsize);
	memset(x,	0,	fsize*4);
	lbfgs_parameter_t	para;	lbfgs_parameter_init(&para);	
	para.m=correction;	para.max_iterations=max_iteration;	para.epsilon=0;
	
	pret=true;
	if(binary)	lbfgs(fsize,	x,	&score,	evaluate_binary,	progress,	this,	&para);
	else	lbfgs(fsize,	x,	&score,	evaluate_real,	progress,	this,	&para);
	cout<<endl;
	pret=false;	
	if(binary)	lbfgs(fsize,	x,	&score,	evaluate_binary,	progress,	this,	&para);
	else	lbfgs(fsize,	x,	&score,	evaluate_real,	progress,	this,	&para);
	cout<<endl;	

	weig.resize(fsize);	memcpy(&weig[0],	x,	fsize*4);	lbfgs_free(x);
	free(dbeta);
	
	double	norm=0;
	for(size_t	i=0;	i<bsize;	i++)	norm+=weig[i]*weig[i];
	norm=1.0f/norm/twist;
	
	memset(&impo[0],	0,	feature*4);
	for(size_t	i=0;	i<bsize;	i++){
		float	w=weig[i]*weig[i]*norm;
		for(size_t	j=0;	j<twist;	j++)	impo[fids[i*twist+j]]+=w;
	}
}

struct	Sort{
	float	w;
	uint32_t	i;
	bool	operator()(Sort	X,	Sort	Y){	return	X.w>Y.w;	}
};

bool	RBR_Model::save_model(const	char	*M){
	ofstream	fo(M);
	if(!fo){	cerr<<"fail to open "<<M<<endl;	return	false;	}
	for(size_t	i=0;	i<feature;	i++)	fo<<impo[i]<<'\t';	fo<<'\n';
	fo<<"\nmean\tprec\n";
	for(size_t	i=0;	i<feature;	i++)	fo<<mean[i]<<'\t'<<prec[i]<<'\n';
	
	fo<<"\nbits\t"<<bsize<<'\n';
	fo<<"twist\t"<<twist<<'\n';
	fo<<"type\t"<<(binary?"binary\n":"real\n");
	fo<<"\nweight\tthreshold";	for(size_t	i=0;	i<twist;	i++)	fo<<"\tfeature"<<i<<"\tweight"<<i;	fo<<'\n';
	vector<Sort>	s(bsize);
	for(size_t	i=0;	i<bsize;	i++){	s[i].w=fabsf(weig[i]);	s[i].i=i;	}
	sort(s.begin(),	s.end(),	Sort());
	for(size_t	k=0;	k<bsize;	k++){
		size_t	i=s[k].i;
		fo<<weig[i]<<'\t'<<thre[i];
		for(size_t	j=0;	j<twist;	j++)	fo<<'\t'<<fids[i*twist+j]<<'\t'<<fwei[i*twist+j];
		fo<<'\n';
	}
	fo<<weig[bsize]<<'\n';
	fo.close();
	return	true;
}

void	RBR_Model::document(void){
	cout<<"Usage:	rbr_model [options] trainX trainY model\n";
	cout<<"\t-b:	number of random bits.		default=10000\n";	
	cout<<"\t-r:	Scaled L2 regularization.	default=0.1\n";
	cout<<"\t-t:	number of features twisted.	default=2\n";
	cout<<"\t-c:	L-BFGS corrections.		default=256\n";
	cout<<"\t-m:	maximum L-BFGS iteration.	default=0 (0=no limit)\n";
	cout<<"\t-T:	number of threads.		default=0 (all)\n";
	cout<<"\t-w:	variable weights file.\n";
	cout<<"\ttrainX:	input files in csv/tsv format without header.\n";
	cout<<"\ttrainY:	input files in plain text.\n";
	cout<<"\tmodel:	model output file.\n";
	exit(0);
}

int	main(int	ac,	char	**av){
	cout<<"********************************\n";
	cout<<"Random Bits Regression\n";
	cout<<"author:	Yi Wang\n";
	cout<<"email:	godspeed.wang@gmail.com\n";
	cout<<"date:	13/Feb/2014\n";
	cout<<"********************************\n";
	
	RBR_Model	m;
	m.bsize=10000;	m.ridge=0.01;	m.twist=2;	
	m.threads=omp_get_num_procs();	m.correction=256;	m.max_iteration=0;
	
	int	opt;
	while((opt=getopt(ac,	av,	"b:r:t:c:m:T:w:"))>=0){
		switch(opt){
		case	'b':	m.bsize=atoi(optarg);	break;
		case	'r':	m.ridge=atof(optarg);	break;
		case	't':	m.twist=atoi(optarg);	break;
		case	'c':	m.correction=atoi(optarg);	break;
		case	'm':	m.max_iteration=atoi(optarg)>0?atoi(optarg):0;	break;
		case	'T':	m.threads=atoi(optarg)>0?atoi(optarg):omp_get_num_procs();	break;
		case	'w':	m.weight_file=optarg;	break;
		default:	m.document();
		}
	}
	if(ac<optind+3)	m.document();
	omp_set_num_threads(m.threads);
	
	if(!m.load_data(av[optind],	av[optind+1]))	return	0;
	m.standardization();	
	if(!m.load_weight())	cerr<<"using uniform feature weights\n";
	else	cerr<<"using custom feature weights\n";
	m.features();
	m.estimate();
	if(!m.save_model(av[optind+2]))	return	0;	
	return	0;
}

