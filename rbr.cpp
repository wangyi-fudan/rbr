#define	USE_SSE
#include	<algorithm>
#include	<iostream>
#include	<stdint.h>
#include	<unistd.h>
#include	<cstdlib>
#include	<fstream>
#include	"lbfgs.c"
#include	<sstream>
#include	"xsa.cpp"
#include	<cfloat>
#include	<string>
#include	<vector>
#include	<cmath>
#include	<ctime>
#include	<omp.h>
using	namespace	std;

typedef	float	v4sf	__attribute__	((__vector_size__	(16)));
union	XMM{
    uint32_t	i[4];
    v4sf	xmm;
}mask[16]={
	{ 0U, 0U, 0U, 0U},	{~0U, 0U, 0U, 0U},	{ 0U,~0U, 0U, 0U},	{~0U,~0U, 0U, 0U},
 	{ 0U, 0U,~0U, 0U},	{~0U, 0U,~0U, 0U},	{ 0U,~0U,~0U, 0U},	{~0U,~0U,~0U, 0U},
 	{ 0U, 0U, 0U,~0U},	{~0U, 0U, 0U,~0U},	{ 0U,~0U, 0U,~0U},	{~0U,~0U, 0U,~0U},
 	{ 0U, 0U,~0U,~0U},	{~0U, 0U,~0U,~0U},	{ 0U,~0U,~0U,~0U},	{~0U,~0U,~0U,~0U},
};

class	RBR{
private:
	vector<float>	weight;
	vector<uint32_t>	fid,	bit,	sum1;
	bool	binary;
	float	*dbeta;
	
	void	set_one(vector<uint32_t>	&V,	uint64_t	I){	V[I>>5]|=1U<<(I&31);	}
	lbfgsfloatval_t	function(const	lbfgsfloatval_t	*x,	lbfgsfloatval_t	*g);
	static	lbfgsfloatval_t	evaluate(void	*instance,	const	lbfgsfloatval_t	*x,	lbfgsfloatval_t	*g,	const int n,	const lbfgsfloatval_t step){
		return	reinterpret_cast<RBR*>(instance)->function(x, g);
	}
    static	int	progress(void	*instance,	const	lbfgsfloatval_t	*x,	const lbfgsfloatval_t *g,	const	lbfgsfloatval_t	fx,	const	lbfgsfloatval_t	xnorm,	const lbfgsfloatval_t gnorm,	const	lbfgsfloatval_t	step,	int n,	int	k,	int	ls){
		cerr<<k<<'\t'<<(2-reinterpret_cast<RBR*>(instance)->binary)*fx/reinterpret_cast<RBR*>(instance)->trainn<<"       \r";
		return	0;
	}
	void	predict(const	float	*X);
public:
	vector<float>	trainy,	testy;
	size_t	trainn,	testn,	feature,	bits,	twist,	correction,	thread;	
	float	penalty;
	
	bool	load_weight(const	char	*F);
	bool	load_matrix(const	char	*F,	vector<float>	&M,	size_t	&R,	size_t	&C);
	void	x2bit(vector<float>	&TrX,	vector<float>	&TeX);
	void	estimate(const	char	*F);
	void	document(void);
};

bool	RBR::load_weight(const	char	*F){
	ifstream	fi(F);
	if(!fi){	cerr<<"fail to open "<<F<<'\n';	return	false;	}	//	W
	float	w;	weight.clear();
	for(fi>>w;	!fi.eof();	fi>>w)	weight.push_back(w);
	fi.close();
	return	true;
}

bool	RBR::load_matrix(const	char	*F,	vector<float>	&M,	uint64_t	&R,	uint64_t	&C) {
	ifstream	fi(F);
	if(!fi) {	cerr<<"fail to open "<<F<<'\n';	return	false;	}
	string	buf;	R=C=0;
	while(getline(fi,buf))	if(buf.size()) {
		char	*p=(char*)buf.data(),	*q;
		for(;;) {	q=p;	float	x=strtod(p,	&p);	if(p!=q)	M.push_back(x);	else	break;	}
		R++;
	}
	fi.close();
	if(M.size()%R) {	cerr<<"unequal column\t"<<F<<'\n';	return	false;	}
	C=M.size()/R;
	cout<<F<<'\t'<<R<<'*'<<C<<'\n';	
	return	true;
}

void	RBR::x2bit(vector<float>	&TrX,	vector<float>	&TeX){
	//	data preprocessing
	size_t	total=trainn+testn,	asize=total%16?(total/16+1)*16:total;
	float	*data;
	if(posix_memalign((void**)&data,	16,	feature*asize*4)){	cerr<<"not enough memory\n";	return;	}	//	G
	#pragma omp parallel for
	for(size_t	i=0;	i<feature;	i++){
		float	*p=data+i*asize;	size_t	n=total;	double	sx=0,	sxx=0;
		for(size_t	j=0;	j<trainn;j++){	p[j]=TrX[j*feature+i];	sx+=(double)p[j];	sxx+=(double)p[j]*(double)p[j];	}
		p=data+i*asize+trainn;
		for(size_t	j=0;	j<testn;	j++){	p[j]=TeX[j*feature+i];	sx+=(double)p[j];	sxx+=(double)p[j]*(double)p[j];	}
		p=data+i*asize;
		sx/=n;	sxx=sxx-sx*sx*n;	sxx=sxx>0?sqrt((n-1)/sxx):0;
		for(size_t	j=0;	j<total;	j++)	p[j]=sxx*(p[j]-sx);
	}
	vector<float>().swap(TrX);	vector<float>().swap(TeX);
	//	random bits generation
	if(twist>feature)	twist=feature;
	bits=bits%32?(bits/32+1)*32:bits;	size_t	block=bits/32;
	fid.resize(bits*twist);	sum1.resize(bits);	bit.resize(total*bits/32);
	if(weight.size()!=feature){	weight.assign(feature,	1);	cout<<"using uniform variable weights\n";	}
	else	cout<<"using custom variable weights\n";
	cout<<"bits\t"<<bits<<'\n';
	float	sumw=0;	for(size_t	i=0;	i<feature;	i++){	sumw+=weight[i];	weight[i]=sumw;	}
	vector<XSA>	rng(thread);	size_t	t0=time(NULL);	for(size_t	i=0;	i<thread;	i++)	rng[i].set(t0+i);
	float	*temp;	if(posix_memalign((void**)&temp,	16,	thread*asize*4)){	cerr<<"not enough memory\n";	return;	}	//	Y
	for(size_t	i=0;	i<total;	i++)	set_one(bit,	i*bits);
	
	#pragma omp parallel for
	for(size_t	b=0;	b<block;	b++)	for(size_t	f=b*32;	f<b*32+32;	f++)	if(f){
		float	*p=temp+omp_get_thread_num()*asize;
		XSA	&r=rng[omp_get_thread_num()];
		uint32_t	*id=&fid[f*twist];	
		for(size_t	i=0;	i<twist;	i++)	do	id[i]=lower_bound(weight.begin(),	weight.end(),	r.uniform_single()*sumw)-weight.begin();	while(find(id,	id+i,	id[i])!=id+i);
		sort(id,	id+twist);
		memset(p,	0,	total*4);
		for(size_t	i=0;	i<twist;	i++){
			float	*q=data+id[i]*asize,	w=r.fast_normal();	v4sf	z={w,w,w,w};
			for(size_t	j=0;	j<asize;	j+=16){
				*(v4sf*)(p+j+0)=__builtin_ia32_addps(*(v4sf*)(p+j+0),	__builtin_ia32_mulps(z,	*(v4sf*)(q+j+0)));
				*(v4sf*)(p+j+4)=__builtin_ia32_addps(*(v4sf*)(p+j+4),	__builtin_ia32_mulps(z,	*(v4sf*)(q+j+4)));
				*(v4sf*)(p+j+8)=__builtin_ia32_addps(*(v4sf*)(p+j+8),	__builtin_ia32_mulps(z,	*(v4sf*)(q+j+8)));
				*(v4sf*)(p+j+12)=__builtin_ia32_addps(*(v4sf*)(p+j+12),	__builtin_ia32_mulps(z,	*(v4sf*)(q+j+12)));
			}	
		}
		float	cut=p[r.get()%total];
		size_t	ones=0;
		for(size_t	i=0;	i<total;	i++)	if(p[i]>=cut){	set_one(bit,	i*bits+f);	ones++;	}
		sum1[f]=ones;
		if(!(f%10000))	cerr<<'=';
	}
	cerr<<'\n';
	free(temp);	free(data);
}

lbfgsfloatval_t	RBR::function(const	lbfgsfloatval_t	*x,	lbfgsfloatval_t	*g){
	double	loss=0;	memset(dbeta,	0,	thread*bits*4);
	#pragma omp parallel for
	for(size_t	i=0;	i<trainn;	i++){
		size_t	k=omp_get_thread_num();
		uint32_t	*p=&bit[i*(bits>>5)];
		float	*q=dbeta+k*bits;
		v4sf	sum={0,0,0,0};
		for(size_t	j=0;	j<bits;	j+=32){
			uint32_t	b=p[j>>5];
			sum=__builtin_ia32_addps(sum,	__builtin_ia32_andps(*(v4sf*)(x+j),	mask[b&0xf].xmm));
			sum=__builtin_ia32_addps(sum,	__builtin_ia32_andps(*(v4sf*)(x+j+4),	mask[(b>>4)&0xf].xmm));
			sum=__builtin_ia32_addps(sum,	__builtin_ia32_andps(*(v4sf*)(x+j+8),	mask[(b>>8)&0xf].xmm));
			sum=__builtin_ia32_addps(sum,	__builtin_ia32_andps(*(v4sf*)(x+j+12),	mask[(b>>12)&0xf].xmm));
			sum=__builtin_ia32_addps(sum,	__builtin_ia32_andps(*(v4sf*)(x+j+16),	mask[(b>>16)&0xf].xmm));
			sum=__builtin_ia32_addps(sum,	__builtin_ia32_andps(*(v4sf*)(x+j+20),	mask[(b>>20)&0xf].xmm));
			sum=__builtin_ia32_addps(sum,	__builtin_ia32_andps(*(v4sf*)(x+j+24),	mask[(b>>24)&0xf].xmm));
			sum=__builtin_ia32_addps(sum,	__builtin_ia32_andps(*(v4sf*)(x+j+28),	mask[(b>>28)&0xf].xmm));
		}
		float	*r=(float*)&sum,	z=r[0]+r[1]+r[2]+r[3],	dl;
		if(binary){
			z=1.0f/(1.0f+expf(-z));
			dl=trainy[i]?-logf(fmaxf(z,	FLT_MIN)):-logf(fmaxf(1.0f-z,	FLT_MIN));
		}
		else	dl=0.5*(z-trainy[i])*(z-trainy[i]);
		z-=trainy[i];
		#pragma omp atomic
		loss+=dl;
		r[0]=r[1]=r[2]=r[3]=z;
		for(size_t	j=0;	j<bits;	j+=32){
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
	for(size_t	j=0;	j<bits;	j++){	g[j]=j?penalty*x[j]:0;	loss+=j?0.5*penalty*x[j]*x[j]:0;	}
	for(size_t	i=0;	i<thread;	i++){
		float	*p=dbeta+i*bits;
		for(size_t	j=0;	j<bits;	j+=4)	*(v4sf*)(g+j)=__builtin_ia32_addps(*(v4sf*)(g+j),	*(v4sf*)(p+j));
	}
	return	loss;
}

void	RBR::predict(const	float	*X){
	testy.resize(testn);
	#pragma omp parallel for
	for(size_t	i=0;	i<testn;	i++){
		uint32_t	*p=&bit[(i+trainn)*(bits>>5)];
		v4sf	sum={0,0,0,0};
		for(size_t	j=0;	j<bits;	j+=32){
			uint32_t	b=p[j>>5];
			sum=__builtin_ia32_addps(sum,	__builtin_ia32_andps(*(v4sf*)(X+j),	mask[b&0xf].xmm));
			sum=__builtin_ia32_addps(sum,	__builtin_ia32_andps(*(v4sf*)(X+j+4),	mask[(b>>4)&0xf].xmm));
			sum=__builtin_ia32_addps(sum,	__builtin_ia32_andps(*(v4sf*)(X+j+8),	mask[(b>>8)&0xf].xmm));
			sum=__builtin_ia32_addps(sum,	__builtin_ia32_andps(*(v4sf*)(X+j+12),	mask[(b>>12)&0xf].xmm));
			sum=__builtin_ia32_addps(sum,	__builtin_ia32_andps(*(v4sf*)(X+j+16),	mask[(b>>16)&0xf].xmm));
			sum=__builtin_ia32_addps(sum,	__builtin_ia32_andps(*(v4sf*)(X+j+20),	mask[(b>>20)&0xf].xmm));
			sum=__builtin_ia32_addps(sum,	__builtin_ia32_andps(*(v4sf*)(X+j+24),	mask[(b>>24)&0xf].xmm));
			sum=__builtin_ia32_addps(sum,	__builtin_ia32_andps(*(v4sf*)(X+j+28),	mask[(b>>28)&0xf].xmm));
		}
		float	*r=(float*)&sum,	z=r[0]+r[1]+r[2]+r[3];
		testy[i]=binary?1.0f/(1.0f+expf(-z)):z;
	}	
}

void	RBR::estimate(const	char	*F){
	binary=true;
	double	sx=0,	sxx=0;
	for(size_t	i=0;	i<trainn;	i++){
		if(trainy[i]!=0&&trainy[i]!=1)	binary=false;
		sx+=trainy[i];	sxx+=(double)trainy[i]*trainy[i];
	}
	sx/=trainn;	sxx=sxx/trainn-sx*sx;	sxx=sxx>0?sqrt(sxx):1;
	if(!binary)	for(size_t	i=0;	i<trainn;	i++)	trainy[i]=(trainy[i]-sx)/sxx;
	penalty=binary?penalty*(bits-1):penalty*(bits-1);
	
	
	if(posix_memalign((void**)&dbeta,	16,	thread*bits*4)){	cerr<<"not enough memory\n";	return;	}	//	I
	lbfgsfloatval_t	score,	*x=lbfgs_malloc(bits);
	memset(x,	0,	bits*4);
	lbfgs_parameter_t	para;	lbfgs_parameter_init(&para);	
	para.m=correction;	para.max_iterations=0;	para.epsilon=0;	para.past=10;	para.delta=1e-4;
	
	cout<<"step\tloss\n";
	lbfgs(bits,	x,	&score,	evaluate,	progress,	this,	&para);
	cout<<endl;	

	free(dbeta);

	predict(x);
	
	weight.assign(feature,	0);
	for(size_t	i=1;	i<bits;	i++){
		float	w=(double)sum1[i]*(trainn+testn-sum1[i])/(trainn+testn)/(trainn+testn)*x[i]*x[i];
		for(size_t	j=0;	j<twist;	j++)	weight[fid[i*twist+j]]+=w;
	}
	lbfgs_free(x);
	
	ofstream	fo(F);
	for(size_t	i=0;	i<testn;	i++)	fo<<(binary?testy[i]:testy[i]*sxx+sx)<<'\n';
	fo.close();
	
	fo.open("variable.importance");
	for(size_t	i=0;	i<feature;	i++)	fo<<weight[i]<<'\n';
	fo.close();
}

void	RBR::document(void){
	cout<<"Usage:	rbr [options] trainX trainY testX testYhat\n";
	cout<<"\t-b:	number of random bits.		default=100000\n";
	cout<<"\t-r:	Scaled L2 regularization.	default=1\n";
	cout<<"\t-t:	number of features twisted.	default=2\n";
	cout<<"\t-c:	L-BFGS corrections.		default=128\n";
	cout<<"\t-T:	number of threads.		default=0 (all)\n";
	cout<<"\t-w:	variable weights file.\n";
	cout<<"\ttrainX,trainY,testX:	input files in csv/tsv format without header.\n";
	cout<<"\ttestYhat:	output files.\n";
	exit(0);
}

int	main(int	ac,	char	**av){
	size_t	t0=time(NULL);
	cout<<"***********************************\n";
	cout<<"* Random Bits Regression          *\n";
	cout<<"* author: Yi Wang                 *\n";
	cout<<"* email:  godspeed_china@yeah.net *\n";
	cout<<"* date:   27/Jun/2014             *\n";
	cout<<"***********************************\n";
	
	RBR	rbr;
	rbr.bits=100000;	rbr.penalty=1;	rbr.twist=2;	rbr.correction=128;	rbr.thread=0;
	int	opt;
	while((opt=getopt(ac,	av,	"b:r:t:c:T:w:"))>=0){
		switch(opt){
		case	'b':	rbr.bits=atoi(optarg);	break;
		case	'r':	rbr.penalty=atof(optarg);	break;
		case	't':	rbr.twist=atoi(optarg);	break;
		case	'c':	rbr.correction=atoi(optarg);	break;
		case	'T':	rbr.thread=atoi(optarg);	break;
		case	'w':	rbr.load_weight(optarg);	break;
		default:	rbr.document();
		}
	}
	if(ac<optind+4)	rbr.document();
	if(!rbr.thread)	rbr.thread=omp_get_num_procs();
	omp_set_num_threads(rbr.thread);
	
	vector<float>	trainx,	testx;
	size_t	trainn,	feature;
	if(!rbr.load_matrix(av[optind],	trainx,	rbr.trainn,	rbr.feature))	return	0;
	if(!rbr.load_matrix(av[optind+1],	rbr.trainy,	trainn,	feature))	return	0;
	if(trainn!=rbr.trainn){	cerr<<av[optind]<<'\t'<<av[optind+1]<<"\thas different number of rows\n";	return	0;	}
	if(feature!=1){	cerr<<av[optind+1]<<"\thas more than one column\n";	return	0;	}
	if(!rbr.load_matrix(av[optind+2],	testx,	rbr.testn,	feature))	return	0;
	if(feature!=rbr.feature){	cerr<<av[optind]<<'\t'<<av[optind+2]<<"\thas different number of columns\n";	return	0;	}
	
	rbr.x2bit(trainx,	testx);
	
	rbr.estimate(av[optind+3]);
	cout<<"time\t"<<time(NULL)-t0+1<<"s\n\n";	
	return	0;
}
