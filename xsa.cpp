#ifndef	XSA_INCLUDED
#define	XSA_INCLUDED
#include	<stdint.h>
#include	<math.h>
union	XSA_union_double{	uint64_t	i;	double	f;	};
union	XSA_union_single{	uint32_t	i;	float	f;	};
class	XSA{
private:
	uint64_t	s[2];
	double	uniform_positive(void){	XSA_union_double	u;	do	u.i=(get()&0xfffffffffffffull)|0x3ff0000000000000ull;	while(u.f==1.0);	return	u.f-1.0;	}
public:
	void	set(uint64_t	S){	s[0]=S*0x5DEECE66Dull+1013904223ull;	s[1]=0x1234567890ABCDEFull;	}
	uint64_t	get(void){	uint64_t	s1=s[0],	s0=s[1];	s[0]=s0;	s1^=s1<<23;	return	(s[1]=(s1^s0^(s1>>17)^(s0>>26)))+s0;	}
	double	uniform(void){	XSA_union_double	u;	u.i=(get()&0xfffffffffffffull)|0x3ff0000000000000ull;	return	u.f-1.0;	}	
	float	uniform_single(void){	XSA_union_single	u;	u.i=(get()&0x7ffffful)|0x3f800000ul;	return	u.f-1.0f;	}
	double	normal(void){	double	x,	y,	r2;	do{	x=-1.0+2.0*uniform_positive();	y=-1.0+2.0*uniform_positive();	r2=x*x+y*y;	}while(r2>1.0||r2==0);	return	y*sqrt(-2.0*log(r2)/r2);	}
	double	fast_normal(void){	return	2.0*(uniform()+uniform()+uniform()-1.5);	}
};
#endif

