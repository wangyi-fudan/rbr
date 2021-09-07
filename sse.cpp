/*
	Name:	SSE vectorization library
	Ver:	0.1
	Author:	Wangyi
	Date:	14-12-13
*/

#ifndef	SSE_Included
#define	SSE_Included

#include	<stdlib.h>
typedef	float	v4sf	__attribute__	((__vector_size__	(16)));
typedef	double	v2df	__attribute__	((__vector_size__	(16)));
const	size_t	L2Cache=512*1024;

//	aligned size of X
inline	size_t	a_size_f(size_t	X){
	return	X%4?(X/4+1)*4:X;
}

inline	size_t	a_size_d(size_t	X){
	return	X%2?(X/2+1)*2:X;
}

//	A.B
inline	double	vec_vec_dot_d(const	double	*A,	const	double	*B,	size_t	N){
	v2df	s={0,0};	double	*r=(double*)&s;
	if(N<8)	switch(N){
		case	0:	return	0;
		case	1:	return	A[0]*B[0];
		case	2:	s=__builtin_ia32_addpd(s,	__builtin_ia32_mulpd(*(v2df*)(A),	*(v2df*)(B)));
					return	r[0]+r[1];
		case	3:	s=__builtin_ia32_addpd(s,	__builtin_ia32_mulpd(*(v2df*)(A),	*(v2df*)(B)));
					return	r[0]+r[1]+A[2]*B[2];
		case	4:	s=__builtin_ia32_addpd(s,	__builtin_ia32_mulpd(*(v2df*)(A),	*(v2df*)(B)));
					s=__builtin_ia32_addpd(s,	__builtin_ia32_mulpd(*(v2df*)(A+2),	*(v2df*)(B+2)));
					return	r[0]+r[1];
		case	5:	s=__builtin_ia32_addpd(s,	__builtin_ia32_mulpd(*(v2df*)(A),	*(v2df*)(B)));
					s=__builtin_ia32_addpd(s,	__builtin_ia32_mulpd(*(v2df*)(A+2),	*(v2df*)(B+2)));
					return	r[0]+r[1]+A[4]*B[4];
		case	6:	s=__builtin_ia32_addpd(s,	__builtin_ia32_mulpd(*(v2df*)(A),	*(v2df*)(B)));
					s=__builtin_ia32_addpd(s,	__builtin_ia32_mulpd(*(v2df*)(A+2),	*(v2df*)(B+2)));
					s=__builtin_ia32_addpd(s,	__builtin_ia32_mulpd(*(v2df*)(A+4),	*(v2df*)(B+4)));
					return	r[0]+r[1];
		case	7:	s=__builtin_ia32_addpd(s,	__builtin_ia32_mulpd(*(v2df*)(A),	*(v2df*)(B)));
					s=__builtin_ia32_addpd(s,	__builtin_ia32_mulpd(*(v2df*)(A+2),	*(v2df*)(B+2)));
					s=__builtin_ia32_addpd(s,	__builtin_ia32_mulpd(*(v2df*)(A+4),	*(v2df*)(B+4)));
					return	r[0]+r[1]+A[6]*B[6];
	}
	size_t	max=(N>>3)<<3;
	for(size_t	i=0;	i<max;	i+=8){
		s=__builtin_ia32_addpd(s,__builtin_ia32_addpd(
		__builtin_ia32_mulpd(*(v2df*)(A+i),	*(v2df*)(B+i)),	
		__builtin_ia32_mulpd(*(v2df*)(A+i+2),	*(v2df*)(B+i+2))));
		s=__builtin_ia32_addpd(s,__builtin_ia32_addpd(
		__builtin_ia32_mulpd(*(v2df*)(A+i+4),	*(v2df*)(B+i+4)),	
		__builtin_ia32_mulpd(*(v2df*)(A+i+6),	*(v2df*)(B+i+6))));
	}
	A+=max;	B+=max;
	switch(N&7){
		case	0:	return	r[0]+r[1];
		case	1:	return	r[0]+r[1]+A[0]*B[0];
		case	2:	s=__builtin_ia32_addpd(s,	__builtin_ia32_mulpd(*(v2df*)(A),	*(v2df*)(B)));
					return	r[0]+r[1];
		case	3:	s=__builtin_ia32_addpd(s,	__builtin_ia32_mulpd(*(v2df*)(A),	*(v2df*)(B)));
					return	r[0]+r[1]+A[2]*B[2];
		case	4:	s=__builtin_ia32_addpd(s,	__builtin_ia32_mulpd(*(v2df*)(A),	*(v2df*)(B)));
					s=__builtin_ia32_addpd(s,	__builtin_ia32_mulpd(*(v2df*)(A+2),	*(v2df*)(B+2)));
					return	r[0]+r[1];
		case	5:	s=__builtin_ia32_addpd(s,	__builtin_ia32_mulpd(*(v2df*)(A),	*(v2df*)(B)));
					s=__builtin_ia32_addpd(s,	__builtin_ia32_mulpd(*(v2df*)(A+2),	*(v2df*)(B+2)));
					return	r[0]+r[1]+A[4]*B[4];
		case	6:	s=__builtin_ia32_addpd(s,	__builtin_ia32_mulpd(*(v2df*)(A),	*(v2df*)(B)));
					s=__builtin_ia32_addpd(s,	__builtin_ia32_mulpd(*(v2df*)(A+2),	*(v2df*)(B+2)));
					s=__builtin_ia32_addpd(s,	__builtin_ia32_mulpd(*(v2df*)(A+4),	*(v2df*)(B+4)));
					return	r[0]+r[1];
		case	7:	s=__builtin_ia32_addpd(s,	__builtin_ia32_mulpd(*(v2df*)(A),	*(v2df*)(B)));
					s=__builtin_ia32_addpd(s,	__builtin_ia32_mulpd(*(v2df*)(A+2),	*(v2df*)(B+2)));
					s=__builtin_ia32_addpd(s,	__builtin_ia32_mulpd(*(v2df*)(A+4),	*(v2df*)(B+4)));
					return	r[0]+r[1]+A[6]*B[6];
	}
	return	0;
}

inline	float	vec_vec_dot_f(const	float	*A,	const	float	*B,	size_t	N){
	v4sf	s={0,0,0,0};	float	*r=(float*)&s;
	if(N<16)	switch(N){
		case	0:	return	0;
		case	1:	return	A[0]*B[0];
		case	2:	return	A[0]*B[0]+A[1]*B[1];
		case	3:	return	A[0]*B[0]+A[1]*B[1]+A[2]*B[2];
		case	4:	s=__builtin_ia32_addps(s,	__builtin_ia32_mulps(*(v4sf*)(A),	*(v4sf*)(B)));
					return	r[0]+r[1]+r[2]+r[3];
		case	5:	s=__builtin_ia32_addps(s,	__builtin_ia32_mulps(*(v4sf*)(A),	*(v4sf*)(B)));
					return	r[0]+r[1]+r[2]+r[3]+A[4]*B[4];
		case	6:	s=__builtin_ia32_addps(s,	__builtin_ia32_mulps(*(v4sf*)(A),	*(v4sf*)(B)));
					return	r[0]+r[1]+r[2]+r[3]+A[4]*B[4]+A[5]*B[5];
		case	7:	s=__builtin_ia32_addps(s,	__builtin_ia32_mulps(*(v4sf*)(A),	*(v4sf*)(B)));
					return	r[0]+r[1]+r[2]+r[3]+A[4]*B[4]+A[5]*B[5]+A[6]*B[6];
		case	8:	s=__builtin_ia32_addps(s,	__builtin_ia32_mulps(*(v4sf*)(A),	*(v4sf*)(B)));
					s=__builtin_ia32_addps(s,	__builtin_ia32_mulps(*(v4sf*)(A+4),	*(v4sf*)(B+4)));
					return	r[0]+r[1]+r[2]+r[3];
		case	9:	s=__builtin_ia32_addps(s,	__builtin_ia32_mulps(*(v4sf*)(A),	*(v4sf*)(B)));
					s=__builtin_ia32_addps(s,	__builtin_ia32_mulps(*(v4sf*)(A+4),	*(v4sf*)(B+4)));
					return	r[0]+r[1]+r[2]+r[3]+A[8]*B[8];
		case	10:	s=__builtin_ia32_addps(s,	__builtin_ia32_mulps(*(v4sf*)(A),	*(v4sf*)(B)));
					s=__builtin_ia32_addps(s,	__builtin_ia32_mulps(*(v4sf*)(A+4),	*(v4sf*)(B+4)));		
					return	r[0]+r[1]+r[2]+r[3]+A[8]*B[8]+A[9]*B[9];
		case	11:	s=__builtin_ia32_addps(s,	__builtin_ia32_mulps(*(v4sf*)(A),	*(v4sf*)(B)));
					s=__builtin_ia32_addps(s,	__builtin_ia32_mulps(*(v4sf*)(A+4),	*(v4sf*)(B+4)));		
					return	r[0]+r[1]+r[2]+r[3]+A[8]*B[8]+A[9]*B[9]+A[10]*B[10];	
		case	12:	s=__builtin_ia32_addps(s,	__builtin_ia32_mulps(*(v4sf*)(A),	*(v4sf*)(B)));
					s=__builtin_ia32_addps(s,	__builtin_ia32_mulps(*(v4sf*)(A+4),	*(v4sf*)(B+4)));
					s=__builtin_ia32_addps(s,	__builtin_ia32_mulps(*(v4sf*)(A+8),	*(v4sf*)(B+8)));
					return	r[0]+r[1]+r[2]+r[3];
		case	13:	s=__builtin_ia32_addps(s,	__builtin_ia32_mulps(*(v4sf*)(A),	*(v4sf*)(B)));
					s=__builtin_ia32_addps(s,	__builtin_ia32_mulps(*(v4sf*)(A+4),	*(v4sf*)(B+4)));
					s=__builtin_ia32_addps(s,	__builtin_ia32_mulps(*(v4sf*)(A+8),	*(v4sf*)(B+8)));
					return	r[0]+r[1]+r[2]+r[3]+A[12]*B[12];
		case	14:	s=__builtin_ia32_addps(s,	__builtin_ia32_mulps(*(v4sf*)(A),	*(v4sf*)(B)));
					s=__builtin_ia32_addps(s,	__builtin_ia32_mulps(*(v4sf*)(A+4),	*(v4sf*)(B+4)));
					s=__builtin_ia32_addps(s,	__builtin_ia32_mulps(*(v4sf*)(A+8),	*(v4sf*)(B+8)));					
					return	r[0]+r[1]+r[2]+r[3]+A[12]*B[12]+A[13]*B[13];
		case	15:	s=__builtin_ia32_addps(s,	__builtin_ia32_mulps(*(v4sf*)(A),	*(v4sf*)(B)));
					s=__builtin_ia32_addps(s,	__builtin_ia32_mulps(*(v4sf*)(A+4),	*(v4sf*)(B+4)));
					s=__builtin_ia32_addps(s,	__builtin_ia32_mulps(*(v4sf*)(A+8),	*(v4sf*)(B+8)));					
					return	r[0]+r[1]+r[2]+r[3]+A[12]*B[12]+A[13]*B[13]+A[14]*B[14];
	}
	
	size_t	max=(N>>4)<<4;
	for(size_t	i=0;	i<max;	i+=16){
		s=__builtin_ia32_addps(s,__builtin_ia32_addps(
		__builtin_ia32_mulps(*(v4sf*)(A+i),	*(v4sf*)(B+i)),	
		__builtin_ia32_mulps(*(v4sf*)(A+i+4),	*(v4sf*)(B+i+4))));
		s=__builtin_ia32_addps(s,__builtin_ia32_addps(
		__builtin_ia32_mulps(*(v4sf*)(A+i+8),	*(v4sf*)(B+i+8)),
		__builtin_ia32_mulps(*(v4sf*)(A+i+12),	*(v4sf*)(B+i+12))));
	}
	A+=max;	B+=max;
	switch(N&15){
		case	0:	return	r[0]+r[1]+r[2]+r[3];
		case	1:	return	r[0]+r[1]+r[2]+r[3]+A[0]*B[0];
		case	2:	return	r[0]+r[1]+r[2]+r[3]+A[0]*B[0]+A[1]*B[1];
		case	3:	return	r[0]+r[1]+r[2]+r[3]+A[0]*B[0]+A[1]*B[1]+A[2]*B[2];
		case	4:	s=__builtin_ia32_addps(s,	__builtin_ia32_mulps(*(v4sf*)(A),	*(v4sf*)(B)));
					return	r[0]+r[1]+r[2]+r[3];
		case	5:	s=__builtin_ia32_addps(s,	__builtin_ia32_mulps(*(v4sf*)(A),	*(v4sf*)(B)));
					return	r[0]+r[1]+r[2]+r[3]+A[4]*B[4];
		case	6:	s=__builtin_ia32_addps(s,	__builtin_ia32_mulps(*(v4sf*)(A),	*(v4sf*)(B)));
					return	r[0]+r[1]+r[2]+r[3]+A[4]*B[4]+A[5]*B[5];
		case	7:	s=__builtin_ia32_addps(s,	__builtin_ia32_mulps(*(v4sf*)(A),	*(v4sf*)(B)));
					return	r[0]+r[1]+r[2]+r[3]+A[4]*B[4]+A[5]*B[5]+A[6]*B[6];
		case	8:	s=__builtin_ia32_addps(s,	__builtin_ia32_mulps(*(v4sf*)(A),	*(v4sf*)(B)));
					s=__builtin_ia32_addps(s,	__builtin_ia32_mulps(*(v4sf*)(A+4),	*(v4sf*)(B+4)));
					return	r[0]+r[1]+r[2]+r[3];
		case	9:	s=__builtin_ia32_addps(s,	__builtin_ia32_mulps(*(v4sf*)(A),	*(v4sf*)(B)));
					s=__builtin_ia32_addps(s,	__builtin_ia32_mulps(*(v4sf*)(A+4),	*(v4sf*)(B+4)));
					return	r[0]+r[1]+r[2]+r[3]+A[8]*B[8];
		case	10:	s=__builtin_ia32_addps(s,	__builtin_ia32_mulps(*(v4sf*)(A),	*(v4sf*)(B)));
					s=__builtin_ia32_addps(s,	__builtin_ia32_mulps(*(v4sf*)(A+4),	*(v4sf*)(B+4)));		
					return	r[0]+r[1]+r[2]+r[3]+A[8]*B[8]+A[9]*B[9];
		case	11:	s=__builtin_ia32_addps(s,	__builtin_ia32_mulps(*(v4sf*)(A),	*(v4sf*)(B)));
					s=__builtin_ia32_addps(s,	__builtin_ia32_mulps(*(v4sf*)(A+4),	*(v4sf*)(B+4)));		
					return	r[0]+r[1]+r[2]+r[3]+A[8]*B[8]+A[9]*B[9]+A[10]*B[10];	
		case	12:	s=__builtin_ia32_addps(s,	__builtin_ia32_mulps(*(v4sf*)(A),	*(v4sf*)(B)));
					s=__builtin_ia32_addps(s,	__builtin_ia32_mulps(*(v4sf*)(A+4),	*(v4sf*)(B+4)));
					s=__builtin_ia32_addps(s,	__builtin_ia32_mulps(*(v4sf*)(A+8),	*(v4sf*)(B+8)));
					return	r[0]+r[1]+r[2]+r[3];
		case	13:	s=__builtin_ia32_addps(s,	__builtin_ia32_mulps(*(v4sf*)(A),	*(v4sf*)(B)));
					s=__builtin_ia32_addps(s,	__builtin_ia32_mulps(*(v4sf*)(A+4),	*(v4sf*)(B+4)));
					s=__builtin_ia32_addps(s,	__builtin_ia32_mulps(*(v4sf*)(A+8),	*(v4sf*)(B+8)));
					return	r[0]+r[1]+r[2]+r[3]+A[12]*B[12];
		case	14:	s=__builtin_ia32_addps(s,	__builtin_ia32_mulps(*(v4sf*)(A),	*(v4sf*)(B)));
					s=__builtin_ia32_addps(s,	__builtin_ia32_mulps(*(v4sf*)(A+4),	*(v4sf*)(B+4)));
					s=__builtin_ia32_addps(s,	__builtin_ia32_mulps(*(v4sf*)(A+8),	*(v4sf*)(B+8)));					
					return	r[0]+r[1]+r[2]+r[3]+A[12]*B[12]+A[13]*B[13];
		case	15:	s=__builtin_ia32_addps(s,	__builtin_ia32_mulps(*(v4sf*)(A),	*(v4sf*)(B)));
					s=__builtin_ia32_addps(s,	__builtin_ia32_mulps(*(v4sf*)(A+4),	*(v4sf*)(B+4)));
					s=__builtin_ia32_addps(s,	__builtin_ia32_mulps(*(v4sf*)(A+8),	*(v4sf*)(B+8)));					
					return	r[0]+r[1]+r[2]+r[3]+A[12]*B[12]+A[13]*B[13]+A[14]*B[14];
	}
	return	0;
}

//	A+=Z*B
inline	void	vec_con_fma_d(double	*A,	const	double	*B,	float	Z,	size_t	N){
	v2df	s={Z,Z};
	size_t	max=(N>>3)<<3;
	for(size_t	i=0;	i<max;	i+=8){
		*(v2df*)(A+i)=__builtin_ia32_addpd(*(v2df*)(A+i),	__builtin_ia32_mulpd(s,	*(v2df*)(B+i)));
		*(v2df*)(A+i+2)=__builtin_ia32_addpd(*(v2df*)(A+i+2),	__builtin_ia32_mulpd(s,	*(v2df*)(B+i+2)));
		*(v2df*)(A+i+4)=__builtin_ia32_addpd(*(v2df*)(A+i+4),	__builtin_ia32_mulpd(s,	*(v2df*)(B+i+4)));
		*(v2df*)(A+i+6)=__builtin_ia32_addpd(*(v2df*)(A+i+6),	__builtin_ia32_mulpd(s,	*(v2df*)(B+i+6)));
	}
	A+=max;	B+=max;
	switch(N&7){
		case	1:	A[0]+=Z*B[0];	break;
		case	2:	*(v2df*)(A)=__builtin_ia32_addpd(*(v2df*)(A),	__builtin_ia32_mulpd(s,	*(v2df*)(B)));
					break;
		case	3:	*(v2df*)(A)=__builtin_ia32_addpd(*(v2df*)(A),	__builtin_ia32_mulpd(s,	*(v2df*)(B)));
					A[2]+=Z*B[2];	break;
		case	4:	*(v2df*)(A)=__builtin_ia32_addpd(*(v2df*)(A),	__builtin_ia32_mulpd(s,	*(v2df*)(B)));
					*(v2df*)(A+2)=__builtin_ia32_addpd(*(v2df*)(A+2),	__builtin_ia32_mulpd(s,	*(v2df*)(B+2)));	
					break;
		case	5:	*(v2df*)(A)=__builtin_ia32_addpd(*(v2df*)(A),	__builtin_ia32_mulpd(s,	*(v2df*)(B)));
					*(v2df*)(A+2)=__builtin_ia32_addpd(*(v2df*)(A+2),	__builtin_ia32_mulpd(s,	*(v2df*)(B+2)));
					A[4]+=Z*B[4];	break;
		case	6:	*(v2df*)(A)=__builtin_ia32_addpd(*(v2df*)(A),	__builtin_ia32_mulpd(s,	*(v2df*)(B)));
					*(v2df*)(A+2)=__builtin_ia32_addpd(*(v2df*)(A+2),	__builtin_ia32_mulpd(s,	*(v2df*)(B+2)));
					*(v2df*)(A+4)=__builtin_ia32_addpd(*(v2df*)(A+4),	__builtin_ia32_mulpd(s,	*(v2df*)(B+4)));
					break;
		case	7:	*(v2df*)(A)=__builtin_ia32_addpd(*(v2df*)(A),	__builtin_ia32_mulpd(s,	*(v2df*)(B)));
					*(v2df*)(A+2)=__builtin_ia32_addpd(*(v2df*)(A+2),	__builtin_ia32_mulpd(s,	*(v2df*)(B+2)));
					*(v2df*)(A+4)=__builtin_ia32_addpd(*(v2df*)(A+4),	__builtin_ia32_mulpd(s,	*(v2df*)(B+4)));
					A[6]+=Z*B[6];	break;
	}
}

inline	void	vec_con_fma_f(float	*A,	const	float	*B,	float	Z,	size_t	N){
	size_t	max=(N>>4)<<4;
	v4sf	s={Z,Z,Z,Z};
	for(size_t	i=0;	i<max;	i+=16){
		*(v4sf*)(A+i)=__builtin_ia32_addps(*(v4sf*)(A+i),	__builtin_ia32_mulps(s,	*(v4sf*)(B+i)));
		*(v4sf*)(A+i+4)=__builtin_ia32_addps(*(v4sf*)(A+i+4),	__builtin_ia32_mulps(s,	*(v4sf*)(B+i+4)));
		*(v4sf*)(A+i+8)=__builtin_ia32_addps(*(v4sf*)(A+i+8),	__builtin_ia32_mulps(s,	*(v4sf*)(B+i+8)));
		*(v4sf*)(A+i+12)=__builtin_ia32_addps(*(v4sf*)(A+i+12),	__builtin_ia32_mulps(s,	*(v4sf*)(B+i+12)));
	}
	A+=max;	B+=max;
	switch(N&15){
		case	1:	A[0]+=Z*B[0];	break;
		case	2:	A[0]+=Z*B[0];	A[1]+=Z*B[1];	break;
		case	3:	A[0]+=Z*B[0];	A[1]+=Z*B[1];	A[2]+=Z*B[2];	break;
		case	4:	*(v4sf*)A=__builtin_ia32_addps(*(v4sf*)A,	__builtin_ia32_mulps(s,	*(v4sf*)B));	
					break;
		case	5:	*(v4sf*)A=__builtin_ia32_addps(*(v4sf*)A,	__builtin_ia32_mulps(s,	*(v4sf*)B));	
					A[4]+=Z*B[4];	break;
		case	6:	*(v4sf*)A=__builtin_ia32_addps(*(v4sf*)A,	__builtin_ia32_mulps(s,	*(v4sf*)B));	
					A[4]+=Z*B[4];	A[5]+=Z*B[5];	break;
		case	7:	*(v4sf*)A=__builtin_ia32_addps(*(v4sf*)A,	__builtin_ia32_mulps(s,	*(v4sf*)B));	
					A[4]+=Z*B[4];	A[5]+=Z*B[5];	A[6]+=Z*B[6];	break;
		case	8:	*(v4sf*)A=__builtin_ia32_addps(*(v4sf*)A,	__builtin_ia32_mulps(s,	*(v4sf*)B));
					*(v4sf*)(A+4)=__builtin_ia32_addps(*(v4sf*)(A+4),	__builtin_ia32_mulps(s,	*(v4sf*)(B+4)));
					break;
		case	9:	*(v4sf*)A=__builtin_ia32_addps(*(v4sf*)A,	__builtin_ia32_mulps(s,	*(v4sf*)B));
					*(v4sf*)(A+4)=__builtin_ia32_addps(*(v4sf*)(A+4),	__builtin_ia32_mulps(s,	*(v4sf*)(B+4)));
					A[8]+=Z*B[8];	break;
		case	10:	*(v4sf*)A=__builtin_ia32_addps(*(v4sf*)A,	__builtin_ia32_mulps(s,	*(v4sf*)B));
					*(v4sf*)(A+4)=__builtin_ia32_addps(*(v4sf*)(A+4),	__builtin_ia32_mulps(s,	*(v4sf*)(B+4)));
					A[8]+=Z*B[8];	A[9]+=Z*B[9];	break;
		case	11:	*(v4sf*)A=__builtin_ia32_addps(*(v4sf*)A,	__builtin_ia32_mulps(s,	*(v4sf*)B));
					*(v4sf*)(A+4)=__builtin_ia32_addps(*(v4sf*)(A+4),	__builtin_ia32_mulps(s,	*(v4sf*)(B+4)));
					A[8]+=Z*B[8];	A[9]+=Z*B[9];	A[10]+=Z*B[10];	break;
		case	12:	*(v4sf*)A=__builtin_ia32_addps(*(v4sf*)A,	__builtin_ia32_mulps(s,	*(v4sf*)B));
					*(v4sf*)(A+4)=__builtin_ia32_addps(*(v4sf*)(A+4),	__builtin_ia32_mulps(s,	*(v4sf*)(B+4)));
					*(v4sf*)(A+8)=__builtin_ia32_addps(*(v4sf*)(A+8),	__builtin_ia32_mulps(s,	*(v4sf*)(B+8)));
					break;
		case	13:	*(v4sf*)A=__builtin_ia32_addps(*(v4sf*)A,	__builtin_ia32_mulps(s,	*(v4sf*)B));
					*(v4sf*)(A+4)=__builtin_ia32_addps(*(v4sf*)(A+4),	__builtin_ia32_mulps(s,	*(v4sf*)(B+4)));
					*(v4sf*)(A+8)=__builtin_ia32_addps(*(v4sf*)(A+8),	__builtin_ia32_mulps(s,	*(v4sf*)(B+8)));					
					A[12]+=Z*B[12];	break;
		case	14:	*(v4sf*)A=__builtin_ia32_addps(*(v4sf*)A,	__builtin_ia32_mulps(s,	*(v4sf*)B));
					*(v4sf*)(A+4)=__builtin_ia32_addps(*(v4sf*)(A+4),	__builtin_ia32_mulps(s,	*(v4sf*)(B+4)));
					*(v4sf*)(A+8)=__builtin_ia32_addps(*(v4sf*)(A+8),	__builtin_ia32_mulps(s,	*(v4sf*)(B+8)));					
					A[12]+=Z*B[12];	A[13]+=Z*B[13];	break;
		case	15:	*(v4sf*)A=__builtin_ia32_addps(*(v4sf*)A,	__builtin_ia32_mulps(s,	*(v4sf*)B));
					*(v4sf*)(A+4)=__builtin_ia32_addps(*(v4sf*)(A+4),	__builtin_ia32_mulps(s,	*(v4sf*)(B+4)));
					*(v4sf*)(A+8)=__builtin_ia32_addps(*(v4sf*)(A+8),	__builtin_ia32_mulps(s,	*(v4sf*)(B+8)));					
					A[12]+=Z*B[12];	A[13]+=Z*B[13];	A[14]+=Z*B[14];	break;
	}
}

//	3X3 matrix multiply
inline	void	mat_mul_3x3(
	const	float	*A0,	const	float	*A1,	const	float	*A2,
	const	float	*B0,	const	float	*B1,	const	float	*B2,
	float	*S00,	float	*S01,	float	*S02,
	float	*S10,	float	*S11,	float	*S12,
	float	*S20,	float	*S21,	float	*S22,
	size_t	N){
	
	size_t	max=(N>>2)<<2;

	v4sf	s00={0,0,0,0},	s01={0,0,0,0},	s02={0,0,0,0};
	v4sf	s10={0,0,0,0},	s11={0,0,0,0},	s12={0,0,0,0};
	v4sf	s20={0,0,0,0},	s21={0,0,0,0},	s22={0,0,0,0};
	for(size_t	k=0;	k<max;	k+=4){
		s00=__builtin_ia32_addps(s00,	__builtin_ia32_mulps(*(v4sf*)(A0+k),	*(v4sf*)(B0+k)));
		s01=__builtin_ia32_addps(s01,	__builtin_ia32_mulps(*(v4sf*)(A0+k),	*(v4sf*)(B1+k)));
		s02=__builtin_ia32_addps(s02,	__builtin_ia32_mulps(*(v4sf*)(A0+k),	*(v4sf*)(B2+k)));
		s10=__builtin_ia32_addps(s10,	__builtin_ia32_mulps(*(v4sf*)(A1+k),	*(v4sf*)(B0+k)));
		s11=__builtin_ia32_addps(s11,	__builtin_ia32_mulps(*(v4sf*)(A1+k),	*(v4sf*)(B1+k)));
		s12=__builtin_ia32_addps(s12,	__builtin_ia32_mulps(*(v4sf*)(A1+k),	*(v4sf*)(B2+k)));
		s20=__builtin_ia32_addps(s20,	__builtin_ia32_mulps(*(v4sf*)(A2+k),	*(v4sf*)(B0+k)));
		s21=__builtin_ia32_addps(s21,	__builtin_ia32_mulps(*(v4sf*)(A2+k),	*(v4sf*)(B1+k)));
		s22=__builtin_ia32_addps(s22,	__builtin_ia32_mulps(*(v4sf*)(A2+k),	*(v4sf*)(B2+k)));
	}
	A0+=max;	A1+=max;	A2+=max;	B0+=max;	B1+=max;	B2+=max;
	float	*r;
	switch(N&3){
	case	0:
		r=(float*)&s00;	*S00=r[0]+r[1]+r[2]+r[3];
		r=(float*)&s01;	*S01=r[0]+r[1]+r[2]+r[3];
		r=(float*)&s02;	*S02=r[0]+r[1]+r[2]+r[3];
		r=(float*)&s10;	*S10=r[0]+r[1]+r[2]+r[3];
		r=(float*)&s11;	*S11=r[0]+r[1]+r[2]+r[3];
		r=(float*)&s12;	*S12=r[0]+r[1]+r[2]+r[3];
		r=(float*)&s20;	*S20=r[0]+r[1]+r[2]+r[3];
		r=(float*)&s21;	*S21=r[0]+r[1]+r[2]+r[3];
		r=(float*)&s22;	*S22=r[0]+r[1]+r[2]+r[3];
		break;
	case	1:
		r=(float*)&s00;	*S00=r[0]+r[1]+r[2]+r[3]+A0[0]*B0[0];
		r=(float*)&s01;	*S01=r[0]+r[1]+r[2]+r[3]+A0[0]*B1[0];
		r=(float*)&s02;	*S02=r[0]+r[1]+r[2]+r[3]+A0[0]*B2[0];
		r=(float*)&s10;	*S10=r[0]+r[1]+r[2]+r[3]+A1[0]*B0[0];
		r=(float*)&s11;	*S11=r[0]+r[1]+r[2]+r[3]+A1[0]*B1[0];
		r=(float*)&s12;	*S12=r[0]+r[1]+r[2]+r[3]+A1[0]*B2[0];
		r=(float*)&s20;	*S20=r[0]+r[1]+r[2]+r[3]+A2[0]*B0[0];
		r=(float*)&s21;	*S21=r[0]+r[1]+r[2]+r[3]+A2[0]*B1[0];
		r=(float*)&s22;	*S22=r[0]+r[1]+r[2]+r[3]+A2[0]*B2[0];
		break;
	case	2:
		r=(float*)&s00;	*S00=r[0]+r[1]+r[2]+r[3]+A0[0]*B0[0]+A0[1]*B0[1];
		r=(float*)&s01;	*S01=r[0]+r[1]+r[2]+r[3]+A0[0]*B1[0]+A0[1]*B1[1];
		r=(float*)&s02;	*S02=r[0]+r[1]+r[2]+r[3]+A0[0]*B2[0]+A0[1]*B2[1];
		r=(float*)&s10;	*S10=r[0]+r[1]+r[2]+r[3]+A1[0]*B0[0]+A1[1]*B0[1];
		r=(float*)&s11;	*S11=r[0]+r[1]+r[2]+r[3]+A1[0]*B1[0]+A1[1]*B1[1];
		r=(float*)&s12;	*S12=r[0]+r[1]+r[2]+r[3]+A1[0]*B2[0]+A1[1]*B2[1];
		r=(float*)&s20;	*S20=r[0]+r[1]+r[2]+r[3]+A2[0]*B0[0]+A2[1]*B0[1];
		r=(float*)&s21;	*S21=r[0]+r[1]+r[2]+r[3]+A2[0]*B1[0]+A2[1]*B1[1];
		r=(float*)&s22;	*S22=r[0]+r[1]+r[2]+r[3]+A2[0]*B2[0]+A2[1]*B2[1];
		break;
	case	3:
		r=(float*)&s00;	*S00=r[0]+r[1]+r[2]+r[3]+A0[0]*B0[0]+A0[1]*B0[1]+A0[2]*B0[2];
		r=(float*)&s01;	*S01=r[0]+r[1]+r[2]+r[3]+A0[0]*B1[0]+A0[1]*B1[1]+A0[2]*B1[2];
		r=(float*)&s02;	*S02=r[0]+r[1]+r[2]+r[3]+A0[0]*B2[0]+A0[1]*B2[1]+A0[2]*B2[2];
		r=(float*)&s10;	*S10=r[0]+r[1]+r[2]+r[3]+A1[0]*B0[0]+A1[1]*B0[1]+A1[2]*B0[2];
		r=(float*)&s11;	*S11=r[0]+r[1]+r[2]+r[3]+A1[0]*B1[0]+A1[1]*B1[1]+A1[2]*B1[2];
		r=(float*)&s12;	*S12=r[0]+r[1]+r[2]+r[3]+A1[0]*B2[0]+A1[1]*B2[1]+A1[2]*B2[2];
		r=(float*)&s20;	*S20=r[0]+r[1]+r[2]+r[3]+A2[0]*B0[0]+A2[1]*B0[1]+A2[2]*B0[2];
		r=(float*)&s21;	*S21=r[0]+r[1]+r[2]+r[3]+A2[0]*B1[0]+A2[1]*B1[1]+A2[2]*B1[2];
		r=(float*)&s22;	*S22=r[0]+r[1]+r[2]+r[3]+A2[0]*B2[0]+A2[1]*B2[1]+A2[2]*B2[2];
		break;
	}
}

//	C=A*B
inline	void	mat_mul(const	float	*A,	const	float	*B,	float	*C,	size_t	M,	size_t	N,	size_t	K){
	if(!M||!N||!K)	return;
	size_t	str=a_size_f(K),	str1=a_size_f(N),	cache=(L2Cache/(str*4*4))*3?(L2Cache/(str*4*4))*3:3;
	size_t	maxa=(M/3)*3,	maxb=(N/3)*3,	maxc=(N/cache)*cache;
	
	if(M*N*K<40000){
	for(size_t	c=0;	c<maxc;	c+=cache){
		for(size_t	i=0;	i<maxa;	i+=3)	for(size_t	j=0;	j<cache;	j+=3)
			mat_mul_3x3(
				A+i*str,	A+(i+1)*str,	A+(i+2)*str,
				B+(c+j)*str,	B+(c+j+1)*str,	B+(c+j+2)*str,
				C+i*str1+c+j,	C+i*str1+c+j+1,	C+i*str1+c+j+2,
				C+(i+1)*str1+c+j,	C+(i+1)*str1+c+j+1,	C+(i+1)*str1+c+j+2,
				C+(i+2)*str1+c+j,	C+(i+2)*str1+c+j+1,	C+(i+2)*str1+c+j+2,
				K
			);
		for(size_t	i=maxa;	i<M;	i++)	for(size_t	j=0;	j<cache;	j++)
			C[i*str1+c+j]=vec_vec_dot_f(A+i*str,	B+(c+j)*str,	K);
	}
	for(size_t	i=0;	i<maxa;	i+=3)	for(size_t	j=maxc;	j<maxb;	j+=3)	
		mat_mul_3x3(
			A+i*str,	A+(i+1)*str,	A+(i+2)*str,
			B+j*str,	B+(j+1)*str,	B+(j+2)*str,
			C+i*str1+j,	C+i*str1+j+1,	C+i*str1+j+2,
			C+(i+1)*str1+j,	C+(i+1)*str1+j+1,	C+(i+1)*str1+j+2,
			C+(i+2)*str1+j,	C+(i+2)*str1+j+1,	C+(i+2)*str1+j+2,
			K
		);
	for(size_t	j=maxc;	j<maxb;	j++)	for(size_t	i=maxa;	i<M;	i++)
		C[i*str1+j]=vec_vec_dot_f(A+i*str,	B+j*str,	K);
	for(size_t	i=0;	i<M;	i++)	for(size_t	j=maxb;	j<N;	j++)	
		C[i*str1+j]=vec_vec_dot_f(A+i*str,	B+j*str,	K);
		return;
	}

	for(size_t	c=0;	c<maxc;	c+=cache){
		#pragma omp parallel for
		for(size_t	i=0;	i<maxa;	i+=3)	for(size_t	j=0;	j<cache;	j+=3)
			mat_mul_3x3(
				A+i*str,	A+(i+1)*str,	A+(i+2)*str,
				B+(c+j)*str,	B+(c+j+1)*str,	B+(c+j+2)*str,
				C+i*str1+c+j,	C+i*str1+c+j+1,	C+i*str1+c+j+2,
				C+(i+1)*str1+c+j,	C+(i+1)*str1+c+j+1,	C+(i+1)*str1+c+j+2,
				C+(i+2)*str1+c+j,	C+(i+2)*str1+c+j+1,	C+(i+2)*str1+c+j+2,
				K
			);
		for(size_t	i=maxa;	i<M;	i++)	for(size_t	j=0;	j<cache;	j++)
			C[i*str1+c+j]=vec_vec_dot_f(A+i*str,	B+(c+j)*str,	K);
	}
	#pragma omp parallel for
	for(size_t	i=0;	i<maxa;	i+=3)	for(size_t	j=maxc;	j<maxb;	j+=3)	
		mat_mul_3x3(
			A+i*str,	A+(i+1)*str,	A+(i+2)*str,
			B+j*str,	B+(j+1)*str,	B+(j+2)*str,
			C+i*str1+j,	C+i*str1+j+1,	C+i*str1+j+2,
			C+(i+1)*str1+j,	C+(i+1)*str1+j+1,	C+(i+1)*str1+j+2,
			C+(i+2)*str1+j,	C+(i+2)*str1+j+1,	C+(i+2)*str1+j+2,
			K
		);
	#pragma omp parallel for	
	for(size_t	j=maxc;	j<maxb;	j++)	for(size_t	i=maxa;	i<M;	i++)
		C[i*str1+j]=vec_vec_dot_f(A+i*str,	B+j*str,	K);
	#pragma omp parallel for
	for(size_t	i=0;	i<M;	i++)	for(size_t	j=maxb;	j<N;	j++)	
		C[i*str1+j]=vec_vec_dot_f(A+i*str,	B+j*str,	K);
}

/*
//	||A-B||
inline	float	vec_vec_dis_f(const	float	*A,	const	float	*B,	size_t	N){
	v4sf	s={0,0,0,0};	
	float	*r=(float*)&s;
	
	if(N<16)	switch(N){
		case	0:	return	0;
		case	1:	return	(A[0]-B[0])*(A[0]-B[0]);
		case	2:	return	(A[0]-B[0])*(A[0]-B[0])+(A[1]-B[1])*(A[1]-B[1]);
		case	3:	return	(A[0]-B[0])*(A[0]-B[0])+(A[1]-B[1])*(A[1]-B[1])+(A[2]-B[2])*(A[2]-B[2]);
		case	4:	s=__builtin_ia32_addps(s,	__builtin_ia32_mulps(__builtin_ia32_subps(*(v4sf*)(A),	*(v4sf*)(B)),	__builtin_ia32_subps(*(v4sf*)(A),	*(v4sf*)(B))));
					return	r[0]+r[1]+r[2]+r[3];
		case	5:	s=__builtin_ia32_addps(s,	__builtin_ia32_mulps(__builtin_ia32_subps(*(v4sf*)(A),	*(v4sf*)(B)),	__builtin_ia32_subps(*(v4sf*)(A),	*(v4sf*)(B))));
					return	r[0]+r[1]+r[2]+r[3]+(A[4]-B[4])*(A[4]-B[4]);
		case	6:	s=__builtin_ia32_addps(s,	__builtin_ia32_mulps(__builtin_ia32_subps(*(v4sf*)(A),	*(v4sf*)(B)),	__builtin_ia32_subps(*(v4sf*)(A),	*(v4sf*)(B))));
					return	r[0]+r[1]+r[2]+r[3]+(A[4]-B[4])*(A[4]-B[4])+(A[5]-B[5])*(A[5]-B[5]);
		case	7:	s=__builtin_ia32_addps(s,	__builtin_ia32_mulps(__builtin_ia32_subps(*(v4sf*)(A),	*(v4sf*)(B)),	__builtin_ia32_subps(*(v4sf*)(A),	*(v4sf*)(B))));
					return	r[0]+r[1]+r[2]+r[3]+(A[4]-B[4])*(A[4]-B[4])+(A[5]-B[5])*(A[5]-B[5])+(A[6]-B[6])*(A[6]-B[6]);
		case	8:	s=__builtin_ia32_addps(s,	__builtin_ia32_mulps(__builtin_ia32_subps(*(v4sf*)(A),	*(v4sf*)(B)),	__builtin_ia32_subps(*(v4sf*)(A),	*(v4sf*)(B))));		
					s=__builtin_ia32_addps(s,	__builtin_ia32_mulps(__builtin_ia32_subps(*(v4sf*)(A+4),	*(v4sf*)(B+4)),	__builtin_ia32_subps(*(v4sf*)(A+4),	*(v4sf*)(B+4))));
					return	r[0]+r[1]+r[2]+r[3];
		case	9:	s=__builtin_ia32_addps(s,	__builtin_ia32_mulps(__builtin_ia32_subps(*(v4sf*)(A),	*(v4sf*)(B)),	__builtin_ia32_subps(*(v4sf*)(A),	*(v4sf*)(B))));		
					s=__builtin_ia32_addps(s,	__builtin_ia32_mulps(__builtin_ia32_subps(*(v4sf*)(A+4),	*(v4sf*)(B+4)),	__builtin_ia32_subps(*(v4sf*)(A+4),	*(v4sf*)(B+4))));
					return	r[0]+r[1]+r[2]+r[3]+(A[8]-B[8])*(A[8]-B[8]);
		case	10:	s=__builtin_ia32_addps(s,	__builtin_ia32_mulps(__builtin_ia32_subps(*(v4sf*)(A),	*(v4sf*)(B)),	__builtin_ia32_subps(*(v4sf*)(A),	*(v4sf*)(B))));		
					s=__builtin_ia32_addps(s,	__builtin_ia32_mulps(__builtin_ia32_subps(*(v4sf*)(A+4),	*(v4sf*)(B+4)),	__builtin_ia32_subps(*(v4sf*)(A+4),	*(v4sf*)(B+4))));
					return	r[0]+r[1]+r[2]+r[3]+(A[8]-B[8])*(A[8]-B[8])+(A[9]-B[9])*(A[9]-B[9]);
		case	11:	s=__builtin_ia32_addps(s,	__builtin_ia32_mulps(__builtin_ia32_subps(*(v4sf*)(A),	*(v4sf*)(B)),	__builtin_ia32_subps(*(v4sf*)(A),	*(v4sf*)(B))));		
					s=__builtin_ia32_addps(s,	__builtin_ia32_mulps(__builtin_ia32_subps(*(v4sf*)(A+4),	*(v4sf*)(B+4)),	__builtin_ia32_subps(*(v4sf*)(A+4),	*(v4sf*)(B+4))));
					return	r[0]+r[1]+r[2]+r[3]+(A[8]-B[8])*(A[8]-B[8])+(A[9]-B[9])*(A[9]-B[9])+(A[10]-B[10])*(A[10]-B[10]);	
		case	12:	s=__builtin_ia32_addps(s,	__builtin_ia32_mulps(__builtin_ia32_subps(*(v4sf*)(A),	*(v4sf*)(B)),	__builtin_ia32_subps(*(v4sf*)(A),	*(v4sf*)(B))));		
					s=__builtin_ia32_addps(s,	__builtin_ia32_mulps(__builtin_ia32_subps(*(v4sf*)(A+4),	*(v4sf*)(B+4)),	__builtin_ia32_subps(*(v4sf*)(A+4),	*(v4sf*)(B+4))));
					s=__builtin_ia32_addps(s,	__builtin_ia32_mulps(__builtin_ia32_subps(*(v4sf*)(A+8),	*(v4sf*)(B+8)),	__builtin_ia32_subps(*(v4sf*)(A+8),	*(v4sf*)(B+8))));
					return	r[0]+r[1]+r[2]+r[3];
		case	13:	s=__builtin_ia32_addps(s,	__builtin_ia32_mulps(__builtin_ia32_subps(*(v4sf*)(A),	*(v4sf*)(B)),	__builtin_ia32_subps(*(v4sf*)(A),	*(v4sf*)(B))));		
					s=__builtin_ia32_addps(s,	__builtin_ia32_mulps(__builtin_ia32_subps(*(v4sf*)(A+4),	*(v4sf*)(B+4)),	__builtin_ia32_subps(*(v4sf*)(A+4),	*(v4sf*)(B+4))));
					s=__builtin_ia32_addps(s,	__builtin_ia32_mulps(__builtin_ia32_subps(*(v4sf*)(A+8),	*(v4sf*)(B+8)),	__builtin_ia32_subps(*(v4sf*)(A+8),	*(v4sf*)(B+8))));
					return	r[0]+r[1]+r[2]+r[3]+(A[12]-B[12])*(A[12]-B[12]);
		case	14:	s=__builtin_ia32_addps(s,	__builtin_ia32_mulps(__builtin_ia32_subps(*(v4sf*)(A),	*(v4sf*)(B)),	__builtin_ia32_subps(*(v4sf*)(A),	*(v4sf*)(B))));		
					s=__builtin_ia32_addps(s,	__builtin_ia32_mulps(__builtin_ia32_subps(*(v4sf*)(A+4),	*(v4sf*)(B+4)),	__builtin_ia32_subps(*(v4sf*)(A+4),	*(v4sf*)(B+4))));
					s=__builtin_ia32_addps(s,	__builtin_ia32_mulps(__builtin_ia32_subps(*(v4sf*)(A+8),	*(v4sf*)(B+8)),	__builtin_ia32_subps(*(v4sf*)(A+8),	*(v4sf*)(B+8))));
					return	r[0]+r[1]+r[2]+r[3]+(A[12]-B[12])*(A[12]-B[12])+(A[13]-B[13])*(A[13]-B[13]);
		case	15:	s=__builtin_ia32_addps(s,	__builtin_ia32_mulps(__builtin_ia32_subps(*(v4sf*)(A),	*(v4sf*)(B)),	__builtin_ia32_subps(*(v4sf*)(A),	*(v4sf*)(B))));		
					s=__builtin_ia32_addps(s,	__builtin_ia32_mulps(__builtin_ia32_subps(*(v4sf*)(A+4),	*(v4sf*)(B+4)),	__builtin_ia32_subps(*(v4sf*)(A+4),	*(v4sf*)(B+4))));
					s=__builtin_ia32_addps(s,	__builtin_ia32_mulps(__builtin_ia32_subps(*(v4sf*)(A+8),	*(v4sf*)(B+8)),	__builtin_ia32_subps(*(v4sf*)(A+8),	*(v4sf*)(B+8))));
					return	r[0]+r[1]+r[2]+r[3]+(A[12]-B[12])*(A[12]-B[12])+(A[13]-B[13])*(A[13]-B[13])+(A[14]-B[14])*(A[14]-B[14]);
	}
	size_t	max=(N>>4)<<4;
	for(size_t	i=0;	i<max;	i+=16){
		s=__builtin_ia32_addps(s,	__builtin_ia32_mulps(__builtin_ia32_subps(*(v4sf*)(A+i),	*(v4sf*)(B+i)),	__builtin_ia32_subps(*(v4sf*)(A+i),	*(v4sf*)(B+i))));
		s=__builtin_ia32_addps(s,	__builtin_ia32_mulps(__builtin_ia32_subps(*(v4sf*)(A+i+4),	*(v4sf*)(B+i+4)),	__builtin_ia32_subps(*(v4sf*)(A+i+4),	*(v4sf*)(B+i+4))));
		s=__builtin_ia32_addps(s,	__builtin_ia32_mulps(__builtin_ia32_subps(*(v4sf*)(A+i+8),	*(v4sf*)(B+i+8)),	__builtin_ia32_subps(*(v4sf*)(A+i+8),	*(v4sf*)(B+i+8))));
		s=__builtin_ia32_addps(s,	__builtin_ia32_mulps(__builtin_ia32_subps(*(v4sf*)(A+i+12),	*(v4sf*)(B+i+12)),	__builtin_ia32_subps(*(v4sf*)(A+i+12),	*(v4sf*)(B+i+12))));
	}
	A+=max;	B+=max;
	switch(N&15){
		case	0:	return	r[0]+r[1]+r[2]+r[3];
		case	1:	return	r[0]+r[1]+r[2]+r[3]+(A[0]-B[0])*(A[0]-B[0]);
		case	2:	return	r[0]+r[1]+r[2]+r[3]+(A[0]-B[0])*(A[0]-B[0])+(A[1]-B[1])*(A[1]-B[1]);
		case	3:	return	r[0]+r[1]+r[2]+r[3]+(A[0]-B[0])*(A[0]-B[0])+(A[1]-B[1])*(A[1]-B[1])+(A[2]-B[2])*(A[2]-B[2]);
		case	4:	s=__builtin_ia32_addps(s,	__builtin_ia32_mulps(__builtin_ia32_subps(*(v4sf*)(A),	*(v4sf*)(B)),	__builtin_ia32_subps(*(v4sf*)(A),	*(v4sf*)(B))));
					return	r[0]+r[1]+r[2]+r[3];
		case	5:	s=__builtin_ia32_addps(s,	__builtin_ia32_mulps(__builtin_ia32_subps(*(v4sf*)(A),	*(v4sf*)(B)),	__builtin_ia32_subps(*(v4sf*)(A),	*(v4sf*)(B))));
					return	r[0]+r[1]+r[2]+r[3]+(A[4]-B[4])*(A[4]-B[4]);
		case	6:	s=__builtin_ia32_addps(s,	__builtin_ia32_mulps(__builtin_ia32_subps(*(v4sf*)(A),	*(v4sf*)(B)),	__builtin_ia32_subps(*(v4sf*)(A),	*(v4sf*)(B))));
					return	r[0]+r[1]+r[2]+r[3]+(A[4]-B[4])*(A[4]-B[4])+(A[5]-B[5])*(A[5]-B[5]);
		case	7:	s=__builtin_ia32_addps(s,	__builtin_ia32_mulps(__builtin_ia32_subps(*(v4sf*)(A),	*(v4sf*)(B)),	__builtin_ia32_subps(*(v4sf*)(A),	*(v4sf*)(B))));
					return	r[0]+r[1]+r[2]+r[3]+(A[4]-B[4])*(A[4]-B[4])+(A[5]-B[5])*(A[5]-B[5])+(A[6]-B[6])*(A[6]-B[6]);
		case	8:	s=__builtin_ia32_addps(s,	__builtin_ia32_mulps(__builtin_ia32_subps(*(v4sf*)(A),	*(v4sf*)(B)),	__builtin_ia32_subps(*(v4sf*)(A),	*(v4sf*)(B))));		
					s=__builtin_ia32_addps(s,	__builtin_ia32_mulps(__builtin_ia32_subps(*(v4sf*)(A+4),	*(v4sf*)(B+4)),	__builtin_ia32_subps(*(v4sf*)(A+4),	*(v4sf*)(B+4))));
					return	r[0]+r[1]+r[2]+r[3];
		case	9:	s=__builtin_ia32_addps(s,	__builtin_ia32_mulps(__builtin_ia32_subps(*(v4sf*)(A),	*(v4sf*)(B)),	__builtin_ia32_subps(*(v4sf*)(A),	*(v4sf*)(B))));		
					s=__builtin_ia32_addps(s,	__builtin_ia32_mulps(__builtin_ia32_subps(*(v4sf*)(A+4),	*(v4sf*)(B+4)),	__builtin_ia32_subps(*(v4sf*)(A+4),	*(v4sf*)(B+4))));
					return	r[0]+r[1]+r[2]+r[3]+(A[8]-B[8])*(A[8]-B[8]);
		case	10:	s=__builtin_ia32_addps(s,	__builtin_ia32_mulps(__builtin_ia32_subps(*(v4sf*)(A),	*(v4sf*)(B)),	__builtin_ia32_subps(*(v4sf*)(A),	*(v4sf*)(B))));		
					s=__builtin_ia32_addps(s,	__builtin_ia32_mulps(__builtin_ia32_subps(*(v4sf*)(A+4),	*(v4sf*)(B+4)),	__builtin_ia32_subps(*(v4sf*)(A+4),	*(v4sf*)(B+4))));
					return	r[0]+r[1]+r[2]+r[3]+(A[8]-B[8])*(A[8]-B[8])+(A[9]-B[9])*(A[9]-B[9]);
		case	11:	s=__builtin_ia32_addps(s,	__builtin_ia32_mulps(__builtin_ia32_subps(*(v4sf*)(A),	*(v4sf*)(B)),	__builtin_ia32_subps(*(v4sf*)(A),	*(v4sf*)(B))));		
					s=__builtin_ia32_addps(s,	__builtin_ia32_mulps(__builtin_ia32_subps(*(v4sf*)(A+4),	*(v4sf*)(B+4)),	__builtin_ia32_subps(*(v4sf*)(A+4),	*(v4sf*)(B+4))));
					return	r[0]+r[1]+r[2]+r[3]+(A[8]-B[8])*(A[8]-B[8])+(A[9]-B[9])*(A[9]-B[9])+(A[10]-B[10])*(A[10]-B[10]);	
		case	12:	s=__builtin_ia32_addps(s,	__builtin_ia32_mulps(__builtin_ia32_subps(*(v4sf*)(A),	*(v4sf*)(B)),	__builtin_ia32_subps(*(v4sf*)(A),	*(v4sf*)(B))));		
					s=__builtin_ia32_addps(s,	__builtin_ia32_mulps(__builtin_ia32_subps(*(v4sf*)(A+4),	*(v4sf*)(B+4)),	__builtin_ia32_subps(*(v4sf*)(A+4),	*(v4sf*)(B+4))));
					s=__builtin_ia32_addps(s,	__builtin_ia32_mulps(__builtin_ia32_subps(*(v4sf*)(A+8),	*(v4sf*)(B+8)),	__builtin_ia32_subps(*(v4sf*)(A+8),	*(v4sf*)(B+8))));
					return	r[0]+r[1]+r[2]+r[3];
		case	13:	s=__builtin_ia32_addps(s,	__builtin_ia32_mulps(__builtin_ia32_subps(*(v4sf*)(A),	*(v4sf*)(B)),	__builtin_ia32_subps(*(v4sf*)(A),	*(v4sf*)(B))));		
					s=__builtin_ia32_addps(s,	__builtin_ia32_mulps(__builtin_ia32_subps(*(v4sf*)(A+4),	*(v4sf*)(B+4)),	__builtin_ia32_subps(*(v4sf*)(A+4),	*(v4sf*)(B+4))));
					s=__builtin_ia32_addps(s,	__builtin_ia32_mulps(__builtin_ia32_subps(*(v4sf*)(A+8),	*(v4sf*)(B+8)),	__builtin_ia32_subps(*(v4sf*)(A+8),	*(v4sf*)(B+8))));
					return	r[0]+r[1]+r[2]+r[3]+(A[12]-B[12])*(A[12]-B[12]);
		case	14:	s=__builtin_ia32_addps(s,	__builtin_ia32_mulps(__builtin_ia32_subps(*(v4sf*)(A),	*(v4sf*)(B)),	__builtin_ia32_subps(*(v4sf*)(A),	*(v4sf*)(B))));		
					s=__builtin_ia32_addps(s,	__builtin_ia32_mulps(__builtin_ia32_subps(*(v4sf*)(A+4),	*(v4sf*)(B+4)),	__builtin_ia32_subps(*(v4sf*)(A+4),	*(v4sf*)(B+4))));
					s=__builtin_ia32_addps(s,	__builtin_ia32_mulps(__builtin_ia32_subps(*(v4sf*)(A+8),	*(v4sf*)(B+8)),	__builtin_ia32_subps(*(v4sf*)(A+8),	*(v4sf*)(B+8))));
					return	r[0]+r[1]+r[2]+r[3]+(A[12]-B[12])*(A[12]-B[12])+(A[13]-B[13])*(A[13]-B[13]);
		case	15:	s=__builtin_ia32_addps(s,	__builtin_ia32_mulps(__builtin_ia32_subps(*(v4sf*)(A),	*(v4sf*)(B)),	__builtin_ia32_subps(*(v4sf*)(A),	*(v4sf*)(B))));		
					s=__builtin_ia32_addps(s,	__builtin_ia32_mulps(__builtin_ia32_subps(*(v4sf*)(A+4),	*(v4sf*)(B+4)),	__builtin_ia32_subps(*(v4sf*)(A+4),	*(v4sf*)(B+4))));
					s=__builtin_ia32_addps(s,	__builtin_ia32_mulps(__builtin_ia32_subps(*(v4sf*)(A+8),	*(v4sf*)(B+8)),	__builtin_ia32_subps(*(v4sf*)(A+8),	*(v4sf*)(B+8))));
					return	r[0]+r[1]+r[2]+r[3]+(A[12]-B[12])*(A[12]-B[12])+(A[13]-B[13])*(A[13]-B[13])+(A[14]-B[14])*(A[14]-B[14]);
	}	
	return	0;
}
*/

/*
//	A+=B*C
inline	void	vec_vec_fma(float	*A,	const	float	*B,	const	float	*C,	size_t	N){
	size_t	max=(N>>4)<<4;
	for(size_t	i=0;	i<max;	i+=16){
		*(v4sf*)(A+i)=__builtin_ia32_addps(*(v4sf*)(A+i),	__builtin_ia32_mulps(*(v4sf*)(B+i),	*(v4sf*)(C+i)));
		*(v4sf*)(A+i+4)=__builtin_ia32_addps(*(v4sf*)(A+i+4),	__builtin_ia32_mulps(*(v4sf*)(B+i+4),	*(v4sf*)(C+i+4)));
		*(v4sf*)(A+i+8)=__builtin_ia32_addps(*(v4sf*)(A+i+8),	__builtin_ia32_mulps(*(v4sf*)(B+i+8),	*(v4sf*)(C+i+8)));
		*(v4sf*)(A+i+12)=__builtin_ia32_addps(*(v4sf*)(A+i+12),	__builtin_ia32_mulps(*(v4sf*)(B+i+12),	*(v4sf*)(C+i+12)));
	}

	float	*p=A+max;	const	float	*q=B+max,	*r=C+max;
	switch(N&15){
		case	1:	p[0]+=q[0]*r[0];	break;
		case	2:	p[0]+=q[0]*r[0];	p[1]+=q[1]*r[1];	break;
		case	3:	p[0]+=q[0]*r[0];	p[1]+=q[1]*r[1];	p[2]+=q[2]*r[2];	break;
		case	4:	*(v4sf*)p=__builtin_ia32_addps(*(v4sf*)p,	__builtin_ia32_mulps(*(v4sf*)q,	*(v4sf*)r));	
					break;
		case	5:	*(v4sf*)p=__builtin_ia32_addps(*(v4sf*)p,	__builtin_ia32_mulps(*(v4sf*)q,	*(v4sf*)r));	
					p[4]+=q[4]*r[4];	break;
		case	6:	*(v4sf*)p=__builtin_ia32_addps(*(v4sf*)p,	__builtin_ia32_mulps(*(v4sf*)q,	*(v4sf*)r));	
					p[4]+=q[4]*r[4];	p[5]+=q[5]*r[5];	break;
		case	7:	*(v4sf*)p=__builtin_ia32_addps(*(v4sf*)p,	__builtin_ia32_mulps(*(v4sf*)q,	*(v4sf*)r));	
					p[4]+=q[4]*r[4];	p[5]+=q[5]*r[5];	p[6]+=q[6]*r[6];	break;
		case	8:	*(v4sf*)p=__builtin_ia32_addps(*(v4sf*)p,	__builtin_ia32_mulps(*(v4sf*)q,	*(v4sf*)r));
					*(v4sf*)(p+4)=__builtin_ia32_addps(*(v4sf*)(p+4),	__builtin_ia32_mulps(*(v4sf*)(q+4),	*(v4sf*)(r+4)));
					break;
		case	9:	*(v4sf*)p=__builtin_ia32_addps(*(v4sf*)p,	__builtin_ia32_mulps(*(v4sf*)q,	*(v4sf*)r));
					*(v4sf*)(p+4)=__builtin_ia32_addps(*(v4sf*)(p+4),	__builtin_ia32_mulps(*(v4sf*)(q+4),	*(v4sf*)(r+4)));
					p[8]+=q[8]*r[8];	break;
		case	10:	*(v4sf*)p=__builtin_ia32_addps(*(v4sf*)p,	__builtin_ia32_mulps(*(v4sf*)q,	*(v4sf*)r));
					*(v4sf*)(p+4)=__builtin_ia32_addps(*(v4sf*)(p+4),	__builtin_ia32_mulps(*(v4sf*)(q+4),	*(v4sf*)(r+4)));
					p[8]+=q[8]*r[8];	p[9]+=q[9]*r[9];	break;
		case	11:	*(v4sf*)p=__builtin_ia32_addps(*(v4sf*)p,	__builtin_ia32_mulps(*(v4sf*)q,	*(v4sf*)r));
					*(v4sf*)(p+4)=__builtin_ia32_addps(*(v4sf*)(p+4),	__builtin_ia32_mulps(*(v4sf*)(q+4),	*(v4sf*)(r+4)));
					p[8]+=q[8]*r[8];	p[9]+=q[9]*r[9];	p[10]+=q[10]*r[10];break;
		case	12:	*(v4sf*)p=__builtin_ia32_addps(*(v4sf*)p,	__builtin_ia32_mulps(*(v4sf*)q,	*(v4sf*)r));
					*(v4sf*)(p+4)=__builtin_ia32_addps(*(v4sf*)(p+4),	__builtin_ia32_mulps(*(v4sf*)(q+4),	*(v4sf*)(r+4)));
					*(v4sf*)(p+8)=__builtin_ia32_addps(*(v4sf*)(p+8),	__builtin_ia32_mulps(*(v4sf*)(q+8),	*(v4sf*)(r+8)));
					break;
		case	13:	*(v4sf*)p=__builtin_ia32_addps(*(v4sf*)p,	__builtin_ia32_mulps(*(v4sf*)q,	*(v4sf*)r));
					*(v4sf*)(p+4)=__builtin_ia32_addps(*(v4sf*)(p+4),	__builtin_ia32_mulps(*(v4sf*)(q+4),	*(v4sf*)(r+4)));
					*(v4sf*)(p+8)=__builtin_ia32_addps(*(v4sf*)(p+8),	__builtin_ia32_mulps(*(v4sf*)(q+8),	*(v4sf*)(r+8)));		
					p[12]+=q[12]*r[12];	break;
		case	14:	*(v4sf*)p=__builtin_ia32_addps(*(v4sf*)p,	__builtin_ia32_mulps(*(v4sf*)q,	*(v4sf*)r));
					*(v4sf*)(p+4)=__builtin_ia32_addps(*(v4sf*)(p+4),	__builtin_ia32_mulps(*(v4sf*)(q+4),	*(v4sf*)(r+4)));
					*(v4sf*)(p+8)=__builtin_ia32_addps(*(v4sf*)(p+8),	__builtin_ia32_mulps(*(v4sf*)(q+8),	*(v4sf*)(r+8)));		
					p[12]+=q[12]*r[12];	p[13]+=q[13]*r[13];	break;
		case	15:	*(v4sf*)p=__builtin_ia32_addps(*(v4sf*)p,	__builtin_ia32_mulps(*(v4sf*)q,	*(v4sf*)r));
					*(v4sf*)(p+4)=__builtin_ia32_addps(*(v4sf*)(p+4),	__builtin_ia32_mulps(*(v4sf*)(q+4),	*(v4sf*)(r+4)));
					*(v4sf*)(p+8)=__builtin_ia32_addps(*(v4sf*)(p+8),	__builtin_ia32_mulps(*(v4sf*)(q+8),	*(v4sf*)(r+8)));		
					p[12]+=q[12]*r[12];	p[13]+=q[13]*r[13];	p[14]+=q[14]*r[14];	break;
	}
	
}
*/

#endif
