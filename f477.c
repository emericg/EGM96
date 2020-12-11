/*translation of the well-known f477.f to C
 actually this source is folklore
d.ineiev<ineiev@yahoo.co.uk> wrote down and put under zlib/libpng license
(see README.txt)*/
/* this program is designed for the calculation of a geoid undulation
 at a point whose latitude and longitude is specified. the program
 is designed to use the potential coefficient model egm96 and a
 set of spherical harmonic coefficients of a correction term.
 the correction term is composed of several different components
 the primary one being the conversion of a height anomaly to a geoid
 undulation. the principles of this procedure were initially
 described in the paper: use of potential coefficient models for geoid
 undulation determination using a spherical harmonic representation
 of the height anomaly/geoid undulation difference by R.H. Rapp,
 Journal of Geodesy, 1996.
 this program is designed to be used with the constants of egm96
 and those of the wgs84(g873) system. the undulation will refer to
 the wgs84 ellipsoid.
 specific details on the undulation computation will be found in the
 joint project report describing the development of egm96.
 this program is a modification of the program described in the
 following report:
 a fortran program for the computation of gravimetric quantities from
 high degree spherical harmonic expansions, Richard H. Rapp,
 report 334, Department of Geodetic Science and Surveying, the Ohio
 State University, Columbus, 1982
 this program was put in this form in Dec 1996.
 rhrapp.f477.nonly
 


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  the input files consist of:

                correction coefficient set ("CORRCOEF") => unit = 1
                    potential coefficient set ("EGM96") => unit = 12
                points at which to compute ("INPUT.dat") => unit = 14

  the output file is:

                    computed geoid heights ("OUTF477") => unit = 20

    file assignment revisions at NIMA, December 1996.

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

   dimensions of p,q,hc,hs must be at least ((maxn+1)*(maxn+2))/2,
   dimensions of sinml,cosml must be at least maxn,
        where maxn is maximum order of computation
 the current dimensions are set for a maximum degree of 360*/
#include<stdio.h>
#include<math.h>
static FILE*f_1,*f_14,*f_12,*f_20;
#define l_value	(65341)
#define _361	(361)
double hundu(unsigned nmax,double p[l_value+1],
 double hc[l_value+1],double hs[l_value+1],
 double sinml[_361+1],double cosml[_361+1],double gr,double re,
 double cc[l_value+1],double cs[l_value+1])
{/*constants for wgs84(g873);gm in units of m**3/s**2*/
 const double gm=.3986004418e15,ae=6378137.;
 double arn,ar,ac,a,b,sum,sumc,sum2,tempc,temp;int k,n,m;
 ar=ae/re;arn=ar;ac=a=b=0;k=3;
 for(n=2;n<=nmax;n++)
 {arn*=ar;k++;sum=p[k]*hc[k];sumc=p[k]*cc[k];sum2=0;
  for(m=1;m<=n;m++)
  {k++;tempc=cc[k]*cosml[m]+cs[k]*sinml[m];
   temp=hc[k]*cosml[m]+hs[k]*sinml[m];sumc+=p[k]*tempc;sum+=p[k]*temp;
  }ac+=sumc;a+=sum*arn;
 }ac+=cc[1]+p[2]*cc[2]+p[3]*(cc[3]*cosml[1]+cs[3]*sinml[1]);
/*add haco=ac/100 to convert height anomaly on the ellipsoid to the undulation
add -0.53m to make undulation refer to the wgs84 ellipsoid.*/
 return a*gm/(gr*re)+ac/100-.53;
}
void dscml(double rlon,unsigned nmax,double sinml[_361+1],double cosml[_361+1])
{double a,b;int m;a=sin(rlon);b=cos(rlon);
 sinml[1]=a;cosml[1]=b;sinml[2]=2*b*a;cosml[2]=2*b*b-1;
 for(m=3;m<=nmax;m++)
 {sinml[m]=2*b*sinml[m-1]-sinml[m-2];cosml[m]=2*b*cosml[m-1]-cosml[m-2];}
}
void dhcsin(unsigned nmax,double hc[l_value+1],double hs[l_value+1])
{int n,m;double j2,j4,j6,j8,j10,c,s,ec,es;
/*the even degree zonal coefficients given below were computed for the
 wgs84(g873) system of constants and are identical to those values
 used in the NIMA gridding procedure. computed using subroutine
 grs written by N.K. PAVLIS*/
 j2=0.108262982131e-2;j4=-.237091120053e-05;j6=0.608346498882e-8;
 j8=-0.142681087920e-10;j10=0.121439275882e-13;m=((nmax+1)*(nmax+2))/2;
 for(n=1;n<=m;n++)hc[n]=hs[n]=0;
 while(6==fscanf(f_12,"%i %i %lf %lf %lf %lf",&n,&m,&c,&s,&ec,&es))
 {if(n>nmax)continue;n=(n*(n+1))/2+m+1;hc[n]=c;hs[n]=s;}
 hc[4]+=j2/sqrt(5);hc[11]+=j4/3;hc[22]+=j6/sqrt(13);hc[37]+=j8/sqrt(17);
 hc[56]+=j10/sqrt(21);
}
void legfdn(unsigned m,double theta,double rleg[_361+1],unsigned nmx)
/*this subroutine computes  all normalized legendre function
in "rleg". order is always
m, and colatitude is always theta  (radians). maximum deg
is  nmx. all calculations in double precision.
ir  must be set to zero before the first call to this sub.
the dimensions of arrays  rleg must be at least equal to  nmx+1.
Original programmer :Oscar L. Colombo, Dept. of Geodetic Science
the Ohio State University, August 1980
ineiev: I removed the derivatives, for they are never computed here*/
{static double drts[1301],dirt[1301],cothet,sithet,rlnn[_361+1];
 static int ir;int nmx1=nmx+1,nmx2p=2*nmx+1,m1=m+1,m2=m+2,m3=m+3,n,n1,n2;
 if(!ir){ir=1;for(n=1;n<=nmx2p;n++){drts[n]=sqrt(n);dirt[n]=1/drts[n];}}
 cothet=cos(theta);sithet=sin(theta);
 /*compute the legendre functions*/
 rlnn[1]=1;rlnn[2]=sithet*drts[3];
 for(n1=3;n1<=m1;n1++)
 {n=n1-1;n2=2*n;rlnn[n1]=drts[n2+1]*dirt[n2]*sithet*rlnn[n];}
 switch(m)
 {case 1:rleg[2]=rlnn[2];rleg[3]=drts[5]*cothet*rleg[2];break;
  case 0:rleg[1]=1;rleg[2]=cothet*drts[3];break;
 }rleg[m1]=rlnn[m1];
 if(m2<=nmx1)
 {rleg[m2]=drts[m1*2+1]*cothet*rleg[m1];
  if(m3<=nmx1)for(n1=m3;n1<=nmx1;n1++)
  {n=n1-1;if((!m&&n<2)||(m==1&&n<3))continue;
   n2=2*n;rleg[n1]=drts[n2+1]*dirt[n+m]*dirt[n-m]*
    (drts[n2-1]*cothet*rleg[n1-1]-drts[n+m-1]*drts[n-m-1]*dirt[n2-3]*rleg[n1-2]);
  }
 }
}
void radgra(double lat,double lon,double*rlat,double*gr,double*re)
/*this subroutine computes geocentric distance to the point,
the geocentric latitude,and
an approximate value of normal gravity at the point based
the constants of the wgs84(g873) system are used*/
{const double a=6378137.,e2=.00669437999013,geqt=9.7803253359,k=.00193185265246;
 double n,t1=sin(lat)*sin(lat),t2,x,y,z;
 n=a/sqrt(1-e2*t1);t2=n*cos(lat);x=t2*cos(lon);y=t2*sin(lon);
 z=(n*(1-e2))*sin(lat);
 *re=sqrt(x*x+y*y+z*z);/*compute the geocentric radius*/
 *rlat=atan(z/sqrt(x*x+y*y));/*compute the geocentric latitude*/
 *gr=geqt*(1+k*t1)/sqrt(1-e2*t1);/*compute normal gravity:units are m/sec**2*/
}
static double cc[l_value+1],cs[l_value+1],hc[l_value+1],hs[l_value+1],
 p[l_value+1],sinml[_361+1],cosml[_361+1],rleg[_361+1];
static int nmax;
double undulation(double lat,double lon,int nmax,int k)
{double rlat,gr,re;int i,j,m;
 radgra(lat,lon,&rlat,&gr,&re);rlat=M_PI/2-rlat;
 for(j=1;j<=k;j++)
 {m=j-1;legfdn(m,rlat,rleg,nmax);for(i=j;i<=k;i++)p[(i-1)*i/2+m+1]=rleg[i];}
 dscml(lon,nmax,sinml,cosml);return hundu(nmax,p,hc,hs,sinml,cosml,gr,re,cc,cs);
}
void init_arrays(void)
{int ig,i,n,m;double t1,t2;
 f_1=fopen("CORCOEF","rb");/*correction coefficient file: 
  modified with 'sed -e"s/D/e/g"' to be read with fscanf*/
 f_12=fopen("EGM96","rb");/*potential coefficient file*/
 nmax=360;for(i=1;i<=l_value;i++)cc[i]=cs[i]=0;
 while(4==fscanf(f_1,"%i %i %lg %lg",&n,&m,&t1,&t2))
 {ig=(n*(n+1))/2+m+1;cc[ig]=t1;cs[ig]=t2;}
/*the correction coefficients are now read in*/
/*the potential coefficients are now read in and the reference
 even degree zonal harmonic coefficients removed to degree 6*/
 dhcsin(nmax,hc,hs);fclose(f_1);fclose(f_12);
}
int main(void)
{const double rad=180/M_PI;double flat,flon,u;init_arrays();
 f_14=fopen("INPUT.DAT","rb");/*input data file*/
 f_20=fopen("OUTF477.DAT","wb");/*output file*/
/*read geodetic latitude,longitude at point undulation is wanted*/
 while(2==fscanf(f_14,"%lg %lg",&flat,&flon))
 {/*compute the geocentric latitude,geocentric radius,normal gravity*/
  u=undulation(flat/rad,flon/rad,nmax,nmax+1);
  /*u is the geoid undulation from the egm96 potential coefficient model
   including the height anomaly to geoid undulation correction term
   and a correction term to have the undulations refer to the
   wgs84 ellipsoid. the geoid undulation unit is meters.*/
  fprintf(f_20,"%14.7f %14.7f %10.7f\n",flat,flon,u);
 }fclose(f_14);fclose(f_20);return 0;
}
