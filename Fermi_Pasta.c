#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "omp.h"


#define PI 3.14159265

int N=64;
double Beta=1.;
int Nt;
double dt=0.005;
double *x,*xtemp;
double *v,*vtemp;
double *F;
double Q1,Q2,Q3;
double Qp1,Qp2,Qp3;
double E1,E2,E3;

double x_n(int n);
double Q(int k,double *x);
double sum_array(double *a, int num_elements);
double wk2(int k);
void lfs(double x,double v, double F);
void leapfrogStep(double *x,double *v,double *F,int size);

int main(int argc, char *argv[]){

  FILE * energia;
  energia=fopen("energia.txt","w");

  int word_size=atoi(argv[1]);
  omp_set_num_threads(word_size);
  Nt=(int)5*pow(N,2.2)/dt;
  x=malloc(N*sizeof(double));
  xtemp=malloc(N*sizeof(double));
  v=malloc(N*sizeof(double));
  vtemp=malloc(N*sizeof(double));
  int sample_size=N/word_size;//tamanio de lo que evaluara cada procesador
  F=malloc((N)*sizeof(double));

  //condicion inicial de x
  int i;
  for (i=0;i<N;i++){
    x[i]=x_n(i);
    if(i==N-1){//ya que PI es una aproximacion, hay que garantizar que x_{N-1} sea 0
      x[i]=0.;
    }
  //de F
  F[0]=0.;
  F[N-1]=0.;
  //de v
  v[0]=0.;
  v[N-1]=0.;
  }
  float w1=wk2(1);
  float w2=wk2(2);
  float w3=wk2(3);

  for (i=0;i<Nt;i++){
    int j;

    //printf("%d",i);
    #pragma omp parallel for private (j),shared(x,v,F)
    for(j=1;j<N-1;j++){
       double x2=x[j+1]-x[j];
       double x3=x[j]-x[j-1];
       F[j]=(x[j+1]-2*x[j]+x[j-1])+Beta*pow(x2,2)-pow(x3,3);
       if(i==0){
          v[j]=dt*F[j]/2.0;
       }
       xtemp[j]=x[j];
       vtemp[j]=v[j];
       //lfs(x[j],v[j],F[j]);       
    }
    if(i<1000){

    Q1=Q(1,x);
    Q2=Q(2,x);
    Q3=Q(3,x);
    Qp1=Q(1,v);
    Qp2=Q(2,v);
    Qp3=Q(3,v);
    E1=0.5*(pow(Qp1,2)+w1*pow(Q1,2));
    E2=0.5*(pow(Qp2,2)+w1*pow(Q2,2));
    E3=0.5*(pow(Qp3,2)+w1*pow(Q3,2));
    //printf("%lf,%.lf,%.lf,%lf\n",E1,E2,E3,i*dt);
    fprintf(energia,"%lf,%.lf,%.lf,%lf\n",E1,E2,E3,i*dt);
    }
    //leapfrogStep(x,v,F,N);
    #pragma omp parallel for private (j),shared(x,v,F,xtemp,vtemp)
    for (j=1;j<N-1;j++){
       v[j]=vtemp[j]+F[j]*dt;
       x[j]=xtemp[j]+v[j]*dt;
    }
  }
}


double wk2(int k){
   double kf=(double)k;
   return 4.0*pow(sin(kf*PI/(2*N+2)),2);
}
double x_n(int n){
  return sin(PI*n/(N-1));
}
double Q(int k,double *x){

  double Qk=sqrt(2.0/(N+1))*sum_array(x,N);
  int i;
  double sines=0.0;
  for (i=0;i<N;i++){
     sines +=sin(k*PI*i/N+1);
  }
  return Qk*sines;
}

double sum_array(double *a, int num_elements)
{
   int i;
   double sum=0;
   for (i=0; i<num_elements; i++)
   {
	 sum = sum + a[i];
   }
   return(sum);
}
void lfs(double x,double v, double F){
   double xi;
   double vi;
   xi=x;
   vi=v;
   xi+=vi*dt;
   vi+=F*dt;
   x=xi;
   vi=v;
}
void leapfrogStep(double *x,double *v,double *F,int size){
printf("entro\n");
   double *x_in;
   double *v_in;
   *x_in=*x;
   *v_in=*v;
   int k;
   #pragma omp parallel for private (k),shared(x_in,v_in,F)
   for (k=0;k<size-1;k++){
      printf("%f\n",v_in[k]);
      x_in[k]+=v_in[k]*dt;/*drift*/
      v_in[k]+=F[k]*dt;/*kick*/

   }
   
   *x=*x_in;
   *v=*v_in;
}

