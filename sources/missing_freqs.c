#include "missing_freqs.h"

/*************************************************
 ngstat 
 
 source code (c) Emanuele Raineri 2012
 emanuele.raineri@gmail.com
 
 implements the formulae contained in
 Luca Ferretti et al: 
 Neutrality tests from high-throughput resequencing...
 and
 Luca Ferretti, Emanuele Raineri, Sebastian RO:
 
 Neutrality tests for sequences with missing data
 
 (the following license has been ripped from
 http://www.mccaughan.org.uk/software/disarm.c-0.11)
 
 
 This file may be distributed and used freely provided:
 * 1. You do not distribute any version that lacks this
 *    copyright notice (exactly as it appears here, extending
 *    from the start to the end of the C-language comment
 *    containing these words)); and,
 * 2. If you distribute any modified version, its source
 *    contains a clear description of the ways in which
 *    it differs from the original version, and a clear
 *    indication that the changes are not mine.
 * There is no restriction on your permission to use and
 * distribute object code or executable code derived from
 * this.
 
 **************************************************/

/*
 
 P_c(j|n_1\ldots n_N)=
 \sum_{i_1=1}^{2}\ldots\sum_{i_N=1}^{2}I\left(j=\sum_{p=1}^N i_p\right)\prod_{q=1}^N \left(2^{-n_q+1}+(i_q-1)(1-2^{-n_q+2})\right)
 
 */

/* I create and store all the configurations which
 add up to j in a bidimensional array.
 I could create and store all the configurations, for 
 different values of j, and store them in a 3d array 
 */

int min3(int a, int b, int c){
	return min(min(a,b),c);
}
int max3(int a, int b, int c){
	return max(max(a,b),c);
}

int myround(float r){
	// XXX check , why is this necessary?
    if (r-floorf(r)<0.5) return (int)floorf(r);
	return (int)ceilf(r);
}

float lnfact(int n){
	/*static float table[1000];*/
	static int max_n=0;
	static float *table;
	if(max_n == 0) {
		max_n = 1000;
		table = (float *)calloc(max_n,sizeof(float));
	}
	if(n+1 > max_n) {
		max_n = n+1;
		table = (float *)realloc(table,max_n*sizeof(float));
	}
	return (table[n]>0)?table[n]:(table[n]=lgamma(n+1)); 
}

float lnbinomial(int n, int k){
	/*static float table[100*100];*/ /*PROBLEM: NO MORE THAN 100 INDIVIDUALS!!!!*/ 
	static int max_n=0;
	static float **table;
	int x;
	if(max_n == 0) {
		max_n = 100;
		table = (float **)calloc(max_n,sizeof(float *));
		for(x=0;x<max_n;x++) {
			table[x] = (float *)calloc(max_n,sizeof(float));
		}
	}
	if(n+1 > max_n) {
		table = (float **)realloc(table,(n+1)*sizeof(float *));
		for(x=max_n;x<n+1;x++) {
			table[x] = (float *)calloc(n+1,sizeof(float));
		}
		for(x=0;x<max_n;x++) {
			table[x] = (float *)realloc(table[x],(n+1)*sizeof(float));
		}
		max_n = n+1;
	}
	if (n<0 || k<0)  {return -INF;}
	if (k>n) {return -INF;}
	//return lnfact(n)-lnfact(k)-lnfact(n-k);
	/*return (table[n*100+k]>0)?table[n*100+k]:(table[n*100+k]=lnfact(n)-lnfact(k)-lnfact(n-k)); */
	return (table[n][k]>0)?table[n][k]:(table[n][k]=lnfact(n)-lnfact(k)-lnfact(n-k)); 
}
float lnmultinomial(int n, int k1, int k2){
	if (n-k1-k2<0  ) {return -INF;}
	if (n <0 || k1<0 || k2<0) {return -INF;}
	return lnfact(n)-lnfact(k1)-lnfact(k2)-lnfact(n-k1-k2);
}

float a(int n){
	if (0==n) {fprintf(stderr,"a:invalid n\n");exit(1);}
	int i;
	float s=0.0;
	for (i=1;i<n;i++){
		s+=1.0/(float)i;
	}
	return s;
}

float beta(int n,int i){
	int d1=n-i;
	if (d1==0) return 0;
	return 2.0*n*(a(n+1)-a(i))/((float)d1*((float)d1+1.0)) - 2.0/(float)d1; 
}

float ps(int k, int l, int n, float theta){
	int mkl=min(k,l);
	float common_term=theta*theta*(beta(n,mkl)-beta(n,mkl+1))/2.0;
	/*float first_term=0.0;*/
	if (k==l) {
		return theta*theta*beta(n,k);
	}
	return common_term;
}
//
float sigma(int i, int j,int n){
	int tmp;
	if (i<j) {tmp=i;i=j;j=tmp;}
	if (i==j){ 
		if ( 2*i < n) return beta(n,i+1);
		if ( 2*i == n ) return 2.0*(a(n)-a(i))/(n-i) - 1.0/(i*i);
		return beta(n,i)-1.0/(i*i);
	}
	if (i+j<n) return (beta(n,i+1)-beta(n,i))/2.0;
	if (i+j==n) return  ((a(n)-a(i))/(n-i))+((a(n)-a(j))/(n-j))-(beta(n,i)+beta(n,j+1))/2.0-1.0/(i*j);
	return (beta(n,j)-beta(n,j+1))/2.0-1.0/(i*j);
}

//
float pe(int k, int l, int n, float theta){
	float res=0.0;
	if (k+l>n) return res;
	if (k+l==n){
		res=theta*theta*(((a(n)-a(k))/(float)(n-k))+((a(n)-a(l))/(float)(n-l))-((beta(n,k)+beta(n,l))/2.0));
	} else{
		res=theta*theta*(1.0/(k*l) - (beta(n,k)-beta(n,k+1)+beta(n,l)-beta(n,l+1))/2.0 );
		//printf("res:%f, %f\n",res,(beta(n,k)-beta(n,k+1)+beta(n,l)-beta(n,l+1)));
	}
	return res;
}


float cs(int i, int j, int k, int l, int nx, int ny, int nxy,float theta){
	if ((l-j<0)||(nx-nxy-l+j<0)||(k-i<0)||(ny-nxy-k+i<0)||(nx+ny-nxy-k<0)) return 0;
	float first_factor = exp(lnbinomial(nx-nxy,l-j)+lnbinomial(ny - nxy,k - i)- lnmultinomial(nx+ny-nxy,l,k-l));
	int lb = max(i-nxy,l-j);
	int ub = min3(i,nx-nxy,k-j); 
	int kx;
	float sum=0.0;
	for(kx=lb;kx<=ub;kx++){
		float b1=lnbinomial(nxy,i-kx);
		float b2=lnbinomial(nx-nxy+j-l,kx+j-l);
		float b3=lnbinomial(k-kx,j);
		sum+=exp(b1+b2+b3); 
	}
	return first_factor*sum;
}


float ce(int i, int j, int k, int l, int nx, int ny, int nxy,float theta){
	if ((l-j<0)||(nx-nxy-l+j<0)||(k-i<0)||(ny-nxy-k+i<0)||(nx+ny-nxy-k-l<0)) return 0.0;
	if (k+l>nx+ny-nxy) return 0.0;
	float first_factor = exp(lnbinomial(nx-nxy,l-j)+lnbinomial(ny - nxy,k - i)- lnmultinomial(nx+ny-nxy,k,l)); 
	//fprintf(stdout,"i,j,k,l:%d,%d,%d,%d first factor:%g\n",i,j,k,l,first_factor);
	int lb = max3(0,i-nxy,k-ny+j);
	int ub = min(i,nx-nxy+j-l);
	int kx;
	float sum=0.0;
	//fprintf(stdout,"lb=%d,ub=%d\n",lb,ub);
	for(kx=lb;kx<=ub;kx++){
		float b1=lnbinomial(nxy,i-kx);
		float b2=lnbinomial(nx-nxy+j-l,kx);
		float b3=lnbinomial(ny-k+kx,j);
		//printf("b1=%f,b2=%f,b3=%f\n",exp(b1),exp(b2),exp(b3));
		sum+=exp(b1+b2+b3); 
	}
	return first_factor*sum;
}

float p(int i, int j, int nx, int ny, int nxy, float theta){
	int lb = 1;
	int ub = nx +ny -nxy -1 ;
	float sum = 0.0;
	float term;
	int k,l;
	float c_s,p_s,c_e,p_e;
	for (l=lb;l<=ub;l++){
		for(k=lb;k<=ub;k++){
			if (k>=l){
				c_s=cs(i,j,k,l,nx,ny,nxy,theta);
				p_s=ps(k,l,nx+ny-nxy,theta);
				c_e=ce(i,j,k,l,nx,ny,nxy,theta);
				p_e=pe(k,l,nx+ny-nxy,theta);
				term=c_s*p_s+c_e*p_e;
			}
			else{
				term=cs(j,i,l,k,ny,nx,nxy,theta)*ps(k,l,nx+ny-nxy,theta)+ce(j,i,l,k,ny,nx,nxy,theta)*pe(k,l,nx+ny-nxy,theta);
			}  
			//if (isnan(term)) {printf("%d %d %d %d\n",i,j,k,l);} 
			sum+=term;
		}
	}
	return sum;  
}
float cov_missing(int i, int j, int nx, int ny, int nxy, float theta){
	
	
	return p(i,j,nx,ny,nxy,theta)- theta*theta/(i*j);
}

float hellman_sum(int *conf,int nconf, int noind, int *r);
int check_conf(int j, int n,int* conf);
float pc(int j,int *r,int n);
//
float hn(int n){
	//computes (rather naively) the
	//n-th harmonic number
	int i;
	float s =0.0;
	for(i=1;i<=n;i++) s+=1.0/i;
	return s;       
}

void print_conf_matrix(int *conf,int nconf, int n){
    int i,j;
    for (i=0;i<nconf;i++){
        for(j=0;j<n;j++){
            fprintf(stdout,"%d ",conf[i*n+j]);     
        }
        fprintf(stdout,"\n");
    }
}

float watterson(int* a, int n){
	float an=hn(n-1);
	//count the number of segregating sites
	int k=0;
	int i;
	for(i=0;i<n;i++) if (a[k]==1) {k=1;break;}
	return (float)k/an;             
}

float corrected_watterson(int*a, int* r, int n){
	int j;
	int s=0;
	for(j=0;j<n;j++) if (a[j]==1) {s=1;break;}
	if (0==s) return 0.0;
	float sum_pc=0.0;
	for(j=2;j<=2*n;j++){
		sum_pc+=pc(j,r,n)*hn(j-1);      
	}       
	return (float)s/sum_pc;
}

void store_configurations(int n,int s ,int* conf,int *row){
	/*  
     all possible legal (i.e. sum to s) 
	 contributions of chromosomes coming from the 
	 individuals being studied, get stored in conf.
	 row counts such configurations 
	 */
    int* d = malloc(n*sizeof(int));
	int i,j;
	for(i=0;i<n;i++) d[i]=1;
	int end=0;
    *row=0;
	while(0==end){
		for(i=0;i<n;i++){
			d[i]=d[i]+1;
			if (d[i]>2) d[i]=1; 
			else break;
			if (i==(n-1)) end=1;                            
		}
		/* here store or print, if passes the test */
		if (check_conf(s,n,d)){
            for(j=0;j<n;j++) conf[*row*n+j]=d[j];
            (*row)++;
        }
	}
    free(d);
}
float pc(int j,int *r,int n){
	//probability that exactly j 
	//chromosomes are sequenced in a given position,
	//given the read depths for each individual
    int nconf;
    int *conf=calloc(CONFSIZE,sizeof(int)); 
	store_configurations(n, j , conf,&nconf);
    float s = hellman_sum(conf, nconf, n, r);
    free(conf);
    return s;
}
int check_conf(int j, int n,int* conf){
    /* indicator function in Hellman's formula */
    int i,sum=0;
    for (i=0;i<n;i++) sum+=conf[i];
    return (sum==j);
}
//
float hellman_prod(int noind, int* r,int* i){
    /* product appearing in hellman's formula */
    int q;
    float prod=1.0;
    for (q=0;q<noind;q++){
		prod*=pow(2,1-r[q])+(i[q]-1)*(1-pow(2,2-r[q]));     
    }    
    return prod; 
}
//
float hellman_sum(int *conf,int nconf, int noind, int *r){
	/*
	 conf : legal configurations
	 nconf : how many legal configurations are there (rows of conf)
	 n : number of individuals (columns of conf)
	 r : coverage for individual 
	 */
	int q;
	float s = 0;
	for(q=0;q<nconf;q++){
		s+=hellman_prod(noind,r,conf+q*noind);
	}
	return s;
}
//
float ck(int r){
	//avg number of chromosomes sequenced 
	// for a given individual
	return 2*(1-2*pow(2.0,-r));
}
//
int ak(int m,int r){
	/* 
	 r read depth
	 m number of minor alleles 
	 */
	if (0==m) return 0;
	if (r==m) return 2;
	return 1;
}
//
float fhat(int* a, float*c, int n,float* var){
	/*note : only one for loop to build both sums
	 (minor) allele frequency estimation
	 it also returns the variance */
	int k;
	float sumc=0.0,sum1=0.0,sum2=0.0,tmp=0.0;
	float f;
	for(k=0;k<n;k++){ /* loop over the individuals */
		tmp=1.0/(3.0-c[k]);
		sum1+=tmp*a[k];
		sum2+=tmp;
		sumc+=c[k];     
	}
	f=sum1/(2.0*sum2);
	*var=f*(1.0-f)*sumc/(2.0*tmp*tmp);
	return f;
}
void _log(char *s){
	fprintf(stderr,"log:%s\n",s);
}
//
float tajima_pi(int *a,float *c, int n){
	int k,l;
	float sum1=0.0,sum2=0.0,sum3=0.0,sum4=0.0;
	for(k=0;k<n;k++){
		sum1+=a[k]*(2-a[k]);
		for (l=0;l<n;l++){
			if (k==l) continue;
			sum2+=c[k]*c[l]*a[k]*(2.0-a[l])/4.0;
			sum3+=c[k]*c[l]/2.0;
		}
		sum4+=c[k];
	}
	//fprintf(stderr,"sum1=%f,sum2=%f,sum3=%f,sum4=%f\n",sum1,sum2,sum3,sum4);
	return (sum1+sum2)/(sum4-n+sum3); 
}

float watterson_variance(float theta, int l, int* nx, int* ny,int* nxy, float* var_d, float* var_h){
	float asum=0.0;
	int x;
	for(x=0;x<l;x++){
		asum+=a(nx[x]);
	}
	float covsum=0.0;
	float checksum=0.0;
	float dsum=0.0;
	float covdsum=0.0;
	float checksum_h=0.0;
	float hsum=0.0;
	float covhsum=0.0;
	int i,j;
	float wxd,wyd,wxh,wyh;
	float tmpcov;
	/**/
	float sigma_=0.0;
	/**/
	for(x=0;x<l;x++){
		/*printf("x=%d\r",x);*/
		checksum=0.0;
		checksum_h=0.0;
		for(i=1;i<=nx[x]-1;i++){
			//wx=(float)(2*i*(nx[x]-i))/(float)(l*nx[x]*(nx[x]-1))-1.0/asum;
			wxh=(float)(2*i*(nx[x]-i))/(float)(l*nx[x]*(nx[x]-1))-(float)i/(float)(l*(nx[x]-1));
			wxd=(float)(2*i*(nx[x]-i))/(float)(l*nx[x]*(nx[x]-1))-1.0/asum;
			//printf("w[%d,%d]=%f\n",i,nx[x],wx);
			checksum+=wxd/(float)i;//debug
			dsum+=wxd*wxd/(float)i;
			checksum_h+=wxh/(float)i;//debug
			hsum+=wxh*wxh/(float)i;
			for(j=1;j<=ny[x]-1;j++){
				wyh=(float)(2*j*(ny[x]-j))/(float)(l*ny[x]*(ny[x]-1))-(float)j/(float)(l*(ny[x]-1));
				wyd=(float)(2*j*(ny[x]-j))/(float)(l*ny[x]*(ny[x]-1))-1.0/asum;
				tmpcov=cov_missing(i,j,nx[x],ny[x],nxy[x],theta);
				/**/
				sigma_=sigma(i,j,nx[x]);
				printf("\ncov[%d,%d]=%f\tsigma[%d,%d]=%f",i,j,tmpcov/(theta*theta),i,j,sigma_);
				/**/
				covsum+=tmpcov;
				covdsum+=wxd*wyd*/**/tmpcov/**//*sigma_*theta*theta*/;
				covhsum+=wxh*wyh*tmpcov;
			}
		}
    }
	*var_d=dsum*theta*(float)l+covdsum*(float)l;
	*var_h=hsum*theta*(float)l+covhsum*(float)l;
	return (theta/asum)+(1.0/(asum*asum))*covsum*l;     
}
