/* Makefile:
bad:
        g++ main.cpp /usr/lib/liblapack.a /usr/lib/libblas.a -lgfortran -z muldefs -pthread
good:
        g++ -DGOOD main.cpp /usr/lib/liblapack.a /usr/lib/libblas.a -lgfortran -z muldefs -pthread
*/
/* Output I see:
$ make good
g++ -DGOOD main.cpp /usr/lib/liblapack.a /usr/lib/libblas.a -lgfortran -z muldefs -pthread
$ ./a.out 
Info1,2 = 0,0
-7050.245589,-8479.237956
-8580.257527,-9209.117716
-6301.087406,-6399.428666
$ ./a.out 
Info1,2 = 0,0
-7050.245589,-8479.237956
-8580.257527,-9209.117716
-6301.087406,-6399.428666
$ make bad
g++ main.cpp /usr/lib/liblapack.a /usr/lib/libblas.a -lgfortran -z muldefs -pthread
$ ./a.out 
Info1,2 = 0,0
-7050.355400,-8478.864949
-8457.115555,-8501.374698
1703.121368,-3150.274644
$ ./a.out 
Info1,2 = 0,0
-7106.927073,-8501.507914
-8580.257527,-9209.117716
-6299.486160,-6398.683343
*/

#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <pthread.h>
#include <complex>

// Global pthreads stuff
pthread_mutex_t g_mutex = PTHREAD_MUTEX_INITIALIZER;
int count = 0;
pthread_cond_t g_cond = PTHREAD_COND_INITIALIZER;

typedef std::complex<double> complex_t;

extern "C" float slamch_(const char *);
extern "C" double dlamch_(const char *);
void lapack_init(){
	slamch_("E");
	dlamch_("E");
}

struct thread_data{
	complex_t *A;
	complex_t *q;
	int n, info;
};
extern "C" void zgeev_(
	const char *jobvl, const char *jobvr,
	const int &n, complex_t *a, const int &lda,
	complex_t *w,
	complex_t *vl, const int &ldvl,
	complex_t *vr, const int &ldvr,
	complex_t *work, const int &lwork, double *rwork, int *info);
void* thread_func(void *data){
	thread_data *T = (thread_data*)data;
	const int n = T->n;
	int lwork = 2*n;

	// Query workspace
	lwork = -1;
	zgeev_("N", "N", n, T->A, n, T->q,
		NULL, 1, NULL, 1, T->q, lwork, NULL, &T->info);
	lwork = (int)T->q[0].real();

	complex_t *work = (complex_t*)malloc(sizeof(complex_t)*lwork);
	double *rwork = (double*)malloc(sizeof(double)*2*n);
	
	volatile int t;
	// Sync threads at this point
	pthread_mutex_lock(&g_mutex);
	++count;
	pthread_mutex_unlock(&g_mutex);
	do{
		pthread_mutex_lock(&g_mutex);
		t = count;
		pthread_mutex_unlock(&g_mutex);
	}while(t < 2);

#ifdef GOOD
	pthread_mutex_lock(&g_mutex);
#endif
	zgeev_("N", "N", n, T->A, n, T->q,
		NULL, 1, NULL, 1, work, lwork, rwork, &T->info);
#ifdef GOOD
	pthread_mutex_unlock(&g_mutex);
#endif

	// Sync threads at this point
	pthread_mutex_lock(&g_mutex);
	++count;
	pthread_mutex_unlock(&g_mutex);
	do{
		pthread_mutex_lock(&g_mutex);
		t = count;
		pthread_mutex_unlock(&g_mutex);
	}while(t < 4);
	
	free(rwork);
	free(work);
	return NULL;
}


int main(int argc, char *argv[]){
	int n = 2*71;
	
	lapack_init();
	
	pthread_t thread1, thread2;
	thread_data data1, data2;
	complex_t *A1, *A2, *q1, *q2;

	data1.n = n;
	data2.n = n;
	data1.A = (complex_t*)malloc(sizeof(complex_t)*n*n);
	data2.A = (complex_t*)malloc(sizeof(complex_t)*n*n);
	data1.q = (complex_t*)malloc(sizeof(complex_t)*n);
	data2.q = (complex_t*)malloc(sizeof(complex_t)*n);

	// Load in data
	FILE *fp;
	fp = fopen("1.dat", "rb");
	fread(data1.A, n*n*sizeof(complex_t), 1, fp);
	fclose(fp);
	fp = fopen("1.dat", "rb");
	fread(data2.A, n*n*sizeof(complex_t), 1, fp);
	fclose(fp);

	pthread_create(&thread1, NULL, thread_func, &data1);
	pthread_create(&thread2, NULL, thread_func, &data2);
	pthread_join(thread1, NULL);
	pthread_join(thread2, NULL);

	printf("Info1,2 = %d,%d\n", (int)data1.info, (int)data2.info);
	for(int i = 0; i < 3; ++i){
		complex_t q = data1.q[i];
		printf("%f,%f\n", q.real(), q.imag());
	}
	
	free(data1.A);
	free(data2.A);
	free(data1.q);
	free(data2.q);
	
	return 0;
}


