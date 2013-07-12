#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <time.h>
#include <setjmp.h>

int main(void) {
	int m, l, N, P, i, j, nproc, myproc, beg, end, mmin, NperProc, last_nonzero, check_freq;
	double start_time;
	unsigned *a, *b, *c;
	MPI_Init(NULL, NULL);
	MPI_Comm_rank(MPI_COMM_WORLD, &myproc);
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);
	mmin = (int) (log(nproc) / log(2) + 0.5);
	mmin = 1 + 3 * mmin / 2;
	mmin += 4 - (mmin % 3);
	for (m = mmin;; m += 3) {
		start_time = MPI_Wtime();
		N = 1 << m;
		NperProc = N / nproc;
		beg = myproc * NperProc;
		end = (myproc + 1) * NperProc - 1;
		P = 3 * N;
		a = (unsigned*) malloc(sizeof(unsigned) * (NperProc + 2));
		b = (unsigned*) malloc(sizeof(unsigned) * (NperProc + 2));
		a[0] = a[NperProc + 1] = 0;
		b[0] = b[NperProc + 1] = 0;
		a -= beg - 1;
		b -= beg - 1;
		i = (m - 1) / 3;
		for (j = beg; j <= end; j++) {
			if (j == 0) {
				a[j] = 1;
			} else {
				a[j] = 0;
			}
		}
		check_freq = P / 10000;
		if (check_freq == 0) {
			check_freq = 1;
		}
		for (l = 0; l < P; l++) {
			if (myproc % 2 == 0) {
				if (myproc - 1 >= 0) {
					MPI_Sendrecv(a + beg, 1, MPI_INTEGER, myproc - 1, 0, a + beg - 1, 1, MPI_INTEGER, myproc - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
				}
				if (myproc + 1 < nproc) {
					MPI_Sendrecv(a + end, 1, MPI_INTEGER, myproc + 1, 0, a + end + 1, 1, MPI_INTEGER, myproc + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
				}
			} else {
				if (myproc + 1 < nproc) {
					MPI_Sendrecv(a + end, 1, MPI_INTEGER, myproc + 1, 0, a + end + 1, 1, MPI_INTEGER, myproc + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
				}
				if (myproc - 1 >= 0) {
					MPI_Sendrecv(a + beg, 1, MPI_INTEGER, myproc - 1, 0, a + beg - 1, 1, MPI_INTEGER, myproc - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
				}
			}
			for (j = beg; j <= end; j++) {
				b[j] = -j * a[j - 1] + (j - 1) * a[j] + a[j + 1];
			}
			c = b;
			b = a;
			a = c;
			if (myproc == 0 && (l % check_freq) == 0) {
				printf("%6.2f%c done\r", 100.0 * (l + 1) / P, '%');
				fflush(stdout);
			}
		}
		if (myproc == 0) {
			printf("m=%i i=%i    \n", m, i);
			for (j = 0; j < 1 << (i + 2); j++) {
				if (j == 0) {
					a[0]--;
				}
				a[j] /= 1 << (m - i);
				a[j] %= 1 << (2 * i + 2);
			}
			last_nonzero = 0;
			for (j = 0; j < 1 << (i + 2); j++) {
				if (a[j] != 0) {
					last_nonzero = j;
				}
			}
			for (j = 0; j <= last_nonzero; j++) {
				printf("%4u ", a[j]);
			}
			printf("\n");
		}
		a += beg - 1;
		b += beg - 1;
		free(a);
		free(b);
		if (myproc == 0) {
			printf("%.2f seconds\n\n",  (MPI_Wtime() - start_time));
		}
	}
	MPI_Finalize();
}
