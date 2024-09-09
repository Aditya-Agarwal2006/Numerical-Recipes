#define WANT_STREAM
#include "nr.h"
#include "nrutil_nr.h"

#ifdef use_namespace
using namespace NEWMAT;
#endif


void print_array(const ColumnVector &a, const int n_per_row, const int col_width)
{
	int i,j=0;

	int n_elements=a.size();
	for (i=0; i < n_elements; i++) {
		if (j == n_per_row) {
			cout << endl;
			j=0;
		}
		cout << setw(col_width) << a[i];
		j++;
	}
	cout << endl;
}
