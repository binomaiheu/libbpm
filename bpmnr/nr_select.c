/**
  @file
  @ingroup nr
*/
#include <bpm/bpm_messages.h>
#include <bpm/bpm_nr.h>

/**
   Find the kth largest element of the array after sorting. Nicked from numerical recipes, C8.5, p342
   See: http://www.library.cornell.edu/nr/cbookcpdf.html
   
   @return The value of the median element
*/
double nr_select( int k, int n, double *org_arr ) {

  unsigned long i,ir,j,l,mid;
  double a,tempr;
  double *arr;

  if( ! org_arr ) {
    bpm_error( "Invalid array in nr_select(...)", 
	       __FILE__, __LINE__ );
    return -DBL_MAX;
  }

  // First, create a copy of the array so this one doesn't get mashed up
  arr = malloc( (n+1) * sizeof(double) );
  memcpy( &(arr[1]), org_arr, n * sizeof(double) );

  // Now, perform the sort to find the selected element
  l=1;
  ir=n;
  for (;;) {

    if (ir <= l+1) {
      if (ir == l+1 && arr[ir] < arr[l]) {
	SWAP(arr[l],arr[ir]);
	  }
      return arr[k];
    } else {
      mid=(l+ir) >> 1;
      SWAP(arr[mid],arr[l+1]);
	if (arr[l] > arr[ir]) {
	  SWAP(arr[l],arr[ir]);
	    }
      if (arr[l+1] > arr[ir]) {
	SWAP(arr[l+1],arr[ir]);
	  }
      if (arr[l] > arr[l+1]) {
	SWAP(arr[l],arr[l+1]);
	  }
      i=l+1;
      j=ir;
      a=arr[l+1];
      for (;;) {
	do i++; while (arr[i] < a);
	do j--; while (arr[j] > a);
	if (j < i) break;
	SWAP(arr[i],arr[j]);
      }
      arr[l+1]=arr[j];
      arr[j]=a;
      if (j >= k) ir=j-1;
      if (j <= k) l=i;
    }
  }

  return 0;
}

	     
