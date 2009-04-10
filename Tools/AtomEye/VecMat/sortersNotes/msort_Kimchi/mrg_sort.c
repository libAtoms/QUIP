/*************************************************************/
/*    M E R G E     S O R T     (non-recursive)              */
/* Send errors to Yechiel M. Kimchi  (Dec. 1997)             */
/*************************************************************/
#include <stdio.h>
#include <stdlib.h>	/* For exit() */
#define N 100

void
FatalError(char *msg)
{ fprintf(stderr,"Fatal error: %s\n", msg);
  exit(1);
}

int
min(int n, int m)
{ return (n < m ? n : m);
}


void
ReadArray(int a[], int *pn)
{ int i;

  if (scanf("%d",pn) < 1) FatalError("Could not read number of elements");
  if (*pn > N) FatalError("You ask me to read too many elements");
  for (i = 0; i < *pn; ++i) {
    if (scanf("%d",&a[i]) < 1) FatalError("Could not read an item");
  }
  return;
}


void
PrintArray(int a[], int n, char *msg)
{  int i;

   for (i = 0; i < n-1; ++i){
      printf("%d,",a[i]);
   }
   printf("%d   %s\n",a[n-1], msg);
   return;
}


void
swap(int *px, int *py)
{  int z = *px;

   *px = *py;
   *py = z;
}


void
CopyArray(int source[], int size, int target[])
{ int indx;

  for (indx = 0; indx < size; indx++) {
    target[indx] = source[indx];
  }
  return;
}
void
Merge(int base[], int n_left, int n_right)
{ int ind = 0, i_l = 0, i_r = 0, tmp[N], *tmp_r = tmp + n_left;

  CopyArray(base, n_left+n_right, tmp);
  while (i_l < n_left && i_r < n_right) {
    base[ind++] = (tmp[i_l] < tmp_r[i_r] ? tmp[i_l++] : tmp_r[i_r++]);
  }
  if (i_l < n_left) {
    CopyArray(tmp+i_l, n_left-i_l, base+ind);
  }else {
    CopyArray(tmp_r+i_r, n_right-i_r, base+ind);
  }
#if 1
  PrintArray(base, n_left+n_right, "sorted part");
#endif
  return;
}


void
MergeSort(int a[], int n)
{ int len, *base;

  for (len = 1; len < n; len *= 2) {
    for (base = a; base+len < a+n; base += 2*len) {
      Merge(base,len, min(len,(a+n - (base+len))));
    }
  }
  return;
}


int
main(void)
{ int n, a[N];

  ReadArray(a, &n);
  PrintArray(a,n, "before sort");
  MergeSort(a,n);
  PrintArray(a,n, "after sort");
  return 0;
}



#if 0
---------------------------- A sample run --------------------------------------

7,4,2,5,9,6,5,3,1,-1,0,-4,-7   before sort
4,7   sorted part
2,5   sorted part
6,9   sorted part
3,5   sorted part
-1,1   sorted part
-4,0   sorted part
2,4,5,7   sorted part
3,5,6,9   sorted part
-4,-1,0,1   sorted part
2,3,4,5,5,6,7,9   sorted part
-7,-4,-1,0,1   sorted part
-7,-4,-1,0,1,2,3,4,5,5,6,7,9   sorted part
-7,-4,-1,0,1,2,3,4,5,5,6,7,9   after sort
#endif
