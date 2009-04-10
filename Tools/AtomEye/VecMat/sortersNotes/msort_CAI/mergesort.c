


/* 
 * Sorting Non-Recursive MergeSort of arrays of any size using the power  
 * of 2 trick 
 */

 #include <stdio.h>
 #include <stdlib.h>
 #include <math.h>
 #define SIZE 23

 void main(void)
 {
    void merge_sort(int key[], int n);
    void merge(int *, int *, int *, int, int);
    int power2(double, long int binary[200]);

    int i, key[SIZE] = {99,9,5,32,15,6,3,21,11,65,12,67,34,4,3,76,34,55,
                        0,17,1,8,7};

    printf("These are the UnSorted keys.\n");
    for(i = 0; i < SIZE; ++i)
        printf("%d", key[i]);

    merge_sort(key, SIZE);

    printf("\nThese are the Sorted keys.\n");
    for(i = 0; i < SIZE; ++i)
        printf("%d", key[i]);
 }

/* 
 * MergeSorts an array by using the merge function and one other
 * array as storage array. Note this version will MergeSort an array of    
 * any size by using the fact that we can write any  integer to the base 
 * 2
 */

 void merge_sort(int key[], int n)
 {
    int i, j, k, m, *temp_store, numb_of_two, shift;
    long int binary[200];
    void merge(int *, int *, int *, int, int);
    int power2(double numb, long int binary[200]);
    
    /* find sums of powers of two and find the number of sums it took */
    numb_of_two = power2((double)n, binary);

    /* Note you shouldn't have to cast here but my compiler made me */
    temp_store = (int*) calloc(n, sizeof(int));

    shift = 0;
    /* produce numb_of_two sub arrays all sorted and ready to be merged */

    for(i = 0; i < numb_of_two; ++i){
    /* the shift variable ensures we start at the correct place when 
       merging */
        if(i > 0) 
            shift += binary[i-1];
        for(k = 1; k < binary[i]; k *= 2){
            for(j = 0; j < binary[i] - k; j += 2 * k)
                merge(key + j + shift, key + shift +j + k,temp_store + j 
                      + shift, k, k);

            for(j = 0; j < binary[i]; ++j)
                key[j+shift] = temp_store[j+shift];

     /* The printf statement and for loop below is to show
        the user what is going  on. It has no effect on the Sort */

            printf("%s%s%d%s%d%s%d%s",
               "\nThese are the keys after MergeSorting. Working with 
                the sub", "arrays starting\nat", shift,"and ending at ", 
                j + shift,"with step sizes of ", k, ".\n");
            for(j = 0; j < SIZE; ++j)
                printf("%d", key[j]);
            }
    }

/*  
 * Do the final MergeSort that brings all the individuaL power of two   
 * sub arrays  together. 
 * Note we only go as far as numb_of_two - 1(numb_of_two is number of 
 * sorted sublists to be merged) because it takes n-1 merges to merge n 
 * sub-arrays e.g. 2 sub-arrays takes 1 merge, 3 take 2 etc.
 */

    shift = 0;
    for(k = 0; k < numb_of_two - 1; ++k){
        shift += binary[k];
        merge(key, key + shift, temp_store, shift, binary[k+1]);
        for(j = 0; j < shift + binary[k+1]; ++j)
            key[j] = temp_store[j];

    /* Note the print statement below has nothing to do with sort but 
       it may help to explain it */
    
    printf("%s%d%s%d\n",
               "\nHave just merged sorted sub-arrays with lengths",
               shift, "and", binary[k+1]);
        for(j = 0; j &ltn; ++j)
            printf("%d", key[j]);
    }
    free(temp_store);
 }

/* 
 * Function to merge array1, length n and array2, length m into
 * destination array dest, without using Sentinels 
 */

 void merge(int array1[], int array2[], int dest[], int m, int n)
 {
    int i, j, k;

    i = 0; j = 0; k = 0;

    /* put together the two arrays in order */
    while(i < m && j < n)
        dest[k++] = (array1[i] < array2[j]) ? array1[i++]: array2[j++];

    /* Insert any remaining elements into c */
    while(i < m)
        dest[k++] = array1[i++];
    while(j < n)
        dest[k++] = array2[j++];
 }

/* 
 * Compute an integer as a summation of powers of two and return
 * the number of powers of two it took. Will work for integers as 
 * large as 999999999. We use the modf function. The modf function
 * returns the fractional part of a floating point number and stores the  
 * integer part in integer 
 */

 int power2(double numb, long int binary[200])
 {
    int n, i, count;
    double fraction, integer;

    i = 0;  n = 0;  count = 0;

    while(numb != 0){
        fraction = modf(numb /= 2, &integer);
        numb = integer;
        if(fraction == 0.5){
            binary[i++] = pow(2.0, (float)n);
            ++count;
            ++n;
        }
        else
            ++n;
    }
    return(count);
 }
