/* msort_Kimchi() */
static int *msort_Kimchi_tmp;
static double *msort_Kimchi_x;
void msort_Kimchi_Merge(int base[], int n_left, int n_right)
{
    register int ind = 0, i_l = 0, i_r = 0;
    register int *c, *d;
    int *msort_Kimchi_tmp_r = msort_Kimchi_tmp+n_left;
    d = base; c = msort_Kimchi_tmp;
loop1:
    REPEAT256( if (d==base+n_left+n_right) goto exit1; *(c++) = *(d++); );
    goto loop1;
exit1:
    while (i_l < n_left && i_r < n_right)
	base[ind++] = (msort_Kimchi_x[msort_Kimchi_tmp[i_l]] <
		       msort_Kimchi_x[msort_Kimchi_tmp_r[i_r]]) ?
	    msort_Kimchi_tmp[i_l++] : msort_Kimchi_tmp_r[i_r++];
    if (i_l < n_left)
    {
	d = msort_Kimchi_tmp+i_l; c = base+ind;
    loop2:
REPEAT256( if (d==msort_Kimchi_tmp+n_left) goto exit2; *(c++) = *(d++); );
goto loop2;
    }
    else
    {
	d = msort_Kimchi_tmp_r+i_r; c = base+ind;
    loop3:
REPEAT256( if (d==msort_Kimchi_tmp_r+n_right) goto exit2; *(c++) = *(d++); );
goto loop3;
    }
exit2:
    return;
}

void msort_Kimchi (int N,  double x[], int idx[], int option)
{
    int len, *base;
    if (option==USE_NEW_IDX) Sequentially_index(N,idx);
    msort_Kimchi_tmp = MImem(N);
    msort_Kimchi_x = x;
    for (len = 1; len < N; len *= 2)
    {
	for (base = idx; base+len < idx+N; base += 2*len)
	{
	    msort_Kimchi_Merge (base,len, MIN(len,(idx+N-(base+len))));
	}
    }
    MFREE (msort_Kimchi_tmp);
    return;
}
/* end msort_Kimchi() */
