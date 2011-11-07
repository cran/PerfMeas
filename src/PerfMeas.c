#include <stdlib.h>

/* 
  prec_recall: Computes precision at all possible recall levels
  INPUT:
  precision: pointer to the precision vector 
  recall: pointer to the recall vector
  order: pointer to the vector storing the order of labels
  labels: pointer to the labels vector
  n : pointer to the number of elements of the above vectors
  OUTPUT
  precision and recall vectors
*/
void  prec_recall(double * precision, double * recall, int * order, int * labels, int * n){
  register int i;
  int TP = 0;
  int TN = 0;
  int FP = 0;
  int FN = 0;
  int np = 0;
  int m;
  
  
  m = *n;
  for (i=0; i < m; i++)
    if (labels[i] == 1) np++;  
  TP = np;
  TN = 0;
  FP = m - np;
  FN = 0;
  
  if((TP + FP) > 0)		
	  precision[m-1] = TP/((double)(TP + FP));
  else
	  precision[m-1] = 0;	
  if((TP + FN) > 0)		
	  recall[m-1] = TP/((double)(TP + FN));
  else
	  recall[m-1] = 0;
  	  
  for (i=(m-2); i >=0; i--) {
    if(labels[order[i+1]-1] == 1 ){
		TP--;
		FN++;			
	} else {
		TN++;
		FP--;
	}
    if(TP + FP != 0)
		precision[i] = TP/((double)(TP + FP));
	else
		precision[i] = 0;	
	if(TP + FN != 0)
		recall[i] = TP/((double)(TP + FN));
	else
		recall[i] = 0;	
  } 
}
