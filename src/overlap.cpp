#include <Rcpp.h>
using namespace Rcpp;

//' Compute the overlap between k ranked lists
//' 
//' @description Computes the overlap (number of items present in all k lists divided by the current rank) for each rank in the k lists
//' @param rankMat A matrix with k columns corresponding to the k ranked lists. Elements of each column are integers between 1 and the length of the lists
//' @return A vector of the same length as the rows in rankMat containing the overlap between the lists for each rank
//' @author Claus Ekstrøm <ekstrom@@sund.ku.dk>
//' @export
// [[Rcpp::export]]
NumericVector overlap(IntegerMatrix rankMat) {

  // The number of lists
  int nLists = rankMat.ncol();
  int nItems = rankMat.nrow();
  NumericVector res(nItems);

  LogicalMatrix items(nItems, nLists);
  IntegerVector pick(nLists);
  IntegerVector lookat;
  
  int i, j, k;
  double count=0;
  bool keep;
 
  for (i=0; i<nItems; i++) {

    // Get the items in the ith row
    for (j=0; j<nLists; j++) {
      pick(j) = rankMat(i, j)-1;	    
      items(pick(j), j) = TRUE;
    }
    // Now run through the unique items that were updated

    lookat = unique(pick);

    for (j=0; j<lookat.length(); j++) {
      keep=TRUE;
      for (k=0; k<nLists; k++) {
	keep = (keep & items(lookat(j), k));
	if (!keep)
	  break;
      }
      if (keep) {
      	count++;
      }
    }

    res(i) = count / (i+1);    
  } 
    
    return(res);
}


