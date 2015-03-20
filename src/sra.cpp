#include <Rcpp.h>
using namespace Rcpp;


// wrapper around R's RNG such that we get a uniform distribution over
// [0,n) as required by the STL algorithm
inline int randWrapper(const int n) { return floor(unif_rand()*n); }

IntegerVector randomShuffle(IntegerVector a) {
    // already added by sourceCpp(), but needed standalone
    RNGScope scope;             

    // clone a into b to leave a alone
    IntegerVector b = clone(a);

    std::random_shuffle(b.begin(), b.end(), randWrapper);

    return b;
}


IntegerVector order_(NumericVector x) {
  // Have removed the 3 lines below because checking should be done elsewhere
  //  if (is_true(any(duplicated(x)))) {
  //   Rf_warning("There are duplicates in 'x'; order not guaranteed to match that of R's base::order");
  // }
  NumericVector sorted = clone(x).sort();
  return match(sorted, x);
}

double square( double x ){
    return x*x;
}


NumericVector SdPerRow(NumericMatrix obj) {

  NumericVector res(obj.nrow());
  
  for (int i=0; i<obj.nrow(); i++) {
    res[i] = var(obj.row(i));
  }
 
  return(res);
}


//[[Rcpp::export]]
NumericVector sra(IntegerMatrix rankMat, int maxlength, int B, IntegerVector cens) {

  // The number of lists
  int nLists = rankMat.ncol();

  // Need a numeric copy of this for the sort function. Doesn't really make sense but
  // keeping it for now 
  NumericMatrix x(maxlength, nLists);
  for (int l=0; l<nLists; l++) {
    for (int i=0; i<std::min(maxlength, rankMat.nrow()); i++) {
      x(i, l) = rankMat(i,l);
    }
  }
  
 
  // Keep the average result and the individual permutations
  NumericVector returnVector(maxlength);
  NumericMatrix result(maxlength, B);

    // Input validation

    // Check that censoring elements is at least one. Or possibly not?

    // Check that the censoring is less than maxlength
  if(max(cens)>rankMat.nrow()) {
    stop("You cannot censor a list later than you have data");
  }
  if(max(cens)>maxlength) {
    stop("You cannot censor a list later than the maximum of the data");
  }


    // Check that dim of x is less than maxlength

    // Check that censoring length is the same as the number of lists
    if (cens.size() != nLists) {
      stop("Number of columns in matrix and length of censoring vector does not match");
    }
    // Check for no duplicates in each list

    LogicalMatrix missingItems(maxlength, nLists);
    NumericMatrix itemRank(maxlength, nLists);
    NumericVector itemMetric(maxlength);
    NumericVector res;

    // Find out which items are missing in the individual lists
    for (int l = 0; l < nLists; l++) {
      // Could be smart and do the filling depending on which tail is shortest
      for (int i = 0; i < cens[l]; i++) {
	missingItems( rankMat(i, l), l) = TRUE;
      }
      missingItems(_, l) = ! missingItems(_, l);
    }


    IntegerVector iv = seq_len(maxlength);

    // Iterate over the number of replications used to average the results from
    // missing data
    for (int b=0; b<B; b++) {

      // The set of variable to keep per depth
      LogicalVector keepDepth(maxlength, FALSE);
      
      // Start by creating a permutation of any missing items for each list
      for (int l = 0; l < nLists; l++) {
	// Check to see if we need to make a permutation for list l
	if (cens[l]==maxlength)
	  continue;

	// Create a permutation if necessary
	IntegerVector permValues = randomShuffle(iv[missingItems(_,l)]);

	// Insert the permuted values
	for (int i = cens[l]; i < maxlength; i++) {
	  x(i,l) = permValues[i-cens[l]];
	}
      }

      for (int l = 0; l < nLists; l++) {
	itemRank(_, l) = order_( x( _ , l) );
      }
      itemMetric = SdPerRow(itemRank);

      // Now compute the agreement for each depth
      for (int depth=0; depth<maxlength; depth++) {
	
	for (int l = 0; l < nLists; l++) {
	  // Add variable to keep
	  // -1 below because lists starts from 0 in C++
	  keepDepth[x(depth, l)-1] = TRUE;	  
	}
		
	// Now weigh the metrics together
	res = itemMetric[keepDepth];
	result(depth, b) = mean(res);
      }
      
    } 
    
    for (int i=0; i<maxlength; i++) {
      returnVector[i] = sum(result(i,_))/B;
    }
    
    return(returnVector);
}


