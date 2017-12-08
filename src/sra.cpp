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

// REMOVED [[Rcpp::export]]
// Fast median computation
double median_rcpp(NumericVector x) {
  NumericVector y = clone(x);
  int n, half;
  double y1, y2;
  n = y.size();
  half = n / 2;
  if(n % 2 == 1) {
    // median for odd length vector
    std::nth_element(y.begin(), y.begin()+half, y.end());
    return y[half];
  } else {
    // median for even length vector
    std::nth_element(y.begin(), y.begin()+half, y.end());
    y1 = y[half];
    std::nth_element(y.begin(), y.begin()+half-1, y.begin()+half);
    y2 = y[half-1];
    return (y1 + y2) / 2.0;
  }
}

// REMOVED [[Rcpp::export]]
// Fast median average distance
double mad_rcpp(NumericVector x, double scale_factor = 1.4826) {
  // scale_factor = 1.4826; default for normal distribution consistent with R
  return median_rcpp(abs(x - median_rcpp(x))) * scale_factor;
}


NumericVector SdPerRow(NumericMatrix obj, int type) {

  NumericVector res(obj.nrow());
  
  for (int i=0; i<obj.nrow(); i++) {
    if (type == 1) {
      // MAD
      // res[i] = mean(abs(obj.row(i) - mean(obj.row(i))));
      res[i] = mad_rcpp(obj.row(i), 1);
    }
    else {
      res[i] = var(obj.row(i));
    }
  }
 
  return(res);
}


//' Compute the sequential rank agreement between k ranked lists
//' 
//' @description Computes the sequential rank agreement (number of items present in all k lists divided by the current rank) for each rank in the k lists
//' @param rankMat A matrix with k columns corresponding to the k ranked lists. Elements of each column are integers between 1 and the length of the lists
//' @param maxlength The maximum depth that are needed XXX
//' @param B The number of resamples to use in the presence of censored lists
//' @param cens A vector of integer values that
//' @param type The type of distance measure to use: 0 (the default) is the variance while 1 is MAD (median absolute deviation)
//' @param epsilon A non-negative numeric vector that contains the minimum limit in proportion of lists that must show the item. Defaults to 0. If a single number is provided then the value will be recycles to the number of items.
//' @return A vector of the same length as the rows in rankMat containing the squared (!) sequential rank agreement between the lists for each depth. If the MAD type was chosen then the sequential MAD values are returned
//' @author Claus Ekstrøm <ekstrom@@sund.ku.dk>
//' @encoding UTF-8 
// [[Rcpp::export]]
NumericVector sracpp(IntegerMatrix rankMat, int maxlength, int B, IntegerVector cens, int type=0, NumericVector epsilon=NumericVector::create(0)) {

  // The number of lists
  int nLists = rankMat.ncol();
  
  NumericVector useEpsilon(maxlength);;
  // Check length  
  if ((epsilon.size() != maxlength) && (epsilon.size() != 1))
    stop("Length of epsilon must be 1 or the number of items");
  if (is_true(any( epsilon<0 )))
    stop("Values for epsilon must be non-negative");

  if (epsilon.size()==1) {
    std::fill(useEpsilon.begin(), useEpsilon.end(), epsilon(0));
  } else {
    useEpsilon = epsilon;
  }

  

  // Sanity check of type
  if (type != 1)
    type = 0;

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
      NumericVector itemWeight(maxlength);
      
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
      itemMetric = SdPerRow(itemRank, type);

      // Now compute the agreement for each depth
      for (int depth=0; depth<maxlength; depth++) {
	
	for (int l = 0; l < nLists; l++) {
	  // Add variable to keep
	  // -1 below because lists starts from 0 in C++
	  // keepDepth[x(depth, l)-1] = TRUE;
	  itemWeight[x(depth, l)-1] += 1.0;
	}

	keepDepth = (itemWeight > useEpsilon*nLists);
	
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


//' Compute the sequential rank agreement between k ranked lists
//' 
//' @description Computes the sequential rank agreement (number of items present in all k lists divided by the current rank) for each rank in the k lists
//' @param rankMat A matrix with k columns corresponding to the k ranked lists. Elements of each column are integers between 1 and the length of the lists
//' @param type The type of distance measure to use: 0 (the default) is the variance while 1 is MAD (mean absolute deviation)
//' @param epsilon A non-negative numeric vector that contains the minimum limit in proportion of lists that must show the item. Defaults to 0. If a single number is provided then the value will be recycles to the number of items.
//' @return A vector of the same length as the rows in rankMat containing the sequential rank agreement between the lists for each depth (squared for type=0)
//' @author Claus Ekstrøm <ekstrom@@sund.ku.dk>
//' @encoding UTF-8 
// [[Rcpp::export]]
List sracppfull(IntegerMatrix rankMat, int type=0, NumericVector epsilon=NumericVector::create(0)) {

  // The number of lists
  int nLists = rankMat.ncol();
  int maxlength=rankMat.nrow();
  NumericVector useEpsilon(maxlength);
  IntegerVector whenIncluded(maxlength, maxlength);

  // Check too big
  // Check length  
  if ((epsilon.size() != maxlength) && (epsilon.size() != 1))
    stop("Length of epsilon must be 1 or the number of items");
  if (is_true(any( epsilon<0 )))
    stop("Values for epsilon must be non-negative");

  if (epsilon.size()==1) {
    std::fill(useEpsilon.begin(), useEpsilon.end(), epsilon(0));
  } else {
    useEpsilon = epsilon;
  }

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

  // Input validation
  if (type != 1)
    type = 0;

  // Check that dim of x is less than maxlength
  
  // Check for no duplicates in each list
  
  NumericMatrix itemRank(maxlength, nLists);
  NumericVector itemMetric(maxlength);
  NumericVector res;
  
  //    IntegerVector iv = seq_len(maxlength);
  
  // The set of variables to keep per depth
  LogicalVector keepDepth(maxlength, FALSE);
  NumericVector itemWeight(maxlength);
  
  for (int l = 0; l < nLists; l++) {
    itemRank(_, l) = order_( x( _ , l) );
  }
  itemMetric = SdPerRow(itemRank, type);
  
  // Now compute the agreement for each depth
  for (int depth=0; depth<maxlength; depth++) {
    
    for (int l = 0; l < nLists; l++) {
      // Add variable to keep
      // -1 below because lists starts from 0 in C++
      //      keepDepth[x(depth, l)-1] = TRUE;
      itemWeight[x(depth, l)-1] += 1.0;  // Should be 1/maxlength but multiply that below
      if (itemWeight[x(depth, l)-1] > useEpsilon[depth]*nLists) {
	keepDepth[x(depth, l)-1] = TRUE;
	whenIncluded[x(depth, l)-1] = std::min(whenIncluded[x(depth, l)-1], depth+1);
      }
    }

    // Rcout << keepDepth << "    -->   for depth " << (depth+1) << std::endl ; 
    
    // Now weigh the metrics together
    res = itemMetric[keepDepth];
    returnVector[depth] = mean(res);
  }

  return List::create(Rcpp::Named("sra")=returnVector,
		      Rcpp::Named("whenIncluded")=whenIncluded
		      );
  
  
  //  return(returnVector);
}



