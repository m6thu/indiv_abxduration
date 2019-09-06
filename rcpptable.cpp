#include <Rcpp.h>

// Save on the typing... 
typedef std::pair<double, int>  ptype;


// A comparison function to rank values in descending order
bool compare_values(const ptype &p1, const ptype &p2)
{
    return p1.second > p2.second;
}

// Get the top number of observations
// [[Rcpp::export]]
Rcpp::List table_cpp(const Rcpp::NumericVector & v, bool sort_data = false)
{
    
    // Create a map
    std::map<double, int> Elt;
    
    Elt.clear();
    
    // Fill the map with occurrences per number.
    for (int i = 0; i != v.size(); ++i) {
        Elt[ v[i] ] += 1;
    }
    
    // Get how many unique elements exist... 
    unsigned int n_obs = Elt.size();
    
    // Switch map to a vector so that we can sort by value
    std::vector< ptype > sorted_Elt(Elt.begin(), Elt.end());
    
    if(sort_data){
        // Perform the sort with a custom sort function.
        std::sort(sorted_Elt.begin(), sorted_Elt.end(), compare_values);
    }
    // Else, return. 
    
    // Stop here if you do not need to import into R.
    // Why? There is no ability to export a set w/ a pair into R. *cries*
    
    // Recast for R using Rcpp::*Vectors to avoid another copy)      
    Rcpp::NumericVector result_keys(n_obs);
    Rcpp::IntegerVector result_vals(n_obs);
    
    unsigned int count = 0;
    
    // Need to use iterators to access objects
    for( std::vector< ptype >::iterator it = sorted_Elt.begin(); it != sorted_Elt.end(); ++it )
    {
        // Move them into split vectors
        result_keys(count) = it->first;
        result_vals(count) = it->second;
        
        count++;
    }
    
    return Rcpp::List::create(Rcpp::Named("lengths") = result_vals,
                              Rcpp::Named("values") = result_keys);
}