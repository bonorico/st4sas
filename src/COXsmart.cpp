#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;
using namespace Rcpp;




mat subdata_arma(const mat& x, const vec& V, const sword vidx){
  mat y = x.rows(find(V == vidx));
  return y ;
  }



// likelihood engine

vec CoxengC(const vec& beta, const mat& data){

  mat sortedd = data.rows(sort_index(data.col(0))) ;
 
  mat x = sortedd.cols(2, sortedd.n_cols-1);

  vec a = exp( x*beta  );

  vec riskset = log( flipud( cumsum( flipud( a ) ) ) );

  vec den = sum( riskset.rows( find( sortedd.col(1)==1 ) ) );

  return den;

    } 
 





// [[Rcpp::export]]

vec lCox(const vec& beta, const vec& y, const mat& x , const vec& str){

 vec str_levs = unique(str);

   int nstr = str_levs.n_elem;

   vec stratasets(nstr);

   for(int sl=0; sl<nstr; ++sl){

     mat Xsl = subdata_arma( x, str, str_levs(sl));

     stratasets(sl) = as_scalar( CoxengC( beta, Xsl ) );

   }

   return  beta.t()*y - sum( stratasets ) ;

  }



//


mat risksetC(const mat& X){

  return flipud(cumsum(flipud(X))) ;

 }



// cox gradient engine



rowvec CoxenggrC( const vec& beta, const mat& data){

 mat sortedd = data.rows(sort_index(data.col(0))) ;
 
  mat x = sortedd.cols(2, sortedd.n_cols-1);

  vec a = exp( x*beta  );

  mat A = risksetC( x.each_col()%a );
  vec B = risksetC( a );

  mat dA = A.each_col()/B;
 
  mat C = dA.rows( find( sortedd.col(1)==1 ) );        
  
    colvec ones(C.n_rows, fill::ones);  

    return ones.t()*C ;

   }



// [[Rcpp::export]]

rowvec grCox(const vec& beta, const vec& y, const mat& x , const vec& str){

 vec str_levs = unique(str);

   int nstr = str_levs.n_elem;

   mat stratasets(nstr, beta.n_elem);

   for(int sl=0; sl<nstr; ++sl){

     mat Xsl = subdata_arma( x, str, str_levs(sl));
  
       stratasets.row(sl) =  CoxenggrC( beta, Xsl ) ;  

   }

     colvec ones(stratasets.n_rows, fill::ones);  
    
     rowvec den =  ones.t()*stratasets;

     return y.t() - den ;
   

    } 



// Cox hessian engine


rowvec as_vec(const mat& x){

  vec res = vectorise(x);
  return res.t();

   }

// 

mat CoxenghessC(const vec& beta, const mat& data){

 mat sortedd = data.rows(sort_index(data.col(0))) ;
 
  mat x = sortedd.cols(2, sortedd.n_cols-1);

  vec a = exp( x*beta  );

  mat Ax = risksetC( x.each_col()%a );
  vec ra = risksetC( a );

  mat Axi(x.n_rows, beta.n_elem*beta.n_elem );

    for(int i=0; i<x.n_rows; ++i ){
     
      double c = as_scalar( exp( dot(beta, x.row(i) ) )  );
      mat z =  c*x.row(i).t()*x.row(i);
      Axi.row(i) = as_vec(z);  // u could in principle speed-up by avoiding vectorization, but you should code a function that takes the elemnt-wise cumulative sum of a matrix first.

     }


    mat dAcs = risksetC( Axi);
   
    mat dAfull = dAcs.each_col()/ra;  // for each cube on a matrix level ??
      
     mat dA = dAfull.rows( find( sortedd.col(1)==1 ) );        

         colvec ones(dA.n_rows, fill::ones);   

     	rowvec Avec = ones.t()*dA;

	mat A = mat(Avec.begin(), beta.n_elem, beta.n_elem);

	mat dBfull = Ax.each_col()/pow(ra,2);
        
      mat dB = dBfull.rows( find( sortedd.col(1)==1 ) );        
      
      mat B = dB.t()*Ax.rows( find( sortedd.col(1)==1 ) );          

      mat den = A - B;
   
      return den;


   }



// [[Rcpp::export]]


mat hessCox(const vec& beta, const vec& y, const mat& x , const vec& str){


 vec str_levs = unique(str);

   int nstr = str_levs.n_elem;

   mat stratasets(nstr, beta.n_elem*beta.n_elem);

   for(int sl=0; sl<nstr; ++sl){

     mat Xsl = subdata_arma( x, str, str_levs(sl));
      
     mat z = CoxenghessC( beta, Xsl );
     stratasets.row(sl) =  as_vec(z) ;  

   }


     colvec ones(stratasets.n_rows, fill::ones);  
    
     rowvec den =  ones.t()*stratasets;

   mat A = mat(den.begin(), beta.n_elem, beta.n_elem);
   
   return (-1)*A;


      }











// // Breslow estimator :only return the cumulative hazard images (no time axis)


// mat BreslowNelsonAalen_estimation(const vec& betahat, const mat& data, const vec& znull){

//   mat sortedd = data.rows(sort_index(data.col(0))) ;
 
//  mat  event_matrix =  sortedd.rows( find( sortedd.col(1) == 1 ) );

//  vec event_times = event_matrix.col(0);
 
//   mat x = sortedd.cols(2, sortedd.n_cols-1);

//   vec a = exp( x*betahat );
 
//   vec risksets = flipud( cumsum( flipud( a ) ) ) ; 

//   vec NelsonAalen_typeHazard =  cumsum( 1 / risksets.rows( find( sortedd.col(1) == 1 ) ) ); 

//   double num = as_scalar( exp( betahat.t()*znull ) );  // Breslow correction
  
//   vec Breslow_typeHazard = num*NelsonAalen_typeHazard ;

//   mat out( event_times.n_elem, 2 );

//   out.col(0) = event_times;

//   out.col(1) =  Breslow_typeHazard;
   
//   return out ;

//     }








