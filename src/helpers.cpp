template <typename T>
T sortByDimNames(const T m){
  List dimNames = m.attr("dimnames");
  CharacterVector colNames = dimNames[1];
  CharacterVector rowNames = dimNames[0];
  int sizeMatr = colNames.size();
  CharacterVector sortedNames(sizeMatr);
  for(int i = 0; i < sizeMatr; i++)
    sortedNames(i) = rowNames(i);
  sortedNames.sort();
  
  NumericVector colIdx(sizeMatr), rowIdx(sizeMatr);
  for(int i = 0; i < sizeMatr; i++){
    for(int j = 0; j < sizeMatr; j++){
      if(colNames(j) == sortedNames(i))
        colIdx(i) = j;
      if(rowNames(j) == sortedNames(i))
        rowIdx(i) = j;
    }
  }
  
  T sortedM(sizeMatr);
  sortedM.attr("dimnames") = List::create(sortedNames, sortedNames);
  
  for(int i = 0; i < sizeMatr; i++)
    for(int j = 0; j < sizeMatr; j++)
      sortedM(i, j) = m(rowIdx(i), colIdx(j));
      
  return sortedM;
}

double lbeta(double p, double q){
  return lgamma(p) + lgamma(q) - lgamma(p + q);
}

template <typename T>
T _transpose(T & mat) {      
  int numRows = mat.nrow();
  int numCols = mat.ncol();
  
  T transpose(numCols, numRows);
  // Assign dim names transposed (rows to cols and viceversa)
  transpose.attr("dimnames") = List::create(colnames(mat), rownames(mat));

  for (int i = 0; i < numCols; ++i) {
    transpose(i, _) = mat(_, i);
  }
  
  return transpose;
}

//****************************************************************************80

double betain ( double x, double p, double q, double beta )

//****************************************************************************80
//
//  Purpose:
//
//    BETAIN computes the incomplete Beta function ratio.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    25 September 2014
//
//  Author:
//
//    Original FORTRAN77 version by KL Majumder, GP Bhattacharjee.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    KL Majumder, GP Bhattacharjee,
//    Algorithm AS 63:
//    The incomplete Beta Integral,
//    Applied Statistics,
//    Volume 22, Number 3, 1973, pages 409-411.
//
//  Parameters:
//
//    Input, double X, the argument, between 0 and 1.
//
//    Input, double P, Q, the parameters, which
//    must be positive.
//
//    Input, double BETA, the logarithm of the complete
//    beta function.
//
//    Output, int &IFAULT, error flag.
//    0, no error.
//    nonzero, an error occurred.
//
//    Output, double BETAIN, the value of the incomplete
//    Beta function ratio.
//
{
  double acu = 0.1E-14;
  double ai;
  double cx;
  bool indx;
  int ns;
  double pp;
  double psq;
  double qq;
  double rx;
  double temp;
  double term;
  double value;
  double xx;

  value = x;
//
//  Special cases.
//
  if ( x == 0.0 || x == 1.0 )
  {
    return value;
  }
//
//  Change tail if necessary and determine S.
//
  psq = p + q;
  cx = 1.0 - x;

  if ( p < psq * x )
  {
    xx = cx;
    cx = x;
    pp = q;
    qq = p;
    indx = true;
  }
  else
  {
    xx = x;
    pp = p;
    qq = q;
    indx = false;
  }

  term = 1.0;
  ai = 1.0;
  value = 1.0;
  ns = ( int ) ( qq + cx * psq );
//
//  Use the Soper reduction formula.
//
  rx = xx / cx;
  temp = qq - ai;
  if ( ns == 0 )
  {
    rx = xx;
  }

  for ( ; ; )
  {
    term = term * temp * rx / ( pp + ai );
    value = value + term;;
    temp = fabs ( term );

    if ( temp <= acu && temp <= acu * value )
    {
      value = value * exp ( pp * log ( xx ) 
      + ( qq - 1.0 ) * log ( cx ) - beta ) / pp;

      if ( indx )
      {
        value = 1.0 - value;
      }
      break;
    }

    ai = ai + 1.0;
    ns = ns - 1;

    if ( 0 <= ns )
    {
      temp = qq - ai;
      if ( ns == 0 )
      {
        rx = xx;
      }
    }
    else
    {
      temp = psq;
      psq = psq + 1.0;
    }
  }

  return value;
}
//****************************************************************************80

double r8_max ( double x, double y )

//****************************************************************************80
//
//  Purpose:
//
//    R8_MAX returns the maximum of two R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, Y, the quantities to compare.
//
//    Output, double R8_MAX, the maximum of X and Y.
//
{
  double value;

  if ( y < x )
  {
    value = x;
  } 
  else
  {
    value = y;
  }
  return value;
}

//****************************************************************************80

double xinbta ( double p, double q, double beta, double alpha )

//****************************************************************************80
//
//  Purpose:
//
//    XINBTA computes inverse of the incomplete Beta function.
//
//  Discussion:
//
//    The accuracy exponent SAE was loosened from -37 to -30, because
//    the code would not otherwise accept the results of an iteration
//    with p = 0.3, q = 3.0, alpha = 0.2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    25 September 2014
//
//  Author:
//
//    Original FORTRAN77 version by GW Cran, KJ Martin, GE Thomas.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    GW Cran, KJ Martin, GE Thomas,
//    Remark AS R19 and Algorithm AS 109:
//    A Remark on Algorithms AS 63: The Incomplete Beta Integral
//    and AS 64: Inverse of the Incomplete Beta Integeral,
//    Applied Statistics,
//    Volume 26, Number 1, 1977, pages 111-114.
//
//  Parameters:
//
//    Input, double P, Q, the parameters of the incomplete
//    Beta function.
//
//    Input, double BETA, the logarithm of the value of
//    the complete Beta function.
//
//    Input, double ALPHA, the value of the incomplete Beta
//    function.  0 <= ALPHA <= 1.
//
//    Output, int &IFAULT, error flag.
//    0, no error occurred.
//    nonzero, an error occurred.
//
//    Output, double XINBTA, the argument of the incomplete
//    Beta function which produces the value ALPHA.
//
//  Local Parameters:
//
//    Local, double SAE, requests an accuracy of about 10^SAE.
//
{
  double a;
  double acu;
  double adj;
  double fpu;
  double g;
  double h;
  int iex;
  bool indx;
  double pp;
  double prev;
  double qq;
  double r;
  double s;
  double sae = -30.0;
  double sq;
  double t;
  double tx;
  double value;
  double w;
  double xin;
  double y;
  double yprev;

  fpu = pow ( 10.0, sae );

  value = alpha;

//
//  If the answer is easy to determine, return immediately.
//
  if ( alpha == 0.0 )
  {
    value = 0.0;
    return value;
  }

  if ( alpha == 1.0 )
  {
    value = 1.0;
    return value;
  }
//
//  Change tail if necessary.
//
  if ( 0.5 < alpha )
  {
    a = 1.0 - alpha;
    pp = q;
    qq = p;
    indx = true;
  }
  else
  {
    a = alpha;
    pp = p;
    qq = q;
    indx = false;
  }
//
//  Calculate the initial approximation.
//
  r = sqrt ( - log ( a * a ) );

  y = r - ( 2.30753 + 0.27061 * r ) 
    / ( 1.0 + ( 0.99229 + 0.04481 * r ) * r );

  if ( 1.0 < pp && 1.0 < qq )
  {
    r = ( y * y - 3.0 ) / 6.0;
    s = 1.0 / ( pp + pp - 1.0 );
    t = 1.0 / ( qq + qq - 1.0 );
    h = 2.0 / ( s + t );
    w = y * sqrt ( h + r ) / h - ( t - s ) 
      * ( r + 5.0 / 6.0 - 2.0 / ( 3.0 * h ) );
    value = pp / ( pp + qq * exp ( w + w ) );
  }
  else
  {
    r = qq + qq;
    t = 1.0 / ( 9.0 * qq );
    t = r * pow ( 1.0 - t + y * sqrt ( t ), 3 );

    if ( t <= 0.0 )
    {
      value = 1.0 - exp ( ( log ( ( 1.0 - a ) * qq ) + beta ) / qq );
    }
    else
    {
      t = ( 4.0 * pp + r - 2.0 ) / t;

      if ( t <= 1.0 )
      {
        value = exp ( ( log ( a * pp ) + beta ) / pp );
      }
      else
      {
        value = 1.0 - 2.0 / ( t + 1.0 );
      }
    }
  }
//
//  Solve for X by a modified Newton-Raphson method,
//  using the function BETAIN.
//
  r = 1.0 - pp;
  t = 1.0 - qq;
  yprev = 0.0;
  sq = 1.0;
  prev = 1.0;

  if ( value < 0.0001 )
  {
    value = 0.0001;
  }

  if ( 0.9999 < value )
  {
    value = 0.9999;
  }

  iex = r8_max ( - 5.0 / pp / pp - 1.0 / pow ( a, 0.2 ) - 13.0, sae );

  acu = pow ( 10.0, iex );
//
//  Iteration loop.
//
  for ( ; ; )
  {
    y = betain ( value, pp, qq, beta );

    xin = value;
    y = ( y - a ) * exp ( beta + r * log ( xin ) + t * log ( 1.0 - xin ) );

    if ( y * yprev <= 0.0 )
    {
      prev = r8_max ( sq, fpu );
    }

    g = 1.0;

    for ( ; ; )
    {
      for ( ; ; )
      {
        adj = g * y;
        sq = adj * adj;

        if ( sq < prev )
        {
          tx = value - adj;

          if ( 0.0 <= tx && tx <= 1.0 )
          {
            break;
          }
        }
        g = g / 3.0;
      }
//
//  Check whether the current estimate is acceptable.
//  The change "VALUE = TX" was suggested by Ivan Ukhov.
//
      if ( prev <= acu || y * y <= acu )
      {
        value = tx;
        if ( indx )
        {
          value = 1.0 - value;
        }
        return value;
      }

      if ( tx != 0.0 && tx != 1.0 )
      {
        break;
      }

      g = g / 3.0;
    }

    if ( tx == value )
    {
      break;
    }

    value = tx;
    yprev = y;
  }

  if ( indx )
  {
    value = 1.0 - value;
  }

  return value;
}
