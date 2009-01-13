/**
 * @file GaussianQuadrature.h
 * @brief Wrapper for dgaus8.
 * @author J. Chiang
 *
 * $Header$
 */

#ifndef genericSources_GaussianQuadrature_h
#define genericSources_GaussianQuadrature_h

#include <cmath>
#include <stdexcept>

namespace {
   inline const double sign(const double & x, const double & y) {
      return y >= 0. ? fabs(x) : -fabs(x);
   }
}

namespace genericSources {

class GaussianQuadrature {

public:

   class dgaus8Exception : public std::runtime_error {
   public:
      dgaus8Exception(const std::string & message, int ierr) 
         : std::runtime_error(message), m_ierr(ierr) {}
      virtual ~dgaus8Exception() throw() {}
      int errCode() const {
         return m_ierr;
      }
   private:
      int m_ierr;
   };

/**
 * @brief This version of dgaus8 has been translated from the Fortran
 * version to Java and thence to C++. 
 * See http://www1.fpl.fs.fed.us/Gaus8.np.java
 *
 * @author P. Nolan
 * @author J. Chiang
 *
 */
   template<typename Functor>
   static double dgaus8(Functor & fun, double a, double b,
                        double & err, int & ierr) {
      const double x1 = 1.83434642495649805E-01;     
      const double x2 = 5.25532409916328986E-01;
      const double x3 = 7.96666477413626740E-01;     
      const double x4 = 9.60289856497536232E-01;
  
      const double w1 = 3.62683783378361983E-01;     
      const double w2 = 3.13706645877887287E-01;
      const double w3 = 2.22381034453374471E-01;     
      const double w4 = 1.01228536290376259E-01;
  
      const double sq2 = 1.41421356E0;
  
      const int nlmn = 1;
      const int kmx = 5000; 
      const int kml = 6;
  
      int k,l,lmn,lmx,mxl,nbits,nib,nlmx;

      double ae,anib,area,c,ce,ee,ef,eps,est,
         gl,glr,tol,vr,x,h,g8xh,ans;

      double aa[61], hh[61], vl[61], gr[61];
      
      int lr[61];
  
      k = 53;

//  r1mach(5) is the log base 10 of 2 = .3010299957

//      anib = .3010299957*k/0.30102000E0;

//      nbits = (int)(anib);

      nbits = 53;
      nlmx = std::min(60, (nbits*5)/8);
      ans = 0.0;
      ierr = 1;
      ce = 0.0;

      if (a == b) {

// 140
         if (err < 0.0) err = ce;
         return ans;
      }
  
      lmx = nlmx;
      lmn = nlmn;
  
      if (b != 0.0) {
    
         if (sign(1.0,b)*a > 0.0) {
      
            c = fabs(1.0 - a/b);
      
            if (c <= 0.1) {
	
               if (c <= 0.0) {
	  
                  // 140
	  
                  if (err < 0.0) err = ce;
	  
                  return ans;
	  
               }
	
               anib = 0.5 - log(c)/0.69314718E0;
               nib = (int)(anib);
               lmx = std::min(nlmx, nbits-nib-7);
               if (lmx < 1) {
	    
                  // 130
                  ierr = -1;

//                   throw std::runtime_error("dgaus8 --- a and b are too "
//                                            "nearly equal to allow normal "
//                                            "integration. "
//                                            "ans is set to 0 and ierr "
//                                            "is set to -1.");
                  std::string message("dgaus8 --- a and b are too "
                                      "nearly equal to allow normal "
                                      "integration. "
                                      "ans is set to 0 and ierr "
                                      "is set to -1.");
                  throw dgaus8Exception(message, ierr);
                  if (err < 0.0) err = ce;
                  return ans;
               }
               lmn = std::min(lmn, lmx);
            }
      
         }
    
      }
  
      // 10
  
      tol = std::max(fabs(err),pow(2.0,5-nbits))/2.0;
  
      if (err == 0.0) {
         // According to the SLATEC documentation 2.22e-16 is the
         // appropriate value for IEEE 754 double precision unit
         // roundoff.
         tol = sqrt(2.22e-16);
      }

      eps = tol;

      hh[1] = (b - a)/4.0;
      aa[1] = a;
      lr[1] = 1;

      l = 1;
  
      x = aa[l] + 2.0*hh[l];
      h = 2.0*hh[l];
  
      g8xh = h*((w1*(fun(x - x1*h) + fun(x + x1*h)) +
                 w2*(fun(x - x2*h) + fun(x + x2*h))) +
                (w3*(fun(x - x3*h) + fun(x + x3*h)) +
                 w4*(fun(x - x4*h) + fun(x + x4*h))));
	    
      est = g8xh;
      k = 8;
      area = fabs(est);
      ef = 0.5;
      mxl = 0;
  
//
//    Compute refined estimates, estimate the error, etc.
//

//   20 loop

      while (true) {
    
         x = aa[l] + hh[l];
         h = hh[l];
    
         g8xh = h*((w1*(fun(x - x1*h) + fun(x + x1*h)) +
                    w2*(fun(x - x2*h) + fun(x + x2*h))) +
                   (w3*(fun(x - x3*h) + fun(x + x3*h)) +
                    w4*(fun(x - x4*h) + fun(x + x4*h))));
    
         gl = g8xh;         
    
         x = aa[l] + 3.0*hh[l];
         h = hh[l];

         g8xh = h*((w1*(fun(x - x1*h) + fun(x + x1*h)) +
                    w2*(fun(x - x2*h) + fun(x + x2*h))) +
                   (w3*(fun(x - x3*h) + fun(x + x3*h)) +
                    w4*(fun(x - x4*h) + fun(x + x4*h))));

         gr[l] = g8xh;         
         k += 16;
         area += (fabs(gl) + fabs(gr[l]) - fabs(est));
         glr = gl + gr[l];
         ee = fabs(est - glr)*ef;
         ae = std::max(eps*area,tol*fabs(glr));
         if (ee - ae > 0.0) {

// 50:

            if (k > kmx) lmx = kml;
            if (l < lmx) {
               l++;
               eps *= 0.5;
               ef /= sq2;
               hh[l] = hh[l-1]*0.5;
               lr[l] = -1;
               aa[l] = aa[l-1];
               est = gl;
            } else {
// 30:
               mxl = 1;
// 40,2:
               ce += (est - glr);
               if (lr[l] <= 0) {
// 60,1:
                  vl[l] = glr;
// 70,1:
                  est = gr[l-1];
                  lr[l] = 1;
                  aa[l] += 4.0*hh[l];
               } else {
// 80,1:
                  vr = glr;
// 90,1:
                  while (true) {
                     if (l <= 1) {
// 120,1:
                        ans = vr;
                        if ((mxl != 0) && (fabs(ce) > 2.0*tol*area)) {
                           ierr = 2;
//                            throw std::runtime_error("ans is probably "
//                                                     "insufficiently accurate");
                           std::string message("ans is probably "
                                               "insufficiently accurate");
                           throw dgaus8Exception(message, ierr);
                        }
// 140,1:
                        if (err < 0.0) err = ce;
                        return ans;
                     }
                     l--;
                     eps *= 2.0;
                     ef *= sq2;
                     if (lr[l] <= 0.0) break;
// 110,1:
                     vr += vl[l+1];
// end of 90,1 loop:
                  }
// 100,1:
                  vl[l] = vl[l+1] + vr;
// 70,3:
                  est = gr[l-1];
                  lr[l] = 1;
                  aa[l] += 4.0*hh[l];
// end of 80,1:
               }
	       // end of 30 branch:
            }

// end of 50 branch:
	    
         } else {

// 40,1:
            ce += (est - glr);
            if (lr[l] <= 0) {
// 60,2:
               vl[l] = glr;
// 70,2:
               est = gr[l-1];
               lr[l] = 1;
               aa[l] += 4.0*hh[l];
            } else {
// 80,2:
               vr = glr;
// 90,2:
               while (true) {
                  if (l <= 1) {
// 120,2:
                     ans = vr;
                     if ((mxl != 0) && (fabs(ce) > 2.0*tol*area)) {
                        ierr = 2;
//                         throw std::runtime_error("ans is probably "
//                                                  "insufficiently accurate.");
                        std::string message("ans is probably "
                                            "insufficiently accurate");
                        throw dgaus8Exception(message, ierr);
                     }
// 140,2:
                     if (err < 0.0) err = ce;
                     return ans;
                  }
                  l--;
                  eps *= 2.0;
                  ef *= sq2;
                  if (lr[l] <= 0.0) break;
// 110,2:
                  vr += vl[l+1];
// end of 90,2 loop:
               }
// 100,2:
               vl[l] = vl[l+1] + vr;
// 70,4:
               est = gr[l-1];
               lr[l] = 1;
               aa[l] += 4.0*hh[l];
// end of 80,2:
            }
// end of 40,1:
         }
// end of 20 loop:
      }
   }

};

} // namespace genericSources

#endif // genericSources_GaussianQuadrature_h
