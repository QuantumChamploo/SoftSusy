void MssmSoftPars::FGMCaseB1(const MssmSusy & xx, double LAMBDA, 
			       double mMess, double epsu, double epsd, double epse, double cgrav) {
  
  // Modified thresholds by JEL 1-26-04 to accomodate numerical infinities
  
  const double epstol = 1.0e-4;
  double x = LAMBDA / mMess;
  double lambda = LAMBDA/(16.0 * sqr(PI));
  
  double f, g;
  
  if(fabs(x) < epstol) { /// hep-ph/9801271
    g = 1.0 + x*x/6.0 + sqr(x*x)/15.0;
    f = 1.0 + x*x/36.0 - 11.0*sqr(x*x)/450.0;
  }
  else if(fabs(x-1.0) < 0.0001) {
    g  =  log(4.0);
    f  = -sqr(PI)/6.0 + log(4.0) + 0.5*sqr(log(4.0));
    g -=  0.0008132638905771205626;
    f -= -0.0049563838821509165200;
  }
  else {
    g = 1.0 / sqr(x) * 
      ((1.0 + x) * log(1.0 + x) + (1.0 - x) * log(1.0 - x));
    f = (1.0 + x) / sqr(x) * 
      (log(1.0 + x) - 2.0 * dilog(x / (1.0 + x)) + 0.5 * 
       dilog(2.0 * x / (1.0 + x))) + 
      (1.0 - x) / sqr(x) * (log(1.0 - x) - 2.0 * dilog(-x / (1.0 - x)) +
			    0.5 * dilog(-2.0 * x / (1.0 - x)));
  }
  
  double n5d = 2;
  
  double g12 = sqr(xx.displayGaugeCoupling(1));
  double g22 = sqr(xx.displayGaugeCoupling(2));
  double g32 = sqr(xx.displayGaugeCoupling(3));
  
  double yt = xx.displayYukawaElement(YU,3,3); 
  double yb = xx.displayYukawaElement(YD,3,3); 
  double ytau = xx.displayYukawaElement(YE,3,3);
  
  double yt2 = sqr(yt);
  double yb2 = sqr(yb);
  double ytau2 = sqr(ytau);
  
  /// There is a relative minus in the mGMSB conditions for gaugino masses,
  /// since these equations are for L=-M/2 gaugino gaugino. See hep-ph/9801271:
  /// BCA 27/7/12
  double m1, m2, m3;
  m1 = n5d * g12 * lambda * g; 
  m2 = n5d * g22 * lambda * g; 
  m3 = n5d * g32 * lambda * g; 
  setGauginoMass(1, m1);   setGauginoMass(2, m2);   setGauginoMass(3, m3);
  
  setM32(2.37e-19 * LAMBDA * mMess * cgrav);
  
  double g1f = sqr(g12);
  double g2f = sqr(g22);
  double g3f = sqr(g32);
  
  double mursq, mdrsq, mersq, mqlsq, mllsq;
  mursq = 2.0 * f * sqr(lambda) * n5d * 
    (4.0 / 3.0 * g3f + 0.6 * 4.0 / 9.0 * g1f);
  mdrsq = 2.0 * f * sqr(lambda) * n5d * 
    (4.0 / 3.0 * g3f + 0.6 * 1.0 / 9.0 * g1f);
  mersq = 2.0 * f * sqr(lambda) * n5d * 
    (0.6 * g1f);
  mqlsq = 2.0 * f * sqr(lambda) * n5d * 
    (4.0 / 3.0 * g3f + 0.75 * g2f + 0.6 * g1f / 36.0);
  mllsq = 2.0 * f * sqr(lambda) * n5d * 
    (0.75 * g2f + 0.6 * 0.25 * g1f) ;
  
  double a0 = 0.,a1 = 0.,a2 = 0.,a3 = 0.,a4 = 0.,a5 = 0.,a6 = 0.,a7 = 0.,a8 = 0.;
  
  
  //Ql
a0= f * sqr(lambda) * (6*sqr(yt2) + 6*sqr(yb2) + 2*yb2*yt2 + yb2*ytau2 - 16./3.*g32*yt2 - 3*g22*yt2 - 13./15.*g12*yt2 - 16./3.*g32*yb2 - 3*g22*yb2 - 7./15.*g12*yb2 + epsu * ( 4./9.*16./3.*g32*yt2 + 4./9.*3.*g22*yt2 + 4./9.*13./15.*g12*yt2 - 4./9.*yb2*yt2 - 32./9.*sqr(yt2) ) + epsd * ( 4./9.*16./3.*g32*yb2 + 4./9.*3*g22*yb2 + 4./9.*7./15.*g12*yb2 - 4./9.*yb2*yt2 - 4./9.*yb2*ytau2 - 32./9.*sqr(yb2) ));
a1= f * sqr(lambda) * (0);
a2= f * sqr(lambda) * (0);
a3= f * sqr(lambda) * (0);
a4= f * sqr(lambda) * (6*sqr(yt2) + 6*sqr(yb2) - 16./3.*g32*yt2 - 3*g22*yt2 - 13./15.*g12*yt2 - 16./3.*g32*yb2 - 3*g22*yb2 - 7./15.*g12*yb2 + 2*yb2*yt2 + epsu*( 10./3.*sqr(yt2) - 4./9.*16./3.*g32*yt2 - 4./9.*3.*g22*yt2 - 4./9.*13./15.*g12*yt2 + 4./9.*yb2*yt2 ) + epsd *(32./9.*sqr(yb2) - 4./9.*16./3.*g32*yb2 - 4./9.*3*g22*yb2 - 4./9.*7./15.*g12*yb2 + 4./9.*yb2*ytau2 + 4./9.*yt2*yb2 ) );
a5= f * sqr(lambda) * ( epsu*(-0.628539361*16./3.*g32*yt2 - 0.628539361*3*g22*yt2 - 0.628539361*13./15.*g12*yt2 + 0.628539361*yb2*yt2 + 3.771236166*sqr(yt2)) + epsd*(-0.628539361*16./3.*g32*yb2 - 0.628539361*3*g22*yb2 - 0.628539361*7./15.*g12*yb2 + 3.771236166*sqr(yb2) + 0.628539361*yb2*ytau2 + 0.628539361*yb2*yt2 ) );
a6= f * sqr(lambda) * (0);
a7= f * sqr(lambda) * ( epsu*(-0.628539361*16./3.*g32*yt2 - 0.628539361*3*g22*yt2 - 0.628539361*13./15.*g12*yt2 + 0.628539361*yb2*yt2 + 3.771236166*sqr(yt2)) + epsd*(-0.628539361*16./3.*g32*yb2 - 0.628539361*3*g22*yb2 - 0.628539361*7./15.*g12*yb2 + 3.771236166*sqr(yb2) + 0.628539361*yb2*ytau2 + 0.628539361*yb2*yt2 ) );
a8= f * sqr(lambda) * (0);
 
  setSoftMassElement(mQl,1,1, mqlsq  + a0);
  setSoftMassElement(mQl,1,2, a1);
  setSoftMassElement(mQl,1,3, a2);
  setSoftMassElement(mQl,2,1, a3);
  setSoftMassElement(mQl,2,2, mqlsq  + a4);
  setSoftMassElement(mQl,2,3, a5);
  setSoftMassElement(mQl,3,1, a6);
  setSoftMassElement(mQl,3,2, a7);
  setSoftMassElement(mQl,3,3, mqlsq  + a8);
  
  //Ur
a0= f * sqr(lambda) * (12*sqr(yt2) + 2*yb2*yt2 - 2.*16./3.*g32*yt2 - 2*3*g22*yt2 - 2.*13./15.*g12*yt2 + epsu*(8./9.*16./3.*g32*yt2 + 8./9.*3.*g22*yt2 + 8./9.*13./15.*g12*yt2 - 8./9.*yb2*yt2 - 64./9.*sqr(yt2) ) );
a1= f * sqr(lambda) * (0);
a2= f * sqr(lambda) * (0);
a3= f * sqr(lambda) * (0);
a4= f * sqr(lambda) * (12*sqr(yt2) + 2*yb2*yt2 - 2.*16./3.*g32*yt2 - 2*3*g22*yt2 - 2.*13./15.*g12*yt2 + epsu*(-8./9.*16./3.*g32*yt2 - 8./9.*3.*g22*yt2 - 8./9.*13./15.*g12*yt2 + 8./9.*yb2*yt2 + 64./9.*sqr(yt2) ) );
a5= f * sqr(lambda) * (epsu*(-1.257078722*16./3.*g32*yt2 - 1.257078722*3*g22*yt2 - 1.257078722*13./15.*g12*yt2 + 1.257078722*yb2*yt2 + 7.542472333*sqr(yt2)));
a6= f * sqr(lambda) * (0);
a7= f * sqr(lambda) * (epsu*(-1.257078722*16./3.*g32*yt2 - 1.257078722*3*g22*yt2 - 1.257078722*13./15.*g12*yt2 + 1.257078722*yb2*yt2 + 7.542472333*sqr(yt2)));
a8= f * sqr(lambda) * (0);
  
  setSoftMassElement(mUr,1,1, mursq  + a0);
  setSoftMassElement(mUr,1,2, a1);
  setSoftMassElement(mUr,1,3, a2);
  setSoftMassElement(mUr,2,1, a3);
  setSoftMassElement(mUr,2,2, mursq  + a4);
  setSoftMassElement(mUr,2,3, a5);
  setSoftMassElement(mUr,3,1, a6);
  setSoftMassElement(mUr,3,2, a7);
  setSoftMassElement(mUr,3,3, mursq  + a8);
  
  //Dr 
a0= f * sqr(lambda) * (-2.*16./3.*g32*yb2 - 2*3*g22*yb2 - 2.*7./15.*g12*yb2 + 12*sqr(yb2) + 2*yb2*yt2 + 2*yb2*ytau2 + epsd*(8./9.*16./3.*g32*yb2 + 8./9.*3.*g22*yb2 + 8./9.*7./15.*g12*yb2 - 64./9.*sqr(yb2) - 8./9.*yb2*yt2 - 8./9.*yb2*ytau2) );
a1= f * sqr(lambda) * (0);
a2= f * sqr(lambda) * (0);
a3= f * sqr(lambda) * (0);
a4= f * sqr(lambda) * (-2.*16./3.*g32*yb2 - 2*3*g22*yb2 - 2.*7./15.*g12*yb2 + 12*sqr(yb2) + 2*yb2*yt2 + 2*yb2*ytau2 + epsd*(-8./9.*16./3.*g32*yb2 - 8./9.*3.*g22*yb2 - 8./9.*7./15.*g12*yb2 + 64./9.*sqr(yb2) + 8./9.*yb2*yt2 + 8./9.*yb2*ytau2) );
a5= f * sqr(lambda) * (epsd*(-1.257078722*16./3.*g32*yb2 - 1.257078722*3.*g22*yb2 - 1.257078722*7./15.*g12*yb2 + 7.542472333*sqr(yb2) + 1.257078722*yb2*yt2 + 1.257078722*yb2*ytau2 ));
a6= f * sqr(lambda) * (0);
a7= f * sqr(lambda) * (epsd*(-1.257078722*16./3.*g32*yb2 - 1.257078722*3.*g22*yb2 - 1.257078722*7./15.*g12*yb2 + 7.542472333*sqr(yb2) + 1.257078722*yb2*yt2 + 1.257078722*yb2*ytau2 ));
a8= f * sqr(lambda) * (0);
  
  setSoftMassElement(mDr,1,1, mdrsq  + a0);
  setSoftMassElement(mDr,1,2, a1);
  setSoftMassElement(mDr,1,3, a2);
  setSoftMassElement(mDr,2,1, a3);
  setSoftMassElement(mDr,2,2, mdrsq  + a4);
  setSoftMassElement(mDr,2,3, a5);
  setSoftMassElement(mDr,3,1, a6);
  setSoftMassElement(mDr,3,2, a7);
  setSoftMassElement(mDr,3,3, mdrsq  + a8);
  
  //Ll
a0= f * sqr(lambda) * (-3*g22*ytau2 - 9./5.*g12*ytau2 + 3*yb2*ytau2 + 4*sqr(ytau2) + epse*(4./9.*3.*g22*ytau2 + 4./9.*9./5.*g12*ytau2 - 4./3.*yb2*ytau2 - 8./3.*sqr(ytau2)));
a1= f * sqr(lambda) * (0);
a2= f * sqr(lambda) * (0);
a3= f * sqr(lambda) * (0);
a4= f * sqr(lambda) * (-3*g22*ytau2 - 9./5.*g12*ytau2 + 3*yb2*ytau2 + 4*sqr(ytau2) + epse*(-4./9.*3.*g22*ytau2 - 4./9.*9./5.*g12*ytau2 + 4./3.*yb2*ytau2 + 8./3.*sqr(ytau2)));
a5= f * sqr(lambda) * (epse*(-0.628539361*3*g22*ytau2 - 0.628539361*9./5.*g12*ytau2 + 1.885618083*yb2*ytau2 + 2.514157444*sqr(ytau2) ));
a6= f * sqr(lambda) * (0);
a7= f * sqr(lambda) * (epse*(-0.628539361*3*g22*ytau2 - 0.628539361*9./5.*g12*ytau2 + 1.885618083*yb2*ytau2 + 2.514157444*sqr(ytau2) ));
a8= f * sqr(lambda) * (0);
  
  setSoftMassElement(mLl,1,1, mllsq  + a0);
  setSoftMassElement(mLl,1,2, a1);
  setSoftMassElement(mLl,1,3, a2);
  setSoftMassElement(mLl,2,1, a3);
  setSoftMassElement(mLl,2,2, mllsq  + a4);
  setSoftMassElement(mLl,2,3, a5);
  setSoftMassElement(mLl,3,1, a6);
  setSoftMassElement(mLl,3,2, a7);
  setSoftMassElement(mLl,3,3, mllsq  + a8);
  
  
  //Er
a0= f * sqr(lambda) * (-2.*3.*g22*ytau2 - 2.*9./5.*g12*ytau2 + 6*yb2*ytau2 + 8*sqr(ytau2) + epse*(8./9.*3.*g22*ytau2 + 8./9.*9./5.*g12*ytau2 - 8./3.*yb2*ytau2 - 16./3.*sqr(ytau2) ) );
a1= f * sqr(lambda) * (0);
a2= f * sqr(lambda) * (0);
a3= f * sqr(lambda) * (0);
a4= f * sqr(lambda) * (-2.*3.*g22*ytau2 - 2.*9./5.*g12*ytau2 + 6*yb2*ytau2 + 8*sqr(ytau2) + epse*(-8./9.*3.*g22*ytau2 - 8./9.*9./5.*g12*ytau2 + 8./3.*yb2*ytau2 + 16./3.*sqr(ytau2) ) );
a5= f * sqr(lambda) * (epse*(-1.257078722*3.*g22*ytau2 - 1.257078722*9./5.*g12*ytau2 + 3.771236166*yb2*ytau2 + 5.028314188*sqr(ytau2) ));
a6= f * sqr(lambda) * (0);
a7= f * sqr(lambda) * (epse*(-1.257078722*3.*g22*ytau2 - 1.257078722*9./5.*g12*ytau2 + 3.771236166*yb2*ytau2 + 5.028314188*sqr(ytau2) ));
a8= f * sqr(lambda) * (0);
  
  setSoftMassElement(mEr,1,1, mersq + a0);
  setSoftMassElement(mEr,1,2, a1);
  setSoftMassElement(mEr,1,3, a2);
  setSoftMassElement(mEr,2,1, a3);
  setSoftMassElement(mEr,2,2, mersq + a4);
  setSoftMassElement(mEr,2,3, a5);
  setSoftMassElement(mEr,3,1, a6);
  setSoftMassElement(mEr,3,2, a7);
  setSoftMassElement(mEr,3,3, mersq + a8);
  
  //higgs masses
  double deltaHu = 0.;
  double deltaHd = 0.;
  
  setMh1Squared(mllsq + deltaHd); //Hd
  setMh2Squared(mllsq + deltaHu); //Hu
  
  //trilinears
  
  setTrilinearElement(UA, 2, 2, -lambda * epsu * (-2./9.*yb2*yt - 2./3.*yt*yt2));
  setTrilinearElement(UA, 2, 3, -lambda * (epsu * (0.628539361*yt*yt2) + epsd * (0.628539361*yb2*yt)) );
  setTrilinearElement(UA, 3, 2, -lambda * epsu * 1.257078722 * yt * yt2);
  
  setTrilinearElement(DA, 2, 2, -lambda * epsd * (-2./3.*yb*yb2 - 2./9.*yb*yt2));
  setTrilinearElement(DA, 2, 3, -lambda * (epsd * (0.628539361*yb*yb2) + epsu * (0.628539361*yb*yt2) );
  setTrilinearElement(DA, 3, 2, -lambda * epsd * 1.257078722 * yb * yb2);
  
  setTrilinearElement(EA, 1, 1, -lambda * epse * (-2./3.) * ytau * ytau2);
  setTrilinearElement(EA, 2, 2, -lambda * epse * 0.628539361 * ytau * ytau2);
  setTrilinearElement(EA, 3, 3, -lambda * espe * 3.771236166 * ytau * ytau2);
  
}