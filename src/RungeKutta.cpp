#include "RungeKutta.hpp"

void Butcher4(struct RungeKuttaCoefficients& RK)
{
  RK.s = 4;

  RK.a.push_back({0.});
  RK.a.push_back({1/2.});
  RK.a.push_back({0., 1/2.});
  RK.a.push_back({0., 0., 1.});

  RK.b.push_back(1/6.);
  RK.b.push_back(1/3.);
  RK.b.push_back(1/3.);
  RK.b.push_back(1/6.);
    
  RK.c.push_back(0.);
  RK.c.push_back(1/2.);
  RK.c.push_back(1/2.);
  RK.c.push_back(1.);
}


void Butcher7(struct RungeKuttaCoefficients& RK, const double nu)
{
  RK.s = 7;

  RK.a.push_back({0.});
  
  RK.a.push_back({nu});
  
  RK.a.push_back({(4.*nu-1.)/(8.*nu), 
                1./(8.*nu)});
  
  RK.a.push_back({(10.*nu-2.)/(27.*nu), 
                2./(27.*nu), 
                (8.*nu)/(27.*nu)});
  
  RK.a.push_back({(56.-77.*nu+(8.-17.*nu)*sqrt(21.))/(392.*nu), 
               (-8.*(7.+sqrt(21.)))/(392.*nu), 
               (48.*(7.+sqrt(21.))*nu)/(392.*nu), 
               (-3.*(21.+sqrt(21.))*nu)/(392.*nu)});
  
  RK.a.push_back({(-5.*(287.*nu-56.-(59.*nu-8.)*sqrt(21.)))/(1960.*nu), 
               (-40.*(7.-sqrt(21.)))/(1960.*nu), 
               (320.*sqrt(21.)*nu)/(1960.*nu), 
               (3.*(21.-121.*sqrt(21.))*nu)/(1960.*nu), 
               (392.*(6.-sqrt(21.))*nu)/(1960.*nu)});
  
  RK.a.push_back({(15.*((30.*nu-8.)-7.*nu*sqrt(21.)))/(180.*nu), 
               120./(180.*nu), 
               (-40.*(5.+7.*sqrt(21.))*nu)/(180.*nu), 
               (63.*(2.+3.*sqrt(21.))*nu)/(180.*nu), 
               (-14.*(49.-9.*sqrt(21.))*nu)/(180.*nu), 
               (70.*(7.+sqrt(21.))*nu)/(180.*nu)});
  
  RK.b.push_back(9/180.);
  RK.b.push_back(0.);
  RK.b.push_back(64/180.);
  RK.b.push_back(0.);
  RK.b.push_back(49/180.);
  RK.b.push_back(49/180.);
  RK.b.push_back(9/180.);

  RK.c.push_back(0.);
  RK.c.push_back(nu);
  RK.c.push_back(1/2.);
  RK.c.push_back(2/3.);
  RK.c.push_back((7.+sqrt(21.))/14.);
  RK.c.push_back((7.-sqrt(21.))/14.);
  RK.c.push_back(1.);
}


void Butcher7(struct RungeKuttaCoefficients& RK, const double lambda, const double mu)
{
  RK.s = 7;

  RK.a.push_back({0.});

  RK.a.push_back({mu});
  
  RK.a.push_back({2/3.-2./(9.*mu),
               2./(9.*mu)});
  
  RK.a.push_back({5/12.-1./(9*mu), 
               1./(9*mu), 
               -1/12.});
  
  RK.a.push_back({17/16.-3./(8.*mu), 
               3./(8.*mu), 
               -3/16., 
               -3/8.});
  
  RK.a.push_back({17/16.-3./(8*mu)+1./(4.*lambda), 
               3./(8.*mu), 
               -3/16.-3./(4.*lambda), 
               -3/8.-3./(2.*lambda), 
               2./lambda});
  
  RK.a.push_back({-27/44.+3./(11.*mu), 
               -3./(11.*mu), 63/44., 
               18/11., 
               4.*((lambda-4.)/11.), 
               -(4.*lambda)/11.});
  
  RK.b.push_back(1/120.);
  RK.b.push_back(0.);
  RK.b.push_back(27/40.);
  RK.b.push_back(27/40.);
  RK.b.push_back((lambda-8.)/15.);
  RK.b.push_back(-lambda/15.);
  RK.b.push_back(11/120.);
  
  RK.c.push_back(0.);
  RK.c.push_back(mu);
  RK.c.push_back(2/3.);
  RK.c.push_back(1/3.);
  RK.c.push_back(1/2.);
  RK.c.push_back(1/2.);
  RK.c.push_back(1.);
}


void Butcher11(struct RungeKuttaCoefficients& RK)
{
  RK.s = 11;

  RK.a.push_back({0.});
  
  RK.a.push_back({1/2.});
  
  RK.a.push_back({1/4., 1/4.});
  
  RK.a.push_back({1/7., (-7.+3.*sqrt(21.))/98., (21.-5.*sqrt(21.))/49.});
  
  RK.a.push_back({(11.-sqrt(21.))/84., 0., (18.-4.*sqrt(21.))/63., (21.+sqrt(21.))/252.});
  
  RK.a.push_back({(5.-sqrt(21.))/48., 0., (9.-sqrt(21.))/36., (-231.-14.*sqrt(21.))/360., (63.+7.*sqrt(21.))/80.});
  
  RK.a.push_back({(10.+sqrt(21.))/42., 0., (-432.-92.*sqrt(21.))/315., (633.+145.*sqrt(21.))/90., (-504.-115.*sqrt(21.))/70., (63.+13.*sqrt(21.))/35.});
  
  RK.a.push_back({1/14., 0., 0., 0., (14.+3.*sqrt(21.))/126., (13.+3.*sqrt(21.))/63., 1/9.});
  
  RK.a.push_back({1/32., 0., 0., 0., (91.+21.*sqrt(21.))/576., 11/72., (-385.+75.*sqrt(21.))/1152., (63.-13.*sqrt(21.))/128.});
  
  RK.a.push_back({1/14., 0., 0., 0., 1/9., (-733.+147.*sqrt(21.))/2205., (515.-111.*sqrt(21.))/504., (-51.+11.*sqrt(21.))/56., (132.-28.*sqrt(21.))/245.});
  
  RK.a.push_back({0., 0., 0., 0., (-42.-7.*sqrt(21.))/18., (-18.-28.*sqrt(21.))/45., (-273.+53.*sqrt(21.))/72., (301.-53.*sqrt(21.))/72., (28.+28.*sqrt(21.))/45., (49.+7.*sqrt(21.))/18.});
  
  RK.b.push_back(1/20.);
  RK.b.push_back(0.);
  RK.b.push_back(0.);
  RK.b.push_back(0.);
  RK.b.push_back(0.);
  RK.b.push_back(0.);
  RK.b.push_back(0.);
  RK.b.push_back(49/180.);
  RK.b.push_back(16/45.);
  RK.b.push_back(49/180.);
  RK.b.push_back(1/20.);
  
  RK.c.push_back(0.);
  RK.c.push_back(1/2.);
  RK.c.push_back(1/2.);
  RK.c.push_back((7.-sqrt(21.))/14.);
  RK.c.push_back((7.-sqrt(21.))/14.);
  RK.c.push_back(1/2.);
  RK.c.push_back((7.+sqrt(21.))/14.);
  RK.c.push_back((7.+sqrt(21.))/14.);
  RK.c.push_back(1/2.);
  RK.c.push_back((7.-sqrt(21.))/14.);
  RK.c.push_back(1.);
  RK.c.push_back(1.);
}
