// ****************************************************************************
// DUT Calibration for RD53A with BDQ53

// Returns the calibration functions for each pixel, where the 
// pixel is identified by the channel: i_col * NumberRows + i_row
// x= col, y = row
// Calibration response function: https://cds.cern.ch/record/2649493/files/CLICdp-Note-2018-008.pdf
// ToT( vcal ) = a vcal + b - \frac{c}{vcal-t}
// Inverse: 
// vcal (ToT ) = (a*t + ToT -b +sqrt((b+a*t-ToT)^2+4*a*c))/(2*a)

// J. Duarte-Campderros (IFCA, Jun 2020)

#include<cmath>
#include<unordered_map>
#include<functional>
#include<stdexcept>
#include<iostream>
#include<fstream>
#include<string>


// The actual function
int __calcurve(float a,float b,float c,float t, int ToT)
{ 
   return std::round((a*t + ToT -b + std::sqrt(std::pow(b+a*t-ToT,2.0)+4.0*a*c))/(2.*a));
}

std::unordered_map<int,std::function<int(int)> > calibration(const std::string & calibration_file, int ntotal_rows)
{
  using namespace std::placeholders;  // for _1, _2, _3...

  std::cout << "Processing gain calibration file: " << calibration_file << std::endl;
  // ASCII text file
  std::ifstream f(calibration_file);
  if( ! f.is_open() )
  {
     throw std::runtime_error(std::string("Invalid file gain: '")+calibration_file+std::string("'"));
  }
  
  // FIXME -- Check valid format
  //
  // FIXME -- Check col-row present??
  
  int col,row;
  float a,b,c,t;

  std::unordered_map<int,std::function<int(int)> >response_vec;
  // Can I get the number of lines, using wc for instance?
  //response_vec.reserve(nlines);
  // Assuming every line is well-formed and consist in
  // col row a b c t
  while(f >> col >> row >> a >> b >> c >> t)
  {
     // FIXME Extract channel from a generic function, to be sure is the same everywhere
     const int channel = col*ntotal_rows + row;
     // The function
     response_vec.emplace( channel,std::bind(__calcurve,a,b,c,t,_1) );      
    
     //std::cout << "Calibration curve for col:" << col << " row: " << row 
     //	<< ":: a= " << a <<" b= " << b << " c= " << c << " t= " << t << std::endl;
  }
  
  return response_vec;
}
