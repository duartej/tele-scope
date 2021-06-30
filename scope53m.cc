
// Daniel Pitzl, DESY, Jun 2018, Mar 2019, Sep 2020
// telescope analysis with RD53A Lin (only)
// module in front

// make scope53m
// scope53m 33095
// needs runs.dat
// needs geo_201x_yy.dat
// needs align_33095.dat from tele
// uses hot_33485.dat from tele
// uses alignDUT_33485.dat
// uses alignMOD_33485.dat
//
// ##########################################
// Adding DUT calibration (RD53A with BDAQ53)
// and other improvements
//
// J. Duarte-Campderros (IFCA, Jun 2020)
//

#include "eudaq/FileReader.hh"
#include "eudaq/PluginManager.hh"

#include <TFile.h>
#include <TH1I.h> // counting
#include <TH1D.h> // weighted counts
#include <TH2I.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TKey.h>

#include <sstream> // stringstream
#include <fstream> // filestream
#include <set>
#include <cmath>
#include <functional>
#include <unordered_map>
#include <stdexcept>
#include <memory>

using namespace std;
using namespace eudaq;

struct pixel {
  int col;
  int row;
  int tot;    // pixel charge
  int frm;    // frame or bunch crossing
  bool pivot; // Mimosa
};

struct cluster {
  vector <pixel> vpix; // pixelsin the cluster
  int size;            // number of pixels (redundant with vpix.size())
  int ncol, nrow;      // number of columns and rows
  int minbc, maxbc;    // earliest and latest pixel bunch cross (clock ticks)
  double col, row;     // center of gravity
  int signal;          // cluster ToT = sum of pixel ToT
  double mindxy;       // isolation from next cluster
};

struct triplet {
  double xm; // mid x [mm]
  double ym; // mid y [mm]
  double zm; // mid z [mm]
  double sx; // slope in x [rad]
  double sy; // slope in y [rad]
  bool lk;       // link flag
  double ttdmin; // distance to next track
  unsigned iA;   // Mimosa cluster index
  unsigned iB;
  unsigned iC;
};

struct CellRelatedHistos
{
    CellRelatedHistos(int cells_shown, float xcell, float ycell) 
    { 
         this->n_cells = cells_shown;
	 this->xcell  = xcell;
	 this->ycell  = ycell;
         this->linnpxvsxmym = new TProfile2D( ("linnpxvsxmym_"+std::to_string(cells_shown)).c_str(), 
		("LIN cluster size vs xmod ymod;x track mod "+std::to_string(int(xcell))+
                   " [#mum];y track mod "+std::to_string(int(ycell))+" [#mum];LIN <cluster size> [pixels]").c_str(),
		50, 0, cells_shown*xcell, 50, 0, cells_shown*ycell, 0, 20 );

 	 this->linqxvsxmym = new TProfile2D( ("linqxvsxmym_"+std::to_string(cells_shown)).c_str(),
		("LIN cluster charge vs xmod ymod;x track mod "+std::to_string(int(xcell))+
                   " [#mum];y track mod "+std::to_string(int(ycell))+" [#mum];LIN <cluster signal> [ToT]").c_str(),
		50, 0, cells_shown*xcell, 50, 0, cells_shown*ycell);

  	 this->effvsxmym = new TProfile2D( ("effvsxmym_"+std::to_string(cells_shown)).c_str(),
		("DUT efficiency vs xmod ymod;x track mod "+std::to_string(int(xcell))+
                   " [#mum];y track mod "+std::to_string(int(ycell))+" [#mum];efficiency").c_str(),
		50, 0, cells_shown*xcell, 50, 0, cells_shown*ycell,-1,2);
    };
    ~CellRelatedHistos() {};

    void add_coordinates(float x,float y) 
    {
         this->xmod = std::fmod(9.0*1e3 + x, this->n_cells*this->xcell);
         this->ymod = std::fmod(9.0*1e3 + y, this->n_cells*this->ycell);
    };
    void fill_clsize(int clsize)
    {
        this->linnpxvsxmym->Fill(this->xmod,this->ymod,clsize);
    };
    void fill_charge(float charge)
    {
         this->linqxvsxmym->Fill(this->xmod,this->ymod,charge);
    };
    void fill_eff(float eff)
    {
         this->effvsxmym->Fill(this->xmod,this->ymod,eff);
    };
	
    int n_cells = 0;
    float xcell = 0;
    float ycell = 0;

    float xmod = 0;
    float ymod = 0;

    TProfile2D * linnpxvsxmym = nullptr;
    TProfile2D * linqxvsxmym = nullptr;
    TProfile2D * effvsxmym = nullptr;
};

// To be use. Get all the cuts here, prepare a constructor
// or a parser to get the cuts from outside
struct Cut
{
    unsigned int bcidmin = 9;
    unsigned int bcidmax = 10;
};

Cut cuts;

//------------------------------------------------------------------------------
// The function returns the calibration functions for each pixel, where the 
// pixel is identified by the channel: i_col * NumberRows + i_row
// x= col, y = row
// Calibration response function: https://cds.cern.ch/record/2649493/files/CLICdp-Note-2018-008.pdf
// ToT( vcal ) = a vcal + b - \frac{c}{vcal-t}
// Inverse: 
// vcal (ToT ) = (a*t + ToT -b +sqrt((b+a*t-ToT)^2+4*a*c))/(2*a)
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

//------------------------------------------------------------------------------
vector <cluster> getClusn( vector <pixel> pb, int fCluCut = 1 ) // 1 = no gap
{
  // returns clusters with pixel coordinates
  // next-neighbour topological clustering (allows fCluCut-1 empty pixels)

  vector <cluster> vc;
  if( pb.size() == 0 ) return vc;

  int* gone = new int[pb.size()];
  for( unsigned i = 0; i < pb.size(); ++i )
    gone[i] = 0;

  unsigned seed = 0;

  while( seed < pb.size() ) {

    // start a new cluster

    cluster c;
    c.vpix.push_back( pb[seed] );
    gone[seed] = 1;

    // let it grow as much as possible:

    int growing;
    do {
      growing = 0;
      for( unsigned i = 0; i < pb.size(); ++i ) {
        if( !gone[i] ){ // unused pixel
          for( unsigned int p = 0; p < c.vpix.size(); ++p ) { // vpix in cluster so far
            int dr = c.vpix.at(p).row - pb[i].row;
            int dc = c.vpix.at(p).col - pb[i].col;
            if( (   dr>=-fCluCut) && (dr<=fCluCut)
		&& (dc>=-fCluCut) && (dc<=fCluCut) ) {
              c.vpix.push_back(pb[i]);
	      gone[i] = 1;
              growing = 1;
              break; // important!
            }
          } // loop over vpix
        } // not gone
      } // loop over all pix
    }
    while( growing );

    // count pixel neighbours:

    for( vector <pixel>::iterator p = c.vpix.begin(); p != c.vpix.end(); ++p ) {
      vector <pixel>::iterator q = p;
      ++q;
      for( ; q != c.vpix.end(); ++q )
	if( abs( p->col - q->col ) <= 1 &&abs( p->row - q->row ) <= 1 ) {
	  ++p->tot;
	  ++q->tot;
	}
    }

    // added all I could. determine position and append it to the list of clusters:

    c.size = c.vpix.size();
    c.col = 0;
    c.row = 0;
    double sumnn = 0;
    int minx = 999;
    int maxx = 0;
    int miny = 999;
    int maxy = 0;
    int minbc = 999;
    int maxbc = 0;

    for( vector<pixel>::iterator p = c.vpix.begin(); p != c.vpix.end(); ++p ) {

      int nn = max( 1, p->tot ); // neighbours
      sumnn += nn;
      c.col += p->col * nn;
      c.row += p->row * nn;

      if( p->col > maxx ) maxx = p->col;
      if( p->col < minx ) minx = p->col;
      if( p->row > maxy ) maxy = p->row;
      if( p->row < miny ) miny = p->row;
      if( p->frm > maxbc ) maxbc = p->frm;
      if( p->frm < minbc ) minbc = p->frm;

    }

    //cout << "(cluster with " << c.vpix.size() << " pixels)" << endl;

    c.col /= sumnn;
    c.row /= sumnn;
    c.signal = sumnn;
    c.ncol = maxx-minx+1;
    c.nrow = maxy-miny+1;
    c.minbc = minbc;
    c.maxbc = maxbc;
    c.mindxy = 999;

    vc.push_back(c); // add cluster to vector

    // look for a new seed = used pixel:

    while( ( ++seed < pb.size() ) && gone[seed] );

  } // while over seeds

  delete[] gone;

  return vc; // vector of clusters

} // getclusn

//------------------------------------------------------------------------------
vector <cluster> getClusq( vector <pixel> pb, int fCluCut = 1 ) // 1 = no gap
{
  // returns clusters with pixel coordinates
  // next-neighbour topological clustering (allows fCluCut-1 empty pixels)

  vector <cluster> vc;
  if( pb.size() == 0 ) return vc;

  bool brck{0};
  if( fCluCut < 0 ) {
    fCluCut *= -1;
    brck = 1;
  }

  int* gone = new int[pb.size()];
  for( unsigned i = 0; i < pb.size(); ++i )
    gone[i] = 0;

  unsigned seed = 0;

  while( seed < pb.size() ) {

    // start a new cluster

    cluster c;
    c.vpix.push_back( pb[seed] );
    gone[seed] = 1;

    // let it grow as much as possible:

    int growing;
    do {
      growing = 0;
      for( unsigned i = 0; i < pb.size(); ++i ) {
        if( !gone[i] ){ // unused pixel
          for( unsigned int p = 0; p < c.vpix.size(); ++p ) { // vpix in cluster so far
            int dr = c.vpix.at(p).row - pb[i].row;
            int dc = c.vpix.at(p).col - pb[i].col;
            if( (   dr>=-fCluCut) && (dr<=fCluCut)
		&& (dc>=-fCluCut) && (dc<=fCluCut) ) {
              c.vpix.push_back(pb[i]);
	      gone[i] = 1;
              growing = 1;
              break; // important!
            }
          } // loop over vpix
        } // not gone
      } // loop over all pix
    }
    while( growing );

    // added all I could. determine position and append it to the list of clusters:

    c.size = c.vpix.size();
    c.col = 0;
    c.row = 0;
    double sumQ = 0;
    int minx = 999;
    int maxx = 0;
    int miny = 999;
    int maxy = 0;
    int minbc = 999;
    int maxbc = 0;

    for( vector<pixel>::iterator p = c.vpix.begin(); p != c.vpix.end(); ++p ) {

      double Qpix = p->tot;

      sumQ += Qpix;

      if( brck && p->row%2 )
	c.col += (p->col+0.5)*Qpix;
      else
	c.col += p->col*Qpix;
      c.row += p->row*Qpix;

      if( p->col > maxx ) maxx = p->col;
      if( p->col < minx ) minx = p->col;
      if( p->row > maxy ) maxy = p->row;
      if( p->row < miny ) miny = p->row;
      if( p->frm > maxbc ) maxbc = p->frm;
      if( p->frm < minbc ) minbc = p->frm;

    }

    //cout << "(cluster with " << c.vpix.size() << " pixels)" << endl;

    if( sumQ > 0 ) {
      c.col /= sumQ;
      c.row /= sumQ;
    }
    else {
      c.col = (*c.vpix.begin()).col;
      c.row = (*c.vpix.begin()).row;
      cout << "GetClus: cluster with signal" << sumQ << endl;
    }

    c.signal = sumQ;
    c.ncol = maxx-minx+1;
    c.nrow = maxy-miny+1;
    c.minbc = minbc;
    c.maxbc = maxbc;
    c.mindxy = 999;

    vc.push_back(c); // add cluster to vector

    // look for a new seed = used pixel:

    while( ( ++seed < pb.size() ) && gone[seed] );

  } // while over seeds

  delete gone;

  return vc; // vector of clusters

} // getclusq

//------------------------------------------------------------------------------
int main( int argc, char* argv[] )
{
  cout << "main " << argv[0] << " called with " << argc << " arguments" << endl;

  if( argc == 1 ) {
    cout << "give run number" << endl;
    return 1;
  }

  // run number = last arg

  string runnum( argv[argc-1] );
  int run = atoi( argv[argc-1] );

  cout << "run " << run << endl;

  if( run < 35034 ) {
    cout << "module at back: should be analysed with scope53.cc" << endl;
    return 2;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // further arguments:

  int fev = 0; // 1st event
  int lev = 999222111; // last event
  bool ldbmod = 0;
  bool use_dut_calibration = true;

  for( int i = 1; i < argc; ++i ) {

    if( !strcmp( argv[i], "-f" ) )
      fev = atoi( argv[++i] ); // 1st event

    if( !strcmp( argv[i], "-l" ) )
      lev = atoi( argv[++i] ); // last event

    if( !strcmp( argv[i], "-m" ) )
      ldbmod = 1; // debug for ref module sync
    
    if( !strcmp( argv[i], "-n" ) )
      use_dut_calibration = false; // no use calibration for DUT

  } // argc

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // runs.dat:

  cout << endl;

  string geoFileName( "geo.dat" );
  double pbeam = 4.8;
  double DUTtilt = 0.5;
  double DUTturn = 0.5; // small turn will not be aligned
  double qwid = 1.5; // The sigma estimated from the Moyal distribution: fitmoyal5('linq0')
  std::string gain_filename_dut;
  double qL = 10; // cutaround Landau peak for best resolution
  double qR = 22;
  int threshold = 0; // offline DUT pixel threshold [TOT]: zero (and 1) = no threshold, 2 or 3 suppress crosstalk
  int minBC = 0; // cut out-of-time noise for honest efficiency
  int maxBC = 31; // look at linpxbc
  int bricked = 0;
  int chip0 = 501;
  int fifty = 0; // default is 100x25
  int rot90 = 0; // default is straight
  int modrun = 0;

  ifstream runsFile( "runs.dat" );

  if( runsFile.bad() || ! runsFile.is_open() ) {
    cout << "Error opening runs.dat" << endl;
    return 1;
  }
  // can there be instructions between if and else ? no

  else {

    cout << "read runs.dat:" << endl;

    string hash( "#" );
    string RUN( "run" );
    string MODRUN( "modrun" );
    string GEO( "geo" );
    string GeV( "GeV" );
    string CHIP( "chip" );
    string QWID( "qsigma_moyal" );
    std::string GAIN_DUT("gain_dut");
    string THRESHOLD( "threshold" );
    string QL("qL");
    string QR("qR");
    string MINBC("minBC");
    string MAXBC("maxBC");
    string TURN( "turn" );
    string TILT( "tilt" );
    string FIFTY( "fifty" );
    string BRICKED( "bricked" );
    string ROT90( "rot90" );
    bool found{0};

    while( ! runsFile.eof() ) {

      string line;
      getline( runsFile, line );

      if( line.empty() ) continue;

      istringstream tokenizer( line );
      string tag;
      tokenizer >> tag; // leading white space is suppressed
      if( tag.substr(0,1) == hash ) // comments start with #
	continue;

      if( tag == RUN )  {
	int ival;
	tokenizer >> ival;
	if( ival == run ) {
	  found = 1;
	  break; // end file reading
	}
      }

      if( tag == MODRUN ) {
	tokenizer >> modrun;
	continue;
      }

      if( tag == TURN ) {
	tokenizer >> DUTturn;
	continue;
      }

      if( tag == TILT ) {
	tokenizer >> DUTtilt;
	continue;
      }
      
      if( tag == GAIN_DUT ) {
	tokenizer >> gain_filename_dut;
	continue;
      }

      if( tag == GEO ) {
	tokenizer >> geoFileName;
	continue;
      }

      if( tag == GeV ) {
	tokenizer >> pbeam;
	continue;
      }

      if( tag == CHIP ) {
	tokenizer >> chip0;
	continue;
      }
      
      if( tag == QWID ) {
	tokenizer >> qwid;
	continue;
      }

      if( tag == FIFTY ) {
	tokenizer >> fifty;
	continue;
      }

      if( tag == BRICKED ) {
	tokenizer >> bricked;
	continue;
      }

      if( tag == ROT90 ) {
	tokenizer >> rot90;
	continue;
      }

      if( tag == QL ) {
	tokenizer >> qL;
	continue;
      }

      if( tag == QR ) {
	tokenizer >> qR;
	continue;
      }

      if( tag == THRESHOLD ) {
	tokenizer >> threshold;
	continue;
      }

      if( tag == MINBC ) {
	tokenizer >> minBC;
	continue;
      }

      if( tag == MAXBC ) {
	tokenizer >> maxBC;
	continue;
      }

      // anything else on the line and in the file gets ignored

    } // while getline

    if( found )
      cout
	<< "  DUT chip " << chip0 << endl
	<< "  modrun " << modrun << endl
	<< "  geo file " << geoFileName << endl
	<< "  nominal DUT turn " << DUTturn << " deg" << endl
	<< "  nominal DUT tilt " << DUTtilt << " deg" << endl
	<< "  DUT chip " << chip0 << endl
        << "  DUT gain file " << gain_filename_dut << std::endl
	<< "  Estimated sigma (from charge global dist.) " << qwid << endl
	<< "  DUT pixel threshold " << threshold << " [ToT or see DUT gain file]" << endl
	<< "  DUT min BC " << minBC << endl
	<< "  DUT max BC " << maxBC << endl
	<< "  DUT qL " << qL << " [ToT or see DUT gain file]" << endl
	<< "  DUT qR " << qR << " [ToT or see DUT gain file]" << endl
	<< "  fifty " << fifty << endl
	<< "  bricked " << bricked << endl
	<< "  rot90 " << rot90 << endl
	<< "  beam " << pbeam << " GeV" << endl
	;
    else {
      cout << "run " << run << " not found in runs.dat" << endl;
      return 1;
    }

  } // runsFile

  runsFile.close();

  const double fTLU = 384E6; // 384 MHz TLU clock

  // 2018 irradiated modules: CERN PS beam spot:

  double xbeam =  3.0; // 509 Lin
  double ybeam = -0.5;
  double rbeam = 3.0;
  /*
    effvsxy->Draw("colz")
    TArc c( -3.0, -0.5, 3.0 ); // circle
    c.SetFillStyle(0);
    c.SetLineWidth(3);
    c.Draw("same");
    c.SetR1(3.0); // ax
    c.SetR2(3.0); // ay
    TArc m( -3.0, -0.5, 0.1 ); // mid
    m.SetFillStyle(1000);
    m.SetFillColor(1);
    m.Draw("same");
  */

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // geometry:

  int nx[9]; // x-pixels per plane
  int ny[9]; // y-pixels per plane
  double sizex[9]; // x size per plane
  double sizey[9]; // y size per plane
  double ptchx[9]; // x-pixel size
  double ptchy[9]; // y-pixel size
  double midx[9]; // x mid
  double midy[9]; // y mid

  double zz[9];

  for( int ipl = 0; ipl < 9; ++ipl )
    nx[ipl] = 0; // missing plane flag

  ifstream geoFile( geoFileName );

  cout << endl;

  if( geoFile.bad() || ! geoFile.is_open() ) {
    cout << "Error opening " << geoFileName << endl;
    return 1;
  }

  cout << "read geometry from " << geoFileName << endl;

  { // open local scope

    string hash( "#" );
    string plane( "plane" );
    string type( "type" );
    string sizexs( "sizex" );
    string sizeys( "sizey" );
    string npixelx( "npixelx" );
    string npixely( "npixely" );
    string zpos( "zpos" );

    int ipl = 0;
    string chiptype;

    while( ! geoFile.eof() ) {

      string line;
      getline( geoFile, line );
      cout << line << endl;

      if( line.empty() ) continue;

      istringstream tokenizer( line );
      string tag;
      tokenizer >> tag; // leading white space is suppressed
      if( tag.substr(0,1) == hash ) // comments start with #
	continue;

      if( tag == plane ) {
	tokenizer >> ipl;
	continue;
      }

      if( ipl < 0 || ipl > 8 ) {
	cout << "geo wrong plane number " << ipl << endl;
	continue;
      }

      if( tag == type ) {
	tokenizer >> chiptype;
	continue;
      }

      if( tag == sizexs ) {
	double val;
	tokenizer >> val;
	sizex[ipl] = val;
	continue;
      }

      if( tag == sizeys ) {
	double val;
	tokenizer >> val;
	sizey[ipl] = val;
	continue;
      }

      if( tag == npixelx ) {
	int val;
	tokenizer >> val;
	nx[ipl] = val;
	continue;
      }

      if( tag == npixely ) {
	int val;
	tokenizer >> val;
	ny[ipl] = val;
	continue;
      }

      if( tag == zpos ) {
	double val;
	tokenizer >> val;
	zz[ipl] = val;
	continue;
      }

      // anything else on the line and in the file gets ignored

    } // while getline

    for( int ipl = 0; ipl < 9; ++ipl ) {
      if( nx[ipl] == 0 ) continue; // missing plane flag
      ptchx[ipl] = sizex[ipl] / nx[ipl]; // pixel size
      ptchy[ipl] = sizey[ipl] / ny[ipl];
      midx[ipl] = 0.5 * sizex[ipl]; // mid of plane
      midy[ipl] = 0.5 * sizey[ipl]; // mid of plane
    }

  } // geo scope

  cout << endl;
  for( int ipl = 1; ipl <= 6; ++ipl )
    cout << ipl << " zz " << zz[ipl] << endl;

  geoFile.close();

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // read Mimosa telescope alignment:

  int aligniteration = 0;
  double alignx[9];
  double aligny[9];
  double alignz[9];
  double rotx[9];
  double roty[9];

  ostringstream alignFileName; // output string stream

  alignFileName << "align_" << run << ".dat";

  ifstream ialignFile( alignFileName.str() );

  cout << endl;

  if( ialignFile.bad() || ! ialignFile.is_open() ) {
    cout << "Error opening " << alignFileName.str() << endl
	 << "  please do: tele -g " << geoFileName << " " << run << endl
	 << endl;
    return 1;
  }
  else {

    cout << "read alignment from " << alignFileName.str() << endl;

    string hash( "#" );
    string iteration( "iteration" );
    string plane( "plane" );
    string shiftx( "shiftx" );
    string shifty( "shifty" );
    string shiftz( "shiftz" );
    string rotxvsy( "rotxvsy" );
    string rotyvsx( "rotyvsx" );

    int ipl = 1;

    while( ! ialignFile.eof() ) {

      string line;
      getline( ialignFile, line );
      cout << line << endl;

      if( line.empty() ) continue;

      istringstream tokenizer( line );
      string tag;
      tokenizer >> tag; // leading white space is suppressed
      if( tag.substr(0,1) == hash ) // comments start with #
	continue;

      if( tag == iteration )
	tokenizer >> aligniteration;

      if( tag == plane )
	tokenizer >> ipl;

      if( ipl < 1 || ipl > 6 ) { // Mimosa
	cout << "align wrong plane number " << ipl << endl;
	continue;
      }

      double val;
      tokenizer >> val;
      if(      tag == shiftx )
	alignx[ipl] = val;
      else if( tag == shifty )
	aligny[ipl] = val;
      else if( tag == shiftz )
	alignz[ipl] = val;
      else if( tag == rotxvsy )
	rotx[ipl] = val;
      else if( tag == rotyvsx )
	roty[ipl] = val;

      // anything else on the line and in the file gets ignored

    } // while getline

  } // alignFile

  ialignFile.close();

  cout << endl;
  for( int ipl = 1; ipl <= 6; ++ipl )
    cout << ipl << " alignz " << alignz[ipl] << endl;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // (re-)create root file:

  ostringstream rootFileName; // output string stream

  rootFileName << "scopeRD" << run << ".root";

  TFile histoFile( rootFileName.str(  ).c_str(  ), "RECREATE" );

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // telescope hot pixels:

  ostringstream hotFileName; // output string stream

  hotFileName << "hot_" << run << ".dat";

  ifstream ihotFile( hotFileName.str() );

  set <int> hotset[9];

  cout << endl;

  if( ihotFile.bad() || ! ihotFile.is_open() ) {
    cout << "no " << hotFileName.str() << " (created by tele)" << endl;
  }
  else {

    cout << "read hot pixel list from " << hotFileName.str() << endl;

    string hash( "#" );
    string plane( "plane" );
    string pix( "pix" );

    int ipl = 0;

    while( ! ihotFile.eof() ) {

      string line;
      getline( ihotFile, line );
      //cout << line << endl;

      if( line.empty() ) continue;

      istringstream tokenizer( line );
      string tag;
      tokenizer >> tag; // leading white space is suppressed
      if( tag.substr(0,1) == hash ) // comments start with #
	continue;

      if( tag == plane )
	tokenizer >> ipl;

      if( ipl < 1 || ipl > 6 ) { // Mimosa
	cout << "hot wrong plane number " << ipl << endl;
	continue;
      }

      if( tag == pix ) {
	int ix, iy;
	tokenizer >> ix;
	tokenizer >> iy;
	int ipx = ix*ny[ipl]+iy;
	hotset[ipl].insert(ipx);
      }

    } // while getline

  } // hotFile

  ihotFile.close();

  for( int ipl = 0; ipl <= 6; ++ipl )
    cout << "  plane " << ipl << ": hot " << hotset[ipl].size() << endl;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // DUT:

  const double wt = atan(1.0) / 45.0; // pi/180 deg

  //double qwid = 1.5; // [ToT] for Moyal in 150 um from x fitmoyal5.C+("linq0")
  //double qxmax = 0.04; // = exp(-qmin/qwid) for qmin = 4.8 ToT lower cutoff ??
  double qxmax = 0.7; // = exp(-80/150) in Vcal units, exp(-800/1500) in ecal units

  int clflag{1}; // default contiguous next-neighbour clustering
  if( bricked ) clflag = -1; // bricked

  int iDUT = 0; // eudaq

  int DUTaligniteration = 0;
  double DUTalignx = 0.0;
  double DUTaligny = 0.0;
  double DUTrot = 0.0;
  double DUTz = 0.5 * ( zz[3] + zz[4] );

  ostringstream DUTalignFileName; // output string stream

  DUTalignFileName << "alignDUT_" << run << ".dat";

  ifstream iDUTalignFile( DUTalignFileName.str() );

  cout << endl;

  if( iDUTalignFile.bad() || ! iDUTalignFile.is_open() ) {
    cout << "no " << DUTalignFileName.str() << ", will bootstrap" << endl;
  }
  else {

    cout << "read DUTalignment from " << DUTalignFileName.str() << endl;

    string hash( "#" );
    string iteration( "iteration" );
    string alignx( "alignx" );
    string aligny( "aligny" );
    string rot( "rot" );
    string tilt( "tilt" );
    string turn( "turn" );
    string dz( "dz" );

    while( ! iDUTalignFile.eof() ) {

      string line;
      getline( iDUTalignFile, line );
      cout << line << endl;

      if( line.empty() ) continue;

      istringstream tokenizer( line );
      string tag;
      tokenizer >> tag; // leading white space is suppressed
      if( tag.substr(0,1) == hash ) // comments start with #
	continue;

      if( tag == iteration )
	tokenizer >> DUTaligniteration;

      double val;
      tokenizer >> val;
      if(      tag == alignx )
	DUTalignx = val;
      else if( tag == aligny )
	DUTaligny = val;
      else if( tag == rot )
	DUTrot = val;
      else if( tag == tilt )
	DUTtilt = val;
      else if( tag == turn )
	DUTturn = val;
      else if( tag == dz )
	DUTz = val + zz[3];

      // anything else on the line and in the file gets ignored

    } // while getline

  } // alignFile

  iDUTalignFile.close();

  double DUTalignx0 = DUTalignx; // at time 0
  double DUTaligny0 = DUTaligny;

  // normal vector on DUT surface:
  // N = ( 0, 0, -1 ) on DUT, towards -z
  // transform into tele system:
  // tilt alpha around x
  // turn omega around y

  const double co = cos( DUTturn*wt );
  const double so = sin( DUTturn*wt );
  const double ca = cos( DUTtilt*wt );
  const double sa = sin( DUTtilt*wt );
  const double cf = cos( DUTrot );
  const double sf = sin( DUTrot );

  const double Nx =-ca*so;
  const double Ny = sa;
  const double Nz =-ca*co;

  const double norm = cos( DUTturn*wt ) * cos( DUTtilt*wt ); // length of Nz

  set <int> deadset;

  ostringstream DUTdeadFileName; // output string stream

  DUTdeadFileName << "dead" << chip0 << ".dat";

  ifstream iDUTdeadFile( DUTdeadFileName.str() );

  if( iDUTdeadFile.bad() || ! iDUTdeadFile.is_open() ) {
    cout << "no " << DUTdeadFileName.str() << endl;
  }
  else {

    cout << "read DUT dead pixel list from " << DUTdeadFileName.str() << endl << endl;

    TH2I * deadxyHisto = new
      TH2I( "deadxy", "DUT dead pixels;x [mm];y [mm];DUT dead pixels",
	    400, 0, 400, 192, 0, 192 ); // bin = pix

    string hash( "#" );
    string COL( "col" );

    while( ! iDUTdeadFile.eof() ) {

      string line;
      getline( iDUTdeadFile, line );
      //cout << line << endl;

      if( line.empty() ) continue;

      istringstream tokenizer( line );

      string tag;
      tokenizer >> tag; // leading white space is suppressed

      if( tag.substr(0,1) == hash ) { // comments start with #
	cout << line << endl;
	continue;
      }

      if( tag == COL ) {

	int col;
	tokenizer >> col;
	cout << "col " << setw(3) << col << ":";

	while( ! tokenizer.eof() ) {

	  int row;
	  tokenizer >> row;
	  cout << "  " << row;

	  // sensor pixels:

	  int ix = col;
	  int iy = row;
	  if( !fifty ) {
	    ix = col/2; // sensor 100
	    if( col%2 == 1 )
	      iy = 2*row + 1; // sensor 25
	    else
	      iy = 2*row + 0;
	  }

	  int ipx = ix * 384 + iy;
	  deadset.insert(ipx);

	  deadxyHisto->Fill( col, row );

	}

	cout << endl;

      } // tag

    } // while getline

  } // deadFile

  iDUTdeadFile.close();

  cout << "DUT dead " << deadset.size() << endl;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // DUT hot pixels:

  cout << endl;

  ostringstream DUThotFileName; // output string stream

  DUThotFileName << "hotDUT_" << run << ".dat";

  ifstream iDUThotFile( DUThotFileName.str() );
 
  bool create_duthotfile = false;
  if( iDUThotFile.bad() || ! iDUThotFile.is_open() ) {
    cout << "no " << DUThotFileName.str() << endl;
    // Calculate the Hot pixels
    create_duthotfile=true;
  }
  else {

    cout << "read DUT hot pixel list from " << DUThotFileName.str() << endl;

    string hash( "#" );
    string pix( "pix" );

    while( ! iDUThotFile.eof() ) {

      string line;
      getline( iDUThotFile, line );
      //cout << line << endl;

      if( line.empty() ) continue;

      stringstream tokenizer( line );
      string tag;
      tokenizer >> tag; // leading white space is suppressed
      if( tag.substr(0,1) == hash ) // comments start with #
	continue;

      if( tag == pix ) {

	int ix, iy;
	tokenizer >> ix; // ROC col
	tokenizer >> iy; // ROC row
	int ipx = ix * ny[iDUT] + iy;
	hotset[iDUT].insert(ipx);

	int col = ix;
	int row = iy;
	if( !fifty ) {
	  col = ix/2; // sensor 100
	  if( ix%2 == 1 )
	    row = 2*iy + 1; // sensor 25
	  else
	    row = 2*iy + 0;
	}
	ipx = col * 384 + row; // sensor
	deadset.insert(ipx);

      }

    } // while getline

  } // hotFile

  iDUThotFile.close();

  cout << "DUT hot " << hotset[iDUT].size() << endl;

  cout << endl;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // MOD:

  string modfileA = "mod/run" + to_string(modrun) + "A.out"; // C++11
  if( run >= 38027 )
    modfileA = "mod/run" + to_string(modrun) + "_ch0.out"; // C++11
  cout << "try to open  " << modfileA;
  ifstream Astream( modfileA.c_str() );
  if( !Astream ) {
    cout << " : failed " << endl;
    modrun = 0; // flag
  }
  else
    cout << " : succeed " << endl;

  string modfileB = "mod/run" + to_string(modrun) + "B.out";
  if( run >= 38027 )
    modfileB = "mod/run" + to_string(modrun) + "_ch1.out"; // C++11
  cout << "try to open  " << modfileB;
  ifstream Bstream( modfileB.c_str() );
  if( !Bstream ) {
    cout << " : failed " << endl;
    modrun = 0;
  }
  else
    cout << " : succeed " << endl;

  string sl;
  int mev = 0;
  if( fev ) cout << "MOD skip " << fev << endl;
  while( mev < fev ) {
    getline( Astream, sl ); // fast forward
    getline( Bstream, sl );
    ++mev;
  }

  int iMOD = 7;

  int MODaligniteration = 0;
  double MODalignx = 0;
  double MODaligny = 0;
  double MODrot = 0;
  double MODtilt =  19.2; // [deg]
  double MODturn = -27.0; // [deg]
  double MODz = -45 + zz[1];

  ostringstream MODalignFileName; // output string stream

  MODalignFileName << "alignMOD_" << run << ".dat";

  ifstream iMODalignFile( MODalignFileName.str() );

  cout << endl;

  if( iMODalignFile.bad() || ! iMODalignFile.is_open() ) {
    cout << "no " << MODalignFileName.str() << ", will bootstrap" << endl;
  }
  else {

    cout << "read MODalignment from " << MODalignFileName.str() << endl;

    string hash( "#" );
    string iteration( "iteration" );
    string alignx( "alignx" );
    string aligny( "aligny" );
    string rot( "rot" );
    string tilt( "tilt" );
    string turn( "turn" );
    string dz( "dz" );

    while( ! iMODalignFile.eof() ) {

      string line;
      getline( iMODalignFile, line );
      cout << line << endl;

      if( line.empty() ) continue;

      istringstream tokenizer( line );
      string tag;
      tokenizer >> tag; // leading white space is suppressed
      if( tag.substr(0,1) == hash ) // comments start with #
	continue;

      if( tag == iteration )
	tokenizer >> MODaligniteration;

      double val;
      tokenizer >> val;
      if(      tag == alignx )
	MODalignx = val;
      else if( tag == aligny )
	MODaligny = val;
      else if( tag == rot )
	MODrot = val;
      else if( tag == tilt )
	MODtilt = val;
      else if( tag == turn )
	MODturn = val;
      else if( tag == dz )
	MODz = val + zz[1];

      // anything else on the line and in the file gets ignored

    } // while getline

  } // alignFile

  iMODalignFile.close();

  // normal vector on MOD surface:
  // N = ( 0, 0, -1 ) on MOD, towards -z
  // transform into tele system:
  // tilt alpha around x
  // turn omega around y

  const double com = cos( MODturn*wt );
  const double som = sin( MODturn*wt );
  const double cam = cos( MODtilt*wt );
  const double sam = sin( MODtilt*wt );
  const double cfm = cos( MODrot );
  const double sfm = sin( MODrot );

  const double Nxm =-cam*som;
  const double Nym = sam;
  const double Nzm =-cam*com;

  const double normm = cos( MODturn*wt ) * cos( MODtilt*wt ); // length of Nz

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // DUT gain:
  int nrows_in_dut = ny[iDUT];
  if( !fifty ) 
  {
     //nrows_in_dut = 2*ny[iDUT];
     nrows_in_dut = 384;
  }
  auto calibration_curves = calibration(gain_filename_dut,nrows_in_dut);
  std::cout << "Loaded calibration curves for DUT. Active pixels: " 
	<< calibration_curves.size() << std::endl;
  histoFile.cd();

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // book histos:
  int qnbins = 16;
  double qvalmax = 16;
  double qscale = 1.0;
  std::string qunit("ToT");
  
  // FIXME: if vcal or electrons or ToT --> from runs.dat (or calibration)  
  std::string charge_units("elec");
  if( use_dut_calibration)
  {
    if( charge_units == "elec" )
    {
       qnbins  = 1000;
       qvalmax = 80000;
       qscale  = 1500;
       qunit = "Electrons";
    }
    else if( charge_units == "vcal" )
    { 
       qnbins  = 1000;
       qvalmax = 8000;
       qscale  = 150;
       qunit = "V_{cal}";
    }
  }


  double f = 4.8/pbeam;

  TH1I t1Histo( "t1", "event time;event time [s];events / 10 ms", 100, 0, 1 );
  TH1I t2Histo( "t2", "event time;event time [s];events / s", 300, 0, 300 );
  TH1I t3Histo( "t3", "event time;event time [s];events / 10 s", 150, 0, 1500 );
  TH1I t4Histo( "t4", "event time;event time [s];events /10 s", 600, 0, 6000 );
  TH1I t5Histo( "t5", "event time;event time [s];events / 60 s", 1100, 0, 66000 );
  TH1I t6Histo( "t6", "event time;event time [h];events / 3 min", 1000, 0, 50 );

  TH1I dtusHisto( "dtus", "time between events;time between events [us];events", 100, 0, 1000 );
  TH1I dtmsHisto( "dtms", "time between events;time between events [ms];events", 100, 0, 1000 );
  TH1I dt373Histo( "dt373", "time between events;time between events mod 373;events", 400, -0.5, 399.5 );
  TH1I dt374Histo( "dt374", "time between events;time between events mod 374;events", 400, -0.5, 399.5 );
  TH1I dt375Histo( "dt375", "time between events;time between events mod 375;events", 400, -0.5, 399.5 );
  TH1I dt376Histo( "dt376", "time between events;time between events mod 376;events", 400, -0.5, 399.5 );
  TH1I dt377Histo( "dt377", "time between events;time between events mod 377;events", 400, -0.5, 399.5 );
  TProfile dt375vsdt( "dt375vsdt", "dt vs dt;time between events [us];<time between events mod 375>",
		      200, 0, 2000, 0, 400 );

  TH1I hpivot[9];
  TH1I hnpx[9];
  TH1I hnpxmsk[9];
  TH1I hnframes[9];
  TH1I hcol[9];
  TH1I hrow[9];
  TH2I * hmap[9];

  TH1I hncl[9];
  TH1I hsiz[9];
  TH1I hncol[9];
  TH1I hnrow[9];
  TH1I hdxy[9];

  for( int ipl = 0; ipl <= 7; ++ipl ) {

    int nbx = 400; // RD53
    int nby = 192;
    int mx = nbx;
    int my = nby;

    if( ipl >= 1 ) {
      nbx = nx[1]/2; // Mimosa
      nby = ny[1]/2; // for next round
      mx = nx[1];
      my = ny[1];
    }
    if( ipl == 7 ) { // MOD
      nbx = 8*(52+2);
      nby = 162; // with big pix
      mx = nbx;
      my = nby;
    }

    hpivot[ipl] = TH1I( Form( "pivot%i", ipl ),
			Form( "%i pivot;pivot;plane %i events", ipl, ipl ),
			1000, 0, 10*1000 );

    hnpx[ipl] = TH1I( Form( "npx%i", ipl ),
		      Form( "%i pixel per event before masking;pixels/event;plane %i events", ipl, ipl ),
		      201, -0.5, 200.5 );
    hnpxmsk[ipl] = TH1I( Form( "npxmsk%i", ipl ),
			 Form( "%i pixel per event after masking;pixels/event;plane %i events", ipl, ipl ),
			 201, -0.5, 200.5 );

    hnframes[ipl] = TH1I( Form( "nframes%i", ipl ),
			  Form( "%i frames per event;frames;plane %i events", ipl, ipl ),
			  51, -0.5, 50.5 );

    hcol[ipl] = TH1I( Form( "col%i", ipl ),
		      Form( "%i col;col;plane %i pixels", ipl, ipl ),
		      nbx, 0, mx );
    hrow[ipl] = TH1I( Form( "row%i", ipl ),
		      Form( "%i row;row;plane %i pixels", ipl, ipl ),
		      nby, 0, my );
    hmap[ipl] = new TH2I( Form( "map%i", ipl ),
			  Form( "%i map;col;row;plane %i pixels", ipl, ipl ),
			  nbx, 0, mx, nby, 0, my );

    hncl[ipl] = TH1I( Form( "ncl%i", ipl ),
		      Form( "plane %i cluster per event;cluster;plane %i events", ipl, ipl ),
		      51, -0.5, 50.5 );
    hsiz[ipl] = TH1I( Form( "clsz%i", ipl ),
		      Form( "%i cluster size;pixels/cluster;plane %i clusters", ipl, ipl ),
		      51, -0.5, 50.5 );
    hncol[ipl] = TH1I( Form( "ncol%i", ipl ),
		       Form( "%i cluster size x;columns/cluster;plane %i clusters", ipl, ipl ),
		       21, -0.5, 20.5 );
    hnrow[ipl] = TH1I( Form( "nrow%i", ipl ),
		       Form( "%i cluster size y;rows/cluster;plane %i clusters", ipl, ipl ),
		       21, -0.5, 20.5 );
    hdxy[ipl] = TH1I( Form( "dxy%i", ipl ),
		      Form( "%i cluster distance;cluster distance [pixels];plane %i clusters", ipl, ipl ),
		      100, 0, 10 );

  } // planes

  TProfile dutnpxvsev( "dutnpxvsev",
		       "DUT pixels vs events;events;DUT pixels / 1000",
		       500, 0, 500E3 );

  hrow[8] = TH1I( Form( "row%i", 8 ),
		  Form( "%i row;row;plane %i pixels", 8, 8 ),
		  192, 0, 192 );

  hmap[8] = new TH2I( Form( "map%i", 8 ),
		      Form( "%i map before masking;col;row;plane %i pixels", 8, 8 ),
		      400, 0, 400, 192, 0, 192 );
  TH1I dutpxq8Histo( "dutpxq8",
		    "DUT pixel signal;DUT pixel signal [ToT];all DUT pixels",
		    16, -0.5, 15.5 );

  TH1I dutpxcol0Histo( "dutpxcol0",
		      "DUT pixel column before threshold;DUT pixel column;pixels",
		      400, 0, 400 );
  TH1I linpxqHisto( "linpxq",
		    "Lin pixel signal;Lin pixel signal [ToT];Lin pixels",
		    qnbins, 0, qvalmax );

  TProfile dutpxqvsx( "dutpxqvsx",
		      "DUT pixel signal vs x;column;<pixel signal> [ToT]",
		      400, 0, 400, 0, qvalmax );
  TProfile2D * dutpxqvsxy = new
    TProfile2D( "dutpxqvsxy",
		"DUT pixel signal map;column;row;<pixel signal> [ToT]",
		400, 0, 400, 192, 0, 192, 0, qvalmax );

  TH1I dutpxbcHisto( "dutpxbc",
		     "DUT pixel BC;DUT pixel BC;DUT pixels without masking",
		     32, -0.5, 31.5 );
  TH1I dutpxbcmHisto( "dutpxbcm",
		      "DUT pixel BC;DUT pixel BC;DUT pixels with masking",
		      32, -0.5, 31.5 );

  TH1I dutpxbc01Histo( "dutpxbc01",
		       "DUT pixel BC;DUT pixel BC;DUT pixels ToT 1",
		       32, -0.5, 31.5 );
  TH1I dutpxbc15Histo( "dutpxbc15",
		       "DUT pixel BC;DUT pixel BC;DUT pixels ToT 15",
		       32, -0.5, 31.5 );
  TH1I dutpxbc08Histo( "dutpxbc08",
		       "DUT pixel BC;DUT pixel BC;DUT pixels ToT mid",
		       32, -0.5, 31.5 );

  TProfile2D * dutpxbcmap = new
    TProfile2D( "dutpxbcmap", "DUT pixel BC;col;row;<DUT pixel BC>",
		136, 127.5, 263.5, 192, -0.5, 191 );
  TH2I * dutpxmapbc7 = new
    TH2I( "dutpxmapbc7", "DUT pixel BC 7;col;row;DUT BC 7 pixel",
		136, 127.5, 263.5, 192, -0.5, 191 );
  TH2I * dutpxmapbc9 = new
    TH2I( "dutpxmapbc9", "DUT pixel BC 9;col;row;DUT BC 9 pixel",
		136, 127.5, 263.5, 192, -0.5, 191 );

  TProfile2D * dutpxbcmap8 = new
    TProfile2D( "dutpxbcmap8", "DUT pixel BC;col%8;row%8;<DUT pixel BC>",
		8, -0.5, 7.5, 8, -0.5, 7.5 );
  TH2I * dutpxmap8 = new
    TH2I( "dutpxmap8", "DUT pixel map;col%8;row%8;cool pixels",
	  8, -0.5, 7.5, 8, -0.5, 7.5 );
  TH2I * dutpxmap8bc10 = new
    TH2I( "dutpxmap8bc10", "DUT pixel map;col%8;row%8;cool pixels",
	  8, -0.5, 7.5, 8, -0.5, 7.5 );

  TProfile dutpxbcvst( "dutpxbcvst", "DUT pixel BC;time [trigger];<piixel BC>",
		       500, 0, 500e3 );

  TH1I linpxbcHisto( "linpxbc",
		    "Lin pixel BC;Lin pixel BC;Lin pixels",
		    32, -0.5, 31.5 );

  TH1I dutpxcol9Histo( "dutpxcol9",
		      "DUT pixel column after threshold;DUT pixel column;pixels",
		      400, 0, 400 );

  TH2I * dutmodxxHisto = new
    TH2I( "dutmodxx", "Mod vs DUT x-x;x_{DUT} [mm];x_{MOD} [mm];cluster pairs",
	  220, -11, 11, 660, -33, 33 );
  TH2I * dutmodyyHisto = new
    TH2I( "dutmodyy", "Mod vs DUT y-y;y_{DUT} [mm];y_{MOD} [mm];cluster pairs",
	  120, -6, 6, 160, -8, 8 );

  TH1I dutdpxHisto( "dutdpx",
		    "DUT cluster distance;cluster distance [pixels];DUT cluster pairs",
		    100, -0.5, 99.5 );

  // triplets:

  TH1I hdx13( "dx13", "1-3 dx;1-3 dx [mm];cluster pairs", 100, -f, f );
  TH1I hdy13( "dy13", "1-3 dy;1-3 dy [mm];cluster pairs", 100, -f, f );

  TH1I htridx( "tridx", "triplet dx;triplet dx [mm];triplets", 100, -0.1, 0.1 );
  TH1I htridy( "tridy", "triplet dy;triplet dy [mm];triplets", 100, -0.1, 0.1 );

  TH1I htridxc( "tridxc", "triplet dx;triplet dx [mm];triplets", 100, -0.05, 0.05 );
  TH1I htridyc( "tridyc", "triplet dy;triplet dy [mm];triplets", 100, -0.05, 0.05 );

  TProfile tridxvsx( "tridxvsx",
		     "triplet dx vs x;triplet x [mm];<triplet #Deltax> [mm]",
		     120, -12, 12, -0.05, 0.05 );
  TProfile trimadxvsx( "trimadxvsx",
		       "triplet MADx vs x;triplet x [mm];triplet MAD(#Deltax) [#mum]",
		       120, -12, 12, 0, 50 );
  TProfile tridxvsy( "tridxvsy",
		     "triplet dx vs y;triplet y [mm];<triplet #Deltax> [mm]",
		     110, -5.5, 5.5, -0.05, 0.05 );
  TProfile tridxvstx( "tridxvstx",
		      "triplet dx vs slope x;triplet slope x [rad];<triplet #Deltax> [mm]",
		      60, -0.003, 0.003, -0.05, 0.05 );
  TProfile tridxvst3( "tridxvst3",
		      "triplet dx vs time;time [s];<triplet #Deltax> [mm] / 10s",
		      300, 0, 3000, -0.05, 0.05 );
  TProfile tridxvst5( "tridxvst5",
		      "triplet dx vs time;time [s];<triplet #Deltax> [mm] / min",
		      1100, 0, 66000, -0.05, 0.05 );

  TProfile tridyvsx( "tridyvsx",
		     "triplet dy vs x;triplet x [mm];<triplet #Deltay> [mm]",
		     110, -11, 11, -0.05, 0.05 );
  TProfile tridyvsty( "tridyvsty",
		      "triplet dy vs slope y;triplet slope y [rad];<triplet #Deltay> [mm]",
		      60, -0.003, 0.003, -0.05, 0.05 );
  TProfile tridyvst3( "tridyvst3",
		      "triplet dy vs time;time [s];<triplet #Deltay> [mm] / 10s",
		      300, 0, 3000, -0.05, 0.05 );
  TProfile tridyvst5( "tridyvst5",
		      "triplet dy vs time;time [h];<triplet #Deltay> [mm] / min",
		      1100, 0, 66000, -0.05, 0.05 );

  TH1I trixHisto( "trix", "triplets x;x [mm];triplets",
		  240, -12, 12 );
  TH1I triyHisto( "triy", "triplets y;y [mm];triplets",
		  120, -6, 6 );
  TH2I * trixyHisto = new
    TH2I( "trixy", "triplets x-y;x [mm];y [mm];triplets",
	  240, -12, 12, 120, -6, 6 );
  TH1I tritxHisto( "tritx", "triplet slope x;slope x [mrad];triplets",
		   100, -5*f, 5*f );
  TProfile tritxwvsx( "tritxwvsx",
		      "triplet slope width vs x;triplet x [mm];MAD(slope_{x}) [mrad]",
		      240, -12, 12, 0, 5*f );
  TH1I trityHisto( "trity", "triplet slope y;slope y [mrad];triplets",
		   100, -5*f, 5*f );

  TH1I ntriHisto( "ntri", "triplets;triplets;events", 51, -0.5, 50.5 );

  // driplets:

  TH1I hdx46( "dx46", "4-6 dx;4-6 dx [mm];cluster pairs", 100, -f, f );
  TH1I hdy46( "dy46", "4-6 dy;4-6 dy [mm];cluster pairs", 100, -f, f );

  TH1I hdridx( "dridx", "driplet dx;driplet dx [mm];driplets", 100, -0.1*f, 0.1*f );
  TH1I hdridy( "dridy", "driplet dy;driplet dy [mm];driplets", 100, -0.1*f, 0.1*f );

  TH1I hdridxc( "dridxc", "driplet dx;driplet dx [mm];driplets", 100, -0.1*f, 0.1*f );
  TH1I hdridyc( "dridyc", "driplet dy;driplet dy [mm];driplets", 100, -0.1*f, 0.1*f );

  TProfile dridxvsy( "dridxvsy",
		     "driplet dx vs y;driplet yB [mm];<driplets #Deltax> [mm]",
		     110, -5.5, 5.5, -0.05*f, 0.05*f );
  TProfile dridyvsx( "dridyvsx",
		     "driplet dy vs x;driplet xB [mm];<driplets #Deltay> [mm]",
		     110, -11, 11, -0.05*f, 0.05*f );

  TProfile dridxvstx( "dridxvstx",
		      "driplet dx vs slope x;driplet slope x [rad];<driplets #Deltax> [mm]",
		      60, -0.003, 0.003, -0.05*f, 0.05*f );
  TProfile dridyvsty( "dridyvsty",
		      "driplet dy vs slope y;driplet slope y [rad];<driplets #Deltay> [mm]",
		      60, -0.003, 0.003, -0.05*f, 0.05*f );
  TProfile dridxvst3( "dridxvst3",
		      "driplet dx vs time;time [s];<driplet #Deltax> [mm] / 10s",
		      300, 0, 3000, -0.05, 0.05 );
  TProfile dridxvst5( "dridxvst5",
		      "driplet dx vs time;time [s];<driplet #Deltax> [mm] / min",
		      1100, 0, 66000, -0.05, 0.05 );
  TProfile dridyvst3( "dridyvst3",
		      "driplet dy vs time;time [s];<driplet #Deltay> [mm] / 10s",
		      300, 0, 3000, -0.05, 0.05 );
  TProfile dridyvst5( "dridyvst5",
		      "driplet dy vs time;time [h];<driplet #Deltay> [mm] / min",
		      1100, 0, 66000, -0.05, 0.05 );

  TH1I drixHisto( "drix", "driplets x;x [mm];driplets",
		  240, -12, 12 );
  TH1I driyHisto( "driy", "driplets x;y [mm];driplets",
		  120, -6, 6 );
  TH2I * drixyHisto = new
    TH2I( "drixy", "driplets x-y;x [mm];y [mm];driplets",
	  240, -12, 12, 120, -6, 6 );
  TH1I dritxHisto( "dritx", "driplet slope x;slope x [rad];driplets",
		   100, -0.005*f, 0.005*f );
  TH1I drityHisto( "drity", "driplet slope y;slope y [rad];driplets",
		   100, -0.005*f, 0.005*f );

  TH1I ndriHisto( "ndri", "driplets;driplets;events", 51, -0.5, 50.5 );

  // MOD vs triplets:

  TH1I modcolHisto( "modcol", "MOD col;MOD col;MOD cluster", 216, 0, 432 );
  TH1I modrowHisto( "modrow", "MOD row;MOD row;MOD cluster", 162, 0, 162 );

  TH1I modxHisto( "modx", "MOD x;MOD cluster x [mm];MOD clusters", 216, -32.4, 32.4 );
  TH1I modyHisto( "mody", "MOD y;MOD cluster y [mm];MOD clusters", 200, -10, 10 );

  TH1I modmindxHisto( "modmindx",
		      "min MOD - triplet x;min MOD cluster - triplet #Deltax [mm];MOD clusters",
		      200, -0.5, 0.5 );
  TH1I modmindyHisto( "modmindy",
		      "min MOD - triplet y;min MOD cluster - triplet #Deltay [mm];MOD clusters",
		      200, -0.5, 0.5 );

  TH1I modxlkHisto( "modxlk",
		    "MOD hits linked to triplet;x in MOD [mm];linked triplets",
		    216, -32.4, 32.4 );
  TH1I modylkHisto( "modylk",
		    "MOD hits linked to triplet;y in MOD [mm];linked triplets",
		    200, -10, 10 );

  TH1I ttdminmod1Histo( "ttdminmod1",
		     "triplet isolation at MOD;triplet at MOD min #Delta_{xy} [mm];triplet pairs",
		     100, 0, 1 );
  TH1I ttdminmod2Histo( "ttdminmod2",
		     "triplet isolation at MOD;triplet at MOD min #Delta_{xy} [mm];triplet pairs",
		     150, 0, 15 );

  TH1I trixmHisto( "trixm",
		   "triplet at MOD x;triplet x at MOD [mm];triplets",
		   216, -32.4, 32.4 );
  TH1I triymHisto( "triym",
		   "triplet at MOD y;triplet y at MOD [mm];triplets",
		   120, -6, 6 );

  TH1I modsxaHisto( "modsxa",
		    "MOD + triplet x;MOD cluster + triplet #Sigmax [mm];MOD clusters",
		    1280, -32, 32 );
  TH1I moddxaHisto( "moddxa",
		    "MOD - triplet x;MOD cluster - triplet #Deltax [mm];MOD clusters",
		    1280, -32, 32 );

  TH1I modsyaHisto( "modsya",
		    "MOD + triplet y;MOD cluster + triplet #Sigmay [mm];MOD clusters",
		    320, -8, 8 );
  TH1I moddyaHisto( "moddya",
		    "MOD - triplet y;MOD cluster - triplet #Deltay [mm];MOD clusters",
		    320, -8, 8 );

  TH1I moddxHisto( "moddx",
		   "MOD - triplet x;MOD cluster - triplet #Deltax [mm];MOD clusters",
		   500, -2.5, 2.5 );
  TH1I moddxcHisto( "moddxc",
		    "MOD - triplet x;MOD cluster - triplet #Deltax [mm];MOD clusters",
		    200, -0.5, 0.5 );
  TProfile moddxvsx( "moddxvsx",
		     "MOD #Deltax vs x;x track [mm];<cluster - triplet #Deltax> [mm]",
		     216, -32.4, 32.4, -2.5, 2.5 );
  TProfile moddxvsy( "moddxvsy",
		     "MOD #Deltax vs y;y track [mm];<cluster - triplet #Deltax> [mm]",
		     160, -8, 8, -2.5, 2.5 );
  TProfile modmadxvsy( "modmadxvsy",
		       "MOD MAD #Deltax vs y;y track [mm];<cluster - triplet MAD #Deltax> [mm]",
		       160, -8, 8, 0, 0.5 );
  TProfile moddxvstx( "moddxvstx",
		      "MOD #Deltax vs #theta_{x};x track slope [rad];<cluster - triplet #Deltax> [mm]",
		      80, -0.002, 0.002, -2.5, 2.5 );
  TProfile moddxvst5( "moddxvst5",
		      "MOD dx vs time;time [s];<MOD #Deltax> [mm] / min",
		      1100, 0, 66000, -0.05, 0.05 );

  TH1I moddyHisto( "moddy",
		   "MOD - triplet y;MOD cluster - triplet #Deltay [mm];MOD clusters",
		   200, -0.5, 0.5 );
  TH1I moddycHisto( "moddyc",
		    "MOD - triplet y;MOD cluster - triplet #Deltay [mm];MOD clusters",
		    200, -0.5, 0.5 );
  TH1I moddycqHisto( "moddycq",
		     "MOD - triplet y Landau peak;MOD cluster - triplet #Deltay [mm];Landau peak MOD clusters",
		     500, -0.5, 0.5 );
  TProfile moddyvsx( "moddyvsx",
		     "MOD #Deltay vs x;x track [mm];<cluster - triplet #Deltay> [mm]",
		     216, -32.4, 32.4, -0.5, 0.5 );
  TProfile moddyvsy( "moddyvsy",
		     "MOD #Deltay vs y;y track [mm];<cluster - triplet #Deltay> [mm]",
		     160, -8, 8, -0.5, 0.5 );
  TProfile modmadyvsy( "modmadyvsy",
		       "MOD MAD#Deltay vs y;y track [mm];<cluster - triplet MAD #Deltay> [mm]",
		       160, -8, 8, 0, 0.25 );
  TProfile moddyvsty( "moddyvsty",
		      "MOD #Deltay vs #theta_{y};y track slope [rad];<cluster - triplet #Deltay> [mm]",
		      80, -0.002, 0.002, -0.5, 0.5 );
  TProfile moddyvst5( "moddyvst5",
		      "MOD dy vs time;time [h];<MOD #Deltay> [mm] / min",
		      1100, 0, 66000, -0.05, 0.05 );

  TH1I modnpxHisto( "modnpx",
		    "MOD linked clusters;MOD cluster size [pixels];linked MOD cluster",
		    20, 0.5, 20.5 );

  TH1I modqHisto( "modq",
		  "MOD linked clusters;MOD cluster charge [ke];linked MOD cluster",
		  80, 0, 80 );
  TH1I modq0Histo( "modq0",
		   "MOD linked clusters;MOD normal cluster charge [ke];linked MOD cluster",
		   80, 0, 80 );

  TProfile2D * modnpxvsxmym = new
    TProfile2D( "modnpxvsxmym",
		"MOD cluster size vs xmod ymod;x track mod 300 [#mum];y track mod 200 [#mum];MOD <cluster size> [pixels]",
		120, 0, 300, 80, 0, 200, 0, 20 );

  TH1I modcollkHisto( "modcollk",
		      "MOD linked col;MOD linked col;linked MOD cluster",
		      216, 0, 432 );
  TH1I modrowlkHisto( "modrowlk",
		      "MOD linked row;MOD linked row;linked MOD cluster",
		      162, 0, 162 );

  TH1I trixmlkHisto( "trixmlk",
		     "linked triplet at MOD x;triplet x at MOD [mm];linked triplets",
		     216, -32.4, 32.4 );
  TH1I triymlkHisto( "triymlk",
		     "linked triplet at MOD y;triplet y at MOD [mm];linked triplets",
		     120, -6, 6 );

  TH1I tripxfrmlkHisto( "tripxfrmlk",
			"triplet linked pixel frame;frame;linked triplet pixels",
			2, -0.5, 1.5 );
  TH1I tripxpivlkHisto( "tripxpivlk",
			"triplet linked pixel pivot;pivot;linked triplet pixels",
			2, -0.5, 1.5 );

  TH1I trixmodHisto( "trixmod",
		     "triplet with MOD x;triplet x at MOD [mm];triplet-Mod links",
		     216, -32.4, 32.4 );
  TH1I triymodHisto( "triymod",
		     "triplet with MOD y;triplet y at MOD [mm];triplet Mod- links",
		     120, -6, 6 );

  TProfile modlkvst1( "modlkvst1",
		      "triplet-MOD links vs time;time [s];triplets with MOD links / s",
		      300, 0, 300, -0.5, 1.5 );
  TProfile modlkvst3( "modlkvst3",
		      "triplet-MOD links vs time;time [s];triplets with MOD links / 10s",
		      150, 0, 1500, -0.5, 1.5 );
  TProfile modlkvst5( "modlkvst5",
		      "triplet-MOD links vs time;time [s];triplets with MOD links / min",
		      1100, 0, 66000, -0.5, 1.5 );
  TProfile modlkvsev( "modlkvsev",
		      "triplet-MOD links vs events;events;triplets with MOD links / 1000",
		      1100, 0, 1.1E6, -0.5, 1.5 );
  TProfile modlkvsev9( "modlkvsev9",
		       "triplet-MOD links vs events;events;triplets with MOD links / 1000",
		       42000, 0, 42E6, -0.5, 1.5 );
  TProfile modlkvsev1( "modlkvsev1",
		       "triplet-MOD links vs events;events;triplets with MOD links / 100",
		       100, 89e3, 99E3, -0.5, 1.5 );
  TProfile modlkvsev2( "modlkvsev2",
		       "triplet-MOD links vs events;events;triplets with MOD links / 100",
		       100, 570e3, 580E3, -0.5, 1.5 );

  TH1I ntrimodHisto( "ntrimod", "triplet - MOD links;triplet - MOD links;events",
		     11, -0.5, 10.5 );

  // DUT clusters:

  TH1I trixcHisto( "trixc", "triplets x at DUT;track x at DUT [mm];triplets",
		  240, -12, 12 );
  TH1I triycHisto( "triyc", "triplets y at DUT;track y at DUT [mm];triplets",
		  120, -6, 6 );
  TH2I * trixycHisto = new
    TH2I( "trixyc", "triplets x-y at DUT;track x at DUT [mm];track y at DUT [mm];triplets",
	  240, -12, 12, 120, -6, 6 );

  TH1I trixcmodHisto( "trixcmod", "MOD-triplets x at DUT;track x at DUT [mm];MOD-triplets",
		  240, -12, 12 );
  TH1I triycmodHisto( "triycmod", "MOD-triplets y at DUT;track y at DUT [mm];MOD-triplets",
		  120, -6, 6 );
  TH2I * trixycmodHisto = new
    TH2I( "trixycmod", "MOD-triplets x-y at DUT;track x at DUT [mm];track y at DUT [mm];MOD-triplets",
	  240, -12, 12, 120, -6, 6 );

  TH1I ttdmin1Histo( "ttdmin1",
		     "telescope triplets isolation;triplet min #Delta_{xy} [mm];triplet pairs",
		     100, 0, 1 );
  TH1I ttdmin2Histo( "ttdmin2",
		     "telescope triplets isolation;triplet min #Delta_{xy} [mm];triplet pairs",
		     150, 0, 15 );

  // driplets - triplets:

  TH1I hsixdx( "sixdx", "six dx;dx [mm];triplet-driplet pairs", 200, -0.2*f, 0.2*f );
  TH1I hsixdy( "sixdy", "six dy;dy [mm];triplet-driplet pairs", 200, -0.2*f, 0.2*f );
  TH1I hsixdxc( "sixdxc", "six dx;dx [mm];triplet-driplet pairs", 200, -0.2*f, 0.2*f );
  TH1I hsixdyc( "sixdyc", "six dy;dy [mm];triplet-driplet pairs", 200, -0.2*f, 0.2*f );

  TProfile sixdxvsx( "sixdxvsx",
		     "six #Deltax vs x;x [mm];<driplet - triplet #Deltax [mm]",
		     220, -11, 11, -0.1, 0.1 );
  TProfile sixmadxvsx( "sixmadxvsx",
		       "six MAD x vs x;x [mm];driplet - triplet MAD #Deltax [mm]",
		       220, -11, 11, 0, 0.1 );
  TProfile sixmadxvsy( "sixmadxvsy",
		       "six MAD x vs y;y [mm];driplet - triplet MAD #Deltax [mm]",
		       110, -5.5, 5.5, 0, 0.1 );
  TProfile sixmadxvstx( "sixmadxvstx",
			"six MAD x vs x;triplet #theta_{x} [rad];driplet - triplet MAD #Deltax [mm]",
			80, -0.002, 0.002, 0, 0.1 );
  TProfile sixmadxvsdtx( "sixmadxvsdtx",
			 "six MAD x vs x;driplet-triplet #Delta#theta_{x} [rad];driplet - triplet MAD #Deltax [mm]",
			 80, -0.002, 0.002, 0, 0.1 );
  TProfile sixdxvsy( "sixdxvsy",
		     "six #Deltax vs y;y [mm];<driplet - triplet #Deltax> [mm]",
		     100, -5, 5, -0.1, 0.1 );
  TProfile sixdxvstx( "sixdxvstx",
		      "six #Deltax vs slope x;slope x [rad];<driplet - triplet #Deltax> [mm]",
		      100, -0.002, 0.002, -0.1, 0.1 );
  TProfile sixdxvsdtx( "sixdxvsdtx",
		       "six #Deltax vs #Delta slope x;#Delta slope x [rad];<driplet - triplet #Deltax> [mm]",
		       100, -0.002, 0.002, -0.1, 0.1 );
  TProfile sixdxvst3( "sixdxvst3",
		      "sixplet dx vs time;time [s];<sixplet #Deltax> [mm] / 10s",
		      300, 0, 3000, -0.05, 0.05 );
  TProfile sixdxvst5( "sixdxvst5",
		      "sixplet dx vs time;time [s];<sixplet #Deltax> [mm] / min",
		      1100, 0, 66000, -0.05, 0.05 );

  TProfile sixdyvsx( "sixdyvsx",
		     "six #Deltay vs x;x [mm];<driplet - triplet #Deltay> [mm]",
		     200, -10, 10, -0.1, 0.1 );
  TProfile sixdyvsy( "sixdyvsy",
		     "six #Deltay vs y;y [mm];<driplet - triplet #Deltay [mm]",
		     110, -5.5, 5.5, -0.1, 0.1 );
  TProfile sixdyvsty( "sixdyvsty",
		      "six #Deltay vs slope y;slope y [rad];<driplet - triplet #Deltay> [mm]",
		      100, -0.002, 0.002, -0.1, 0.1 );
  TProfile sixdyvsdty( "sixdyvsdty",
		       "six #Deltay vs #Delta slope y;#Delta slope y [rad];<driplet - triplet #Deltay> [mm]",
		       100, -0.002, 0.002, -0.1, 0.1 );
  TProfile sixdyvst3( "sixdyvst3",
		      "sixplet dy vs time;time [s];<sixplet #Deltay> [mm] / 10s",
		      300, 0, 3000, -0.05, 0.05 );
  TProfile sixdyvst5( "sixdyvst5",
		      "sixplet dy vs time;time [s];<sixplet #Deltay> [mm] / min",
		      1100, 0, 66000, -0.05, 0.05 );
  TProfile sixmadyvsx( "sixmadyvsx",
		       "six MAD y vs x;x [mm];driplet - triplet MAD #Deltay [mm]",
		       220, -11, 11, 0, 0.1 );
  TProfile sixmadyvsy( "sixmadyvsy",
		       "six MAD y vs y;y [mm];driplet - triplet MAD #Deltay [mm]",
		       110, -5.5, 5.5, 0, 0.1 );
  TProfile sixmadyvsty( "sixmadyvsty",
			"six MAD y vs #theta_{y};triplet #theta_{y} [rad];driplet - triplet MAD #Deltay [mm]",
			80, -0.002, 0.002, 0, 0.1 );
  TProfile sixmadyvsdty( "sixmadyvsdty",
			 "six MAD y vs #Delta#theta_{y};driplet-triplet #Delta#theta_{y} [rad];driplet - triplet MAD #Deltay [mm]",
			 80, -0.002, 0.002, 0, 0.1 );

  TProfile sixdtvsxav( "sixdtvsxav",
		       "driplet - triplet kink_{xy} vs x;x_{avg} [mm];<sqrt(#Delta#theta_{x}^{2}+#Delta#theta_{y}^{2})> [rad]",
		       240, -12, 12, 0, 0.1 );
  TProfile sixdtvsyav( "sixdtvsyav",
		       "driplet - triplet kink_{xy} vs y;y_{avg} [mm];<sqrt(#Delta#theta_{x}^{2}+#Delta#theta_{y}^{2})> [rad]",
		       120, -6, 6, 0, 0.1 );
  TProfile2D * sixdtvsxyav = new
    TProfile2D( "sixdtvsxyav",
		"driplet - triplet kink_{xy} vs x-y;x_{avg} [mm];y_{avg} [mm];<sqrt(#Delta#theta_{x}^{2}+#Delta#theta_{y}^{2})> [rad]",
		120, -12, 12, 60, -6, 6, 0, 0.1 );

  TH1I hsixdtx( "sixdtx",
		"driplet slope x - triplet slope x;driplet slope x - triplet slope x;driplet-triplet pairs",
		100, -0.005*f, 0.005*f );
  TH1I hsixdty( "sixdty",
		"driplet slope y - triplet slope y;driplet slope y - triplet slope y;driplet-triplet pairs",
		100, -0.005*f, 0.005*f );

  TProfile sixdtvsx( "sixdtvsx",
		     "driplet - triplet kink_{xy} vs x;x_{mid} [mm];<sqrt(#Delta#theta_{x}^{2}+#Delta#theta_{y}^{2})> [rad]",
		     110, -11, 11, 0, 0.1 );
  TProfile2D * sixdtvsxy = new
    TProfile2D( "sixdtvsxy",
		"driplet - triplet kink_{xy} vs x-y;x_{DUT} [mm];y_{DUT} [mm];<sqrt(#Delta#theta_{x}^{2}+#Delta#theta_{y}^{2})> [rad]",
		110, -11, 11, 55, -5.5, 5.5, 0, 0.1 );

  TProfile sixdtvsxm( "sixdtvsxm",
		      "driplet - triplet kink_{xy} vs xmod;track x mod 100 [#mum];<sqrt(#Delta#theta_{x}^{2}+#Delta#theta_{y}^{2})> [rad]",
		      50, 0, 100, 0, 0.1 );
  TProfile sixdtvsym( "sixdtvsym",
		      "driplet - triplet kink_{xy} vs ymod;track y mod 100 [#mum];<sqrt(#Delta#theta_{x}^{2}+#Delta#theta_{y}^{2})> [rad]",
		      50, 0, 100, 0, 0.1 );
  TProfile2D * sixdtvsxmym = new // bump bonds ?
    TProfile2D( "sixdtvsxmym",
		"driplet - triplet kink_{xy} vs xmod ymod;track x mod 100 [#mum];track y mod 100 [#mum];<sqrt(#Delta#theta_{x}^{2}+#Delta#theta_{y}^{2})> [rad]",
		50, 0, 100, 50, 0, 100, 0, 0.1 );

  TH2I * sixxyHisto = new
    TH2I( "sixxy", "sixplets at z DUT;x [mm];y [mm];sixplets",
	  240, -12, 12, 120, -6, 6 );

  TH1I nsixHisto( "nsix", "sixplets;sixplets;events", 51, -0.5, 50.5 );

  TH2I * dutxyHisto = new
    TH2I( "dutxy", "tracks at DUT;x [mm];y [mm];tracks",
	  240, -12, 12, 120, -6, 6 );

  // DUT vs triplets:

  TH1I dutdxaHisto( "dutdxa",
		    "DUT - track dx;DUT cluster - track #Deltax [mm];DUT clusters",
		    1000, -50, 50 );
  TH1I dutsxaHisto( "dutsxa",
		    "DUT + track sx;DUT cluster + track #Deltax [mm];DUT clusters",
		    1000, -50, 50 );
  TH1I dutdyaHisto( "dutdya",
		    "DUT - track dy;DUT cluster - track #Deltay [mm];DUT clusters",
		    500, -10, 10 );
  TH1I dutsyaHisto( "dutsya",
		    "DUT + track sy;DUT cluster + track #Deltay [mm];DUT clusters",
		    500, -10, 10 );

  TH1I dutdxHisto( "dutdx",
		   "DUT - track dx;DUT cluster - track #Deltax [mm];DUT clusters",
		   500, -2.5, 2.5 );
  TH1I dutdyHisto( "dutdy",
		   "DUT - track dy;DUT cluster - track #Deltay [mm];DUT clusters",
		   500, -0.25, 0.25 );

  TH2I * dutxxHisto = new
    TH2I( "dutxx", "tracks vs DUT in x;track x [mm];DUT cluster x [mm];track-cluster pairs",
	  500, -5, 5, 500, -5, 5 );

  TH1I dutdxcHisto( "dutdxc",
		    "DUT - track dx;DUT cluster - track #Deltax [mm];DUT clusters",
		    500, -2.5, 2.5 );
  TH1F dutresxHisto( "dutresx",
		    "DUT - track dx;DUT cluster - track #Deltax [mm];DUT clusters",
		    500, -0.1, 0.1 );
  TH1I dutdxccHisto( "dutdxcc",
		    "DUT - track dx crack;DUT cluster - track #Deltax [mm];DUT crack clusters",
		    500, -0.25, 0.25 );
  TH1I dutdxcmHisto( "dutdxcm",
		    "DUT - track dx mid pix;DUT cluster - track #Deltax [mm];DUT mid pix clusters",
		    500, -0.25, 0.25 );
  TH1I lindxcHisto( "lindxc",
		    "Lin - track dx;Lin cluster - track #Deltax [mm];Lin clusters",
		    200, -0.25, 0.25 );
  TH1I dutdxc1Histo( "dutdxc1",
		     "DUT - track dx;DUT cluster - track #Deltax [mm];DUT 1-px clusters",
		     200, -0.25, 0.25 );
  TH1I dutdxc2Histo( "dutdxc2",
		     "DUT - track dx;DUT cluster - track #Deltax [mm];DUT 2-px clusters",
		     200, -0.25, 0.25 );
  TH1I dutdxcqHisto( "dutdxcq",
		     "DUT - track dx Landau peak;DUT cluster - track #Deltax [mm];Landau peak DUT clusters",
		     500, -0.25, 0.25 );

  double limx = 0.05 * ( 1 + fifty ) * ( 1 + (1-rot90)*(1-fifty) );
  if( fabs(DUTturn) > 44 )
    limx = 0.5;
  if( fabs(DUTturn) > 44 )
    limx = 0.5;

  TProfile dutdxvsx( "dutdxvsx",
		     "DUT #Deltax vs x;x track [mm];<cluster - track #Deltax> [mm]",
		     200, -10, 10, -limx, limx );
  TProfile dutdxvsy( "dutdxvsy",
		     "DUT #Deltax vs y;y track [mm];<cluster - track #Deltax> [mm]",
		     160, -8, 8, -limx, limx );
  TProfile dutdxvstx( "dutdxvstx",
		      "DUT #Deltax vs #theta_{x};x track slope [rad];<cluster - track #Deltax> [mm]",
		      80, -0.002, 0.002, -limx, limx );
  TProfile dutdxvstx0( "dutdxvstx0",
		      "DUT #Deltax vs #theta_{x}, x < 0;x track slope [rad];<cluster - track #Deltax> [mm]",
		      80, -0.002, 0.002, -limx, limx );
  TProfile dutdxvstx1( "dutdxvstx1",
		      "DUT #Deltax vs #theta_{x}, x > 0;x track slope [rad];<cluster - track #Deltax> [mm]",
		      80, -0.002, 0.002, -limx, limx );
  TProfile dutdxvsxm( "dutdxvsxm",
		      "DUT #Deltax vs xmod;x track mod 100 [#mum];<cluster - track #Deltax> [mm]",
		      50, 0, 100, -limx, limx );
  TProfile2D * dutdxvsxmym = new TProfile2D( "dutdxvsxmym",
		"DUT #Deltax vs xmod ymod;x track mod 100 [#mum];y track mod 100",
 		 50, 0, 100, 50, 0, 100, -0.1, 0.1 );
  TProfile dutdxvst2( "dutdxvst2",
		      "DUT #Deltax vs time;time [s];<DUT #Deltax> [mm/5s]",
		      200, 0, 1000, -limx, limx );
  TProfile dutdxvst5( "dutdxvst5",
		      "DUT #Deltax vs time;time [s];<DUT #Deltax> [mm]",
		      1100, 0, 66000, -limx, limx );

  TProfile dutmadxvsx( "dutmadxvsx",
		       "DUT MAD(#Deltax) vs x;x track [mm];MAD(#Deltax) [mm]",
		       200, -10, 10, 0, limx );
  TProfile dutmadxvsxm( "dutmadxvsxm",
			"DUT MAD(#Deltax) vs xmod;x track mod 100 [#mum];MAD(#Deltax) [mm]",
			50, 0, 100, 0, limx );
  TProfile dutmadxvsxm5( "dutmadxvsxm5",
			"DUT MAD(#Deltax) vs xmod;x track mod 50 [#mum];MAD(#Deltax) [mm]",
			50, 0, 50, 0, limx );
  TProfile dutmadxvstx( "dutmadxvstx",
			"DUT MAD(#Deltax) vs #theta_{x};x track slope [mrad];MAD(#Deltax) [mm]",
			80, -2, 2, 0, limx );
  TProfile dutmadxvsq( "dutmadxvsq",
		       "DUT MAD(#Deltax) vs Q;cluster signal [ToT];MAD(#Deltax) [mm]",
		       80, 0, 80, 0, limx );

  TH2I * dutyyHisto = new
    TH2I( "dutyy", "tracks vs DUT in y;track y [mm];DUT cluster y [mm];track-cluster pairs",
	  500, -5, 5, 500, -5, 5 );

  TH1I dutdy0Histo( "dutdy0",
		    "DUT - track dy;DUT cluster - track #Deltay [mm];DUT clusters",
		    500, -5, 5 ); // binning like dutyy
  TH1I dutdycHisto( "dutdyc",
		    "DUT - track dy;DUT cluster - track #Deltay [mm];DUT clusters",
		    500, -0.5, 0.5 );
  TH1F dutresyHisto( "dutresy",
		    "DUT - track dy;DUT cluster - track #Deltay [mm];DUT clusters",
		    500, -0.1, 0.1 );
  TProfile dutmadyvsq( "dutmadyvsq",
		       "DUT MAD(#Deltay) vs Q;cluster signal [ToT];MAD(#Deltay) [mm]",
		       80, 0, 80, 0, 0.1 );
  TH1I dutdycqHisto( "dutdycq",
		     "DUT - track dy Landau peak;DUT cluster - track #Deltay [mm];Landau peak DUT clusters",
		     500, -0.5, 0.5 );

  TProfile dutdyvsx( "dutdyvsx",
		     "DUT #Deltay vs x;x track [mm];<cluster - track #Deltay> [mm]",
		     200, -10, 10, -0.1, 0.1 );
  TProfile dutdyvsy( "dutdyvsy",
		     "DUT #Deltay vs y;y track [mm];<cluster - track #Deltay> [mm]",
		     160, -8, 8, -0.1, 0.1 );
  TProfile dutdyvsty( "dutdyvsty",
		      "DUT #Deltay vs #theta_{y};y track slope [rad];<cluster - track #Deltay> [mm]",
		      80, -0.002, 0.002, -0.1, 0.1 );
  TProfile dutdyvsxm( "dutdyvsxm",
		      "DUT #Deltay vs xmod;x track mod 100 [#mum];<cluster - track #Deltay> [mm]",
		      50, 0, 100, -0.1, 0.1 );
  TProfile dutdyvsym( "dutdyvsym",
		      "DUT #Deltay vs ymod;y track mod 100 [#mum];<cluster - track #Deltay> [mm]",
		      50, 0, 100, -0.1, 0.1 );
  TProfile2D * dutdyvsxmym = new
    TProfile2D( "dutdyvsxmym",
		"DUT #Deltay vs xmod ymod;x track mod 100 [#mum];y track mod 100 [#mum];<#Deltay> [mm]",
		50, 0, 100, 50, 0, 100, -0.1, 0.1 );
  TProfile dutdyvst2( "dutdyvst2",
		      "DUT #Deltay vs time;time [s];<DUT #Deltay> [mm/5s]",
		      200, 0, 1000, -0.1, 0.1 );
  TProfile dutdyvst5( "dutdyvst5",
		      "DUT #Deltay vs time;time [s];<DUT #Deltay> [mm]",
		      1100, 0, 66000, -0.1, 0.1 );

  TProfile dutmadyvsx( "dutmadyvsx",
		       "DUT MAD(#Deltay) vs x;x track [mm];MAD(#Deltay) [mm]",
		       200, -10, 10, 0, 0.1 );
  TProfile dutmadyvsy( "dutmadyvsy",
		       "DUT MAD(#Deltay) vs y;y track [mm];MAD(#Deltay) [mm]",
		       160, -8, 8, 0, 0.1 );
  TProfile dutmadyvstx( "dutmadyvstx",
			"DUT MAD(#Deltay) vs #theta_{x};x track slope [rad];MAD(#Deltay) [mm]",
			80, -0.002, 0.002, 0, 0.1 );
  TProfile dutmadyvsty( "dutmadyvsty",
			"DUT MAD(#Deltay) vs #theta_{y};y track slope [rad];MAD(#Deltay) [mm]",
			80, -0.002, 0.002, 0, 0.1 );
  TProfile dutmadyvsxm( "dutmadyvsxm",
			"DUT MAD(#Deltay) vs xmod;x track mod 100 [#mum];MAD(#Deltay) [mm]",
			50, 0, 100, 0, 0.1 );
  TProfile dutmadyvsym( "dutmadyvsym",
			"DUT MAD(#Deltay) vs ymod;y track mod 100 [#mum];MAD(#Deltay) [mm]",
			50, 0, 0.1, 0, 0.1 );
  TProfile dutmadyvst( "dutmadyvst",
		       "DUT MAD(#Deltay) vs time;time [s];MAD(#Deltay) [mm]",
		       100, 0, 1000, 0, 0.1 );

  TH1I dutncolHisto( "dutncol",
		     "DUT linked cluster cols;DUT cluster size [columns];linked DUT clusters",
		     40, 0.5, 40.5 );
  TH1I dutncolfHisto( "dutncolf",
		      "DUT fiducial linked cluster cols;DUT cluster size [columns];fiducial linked DUT clusters",
		      40, 0.5, 40.5 );
  TProfile dutncolvsx( "dutncolvsx",
		       "DUT cluster length vs x;x track [mm];<cluster length> [columns]",
		       200, -10, 10, 0, 99 );

  TH2I * trixy4shrHisto = new
    TH2I( "trixy4shr", "short cluster triplets x-y;x [mm];y [mm];short cluster triplets",
	  240, -12, 12, 120, -6, 6 );
  TH1I trix4shrHisto( "trix4shr",
		    "short cluster track at DUT x;track x at DUT [mm];short cluster tracks",
		    240, -12, 12 );
  TH1I triy4shrHisto( "triy4shr",
		    "short cluster track at DUT y;track y at DUT [mm];short cluster tracks",
		    120, -6, 6 );

  TH1I minbcHisto( "minbc",
		   "DUT cluster min BC;DUT cluster min BC;DUT clusters on tracks",
		   32, -0.5, 31.5 );
  TH1I maxbcHisto( "maxbc",
		   "DUT cluster max BC;DUT cluster max BC;DUT clusters on tracks",
		   32, -0.5, 31.5 );
  TH1I nbcHisto( "linnbc",
		    "DUT cluster nBC;DUT cluster nBC;DUT clusters on tracks",
		    30, 0.5, 30.5 );

  TH1I linqHisto( "linq",
		  "LIN linked clusters;LIN cluster signal [ToT];linked LIN clusters",
		  100, 0, 80*qscale );
  TH1I linq0Histo( "linq0",
		   "LIN linked clusters;LIN normal cluster signal [ToT];linked LIN clusters",
		   500, 0, 80*qscale );
  TProfile linqxvsx( "linqxvsx",
		    "LIN cluster signal vs x;x track [mm];LIN <cluster signal> [ToT]",
		    200, -10, 10, 0, qxmax );
  TProfile linqxvsy( "linqxvsy",
		    "LIN cluster signal vs y;y track [mm];LIN <cluster signal> [ToT]",
		    500, -5, 5, 0, qxmax );
  TProfile2D * linqxvsxy = new
    TProfile2D( "linqxvsxy",
		"LIN cluster signal vs xy;x track [mm];y track [mm];LIN <cluster signal> [ToT]",
		200, -10, 10, 100, -5, 5, 0, qxmax );
  TProfile linqxvsr( "linqxvsr",
		    "LIN cluster signal vs Rbeam;R track from beam [mm];LIN <cluster signal> [ToT]",
		    100, 0, 10, 0, qxmax );

  TH1I linnpxHisto( "linnpx",
		    "LIN linked cluster size;LIN cluster size [pixels];linked LIN clusters",
		    20, 0.5, 20.5 );
  TH1I linncolHisto( "linncol",
		     "LIN linked cluster cols;LIN cluster size [columns];linked LIN clusters",
		     20, 0.5, 20.5 );
  TH1I linnrowHisto( "linnrow",
		     "LIN linked cluster rows;LIN cluster size [rows];linked LIN clusters",
		     20, 0.5, 20.5 );

  TH1I linnrow1Histo( "linnrow1",
		      "LIN linked 1-col cluster rows;LIN cluster size [rows];linked LIN 1-col clusters",
		      20, 0.5, 20.5 );
  TH1I linnrow1eveHisto( "linnrow1eve",
			 "LIN linked even seed row 1-col cluster rows;LIN cluster size [rows];linked LIN even seed row 1-col clusters",
			 20, 0.5, 20.5 );
  TH1I linnrow1oddHisto( "linnrow1odd",
			 "LIN linked odd seed row 1-col cluster rows;LIN cluster size [rows];linked LIN odd seed row 1-col clusters",
			 20, 0.5, 20.5 );

  TH1I linnpxfHisto( "linnpxf",
		    "LIN linked 1-BC cluster size;LIN cluster size [pixels];linked 1-BC LIN clusters",
		     20, 0.5, 20.5 );
  TH1I linncolfHisto( "linncolf",
		      "LIN linked 1-BC cluster cols;LIN cluster size [columns];linked 1-BC LIN clusters",
		      20, 0.5, 20.5 );
  TH1I linnrowfHisto( "linnrowf",
		      "LIN linked 1-BC cluster rows;LIN cluster size [rows];linked 2-BC LIN clusters",
		      20, 0.5, 20.5 );

  TH1I linnpxffHisto( "linnpxff",
		      "LIN linked 2-BC cluster size;LIN cluster size [pixels];linked 2-BC LIN clusters",
		      20, 0.5, 20.5 );
  TH1I linncolffHisto( "linncolff",
		       "LIN linked 2-BC cluster cols;LIN cluster size [columns];linked 2-BC LIN clusters",
		       20, 0.5, 20.5 );
  TH1I linnrowffHisto( "linnrowff",
		       "LIN linked 2-BC cluster rows;LIN cluster size [rows];linked 2-BC LIN clusters",
		       20, 0.5, 20.5 );

  TProfile2D * linnpxvsxmym = new
    TProfile2D( "linnpxvsxmym",
		"LIN cluster size vs xmod ymod;x track mod 100 [#mum];y track mod 100 [#mum];LIN <cluster size> [pixels]",
		50, 0, 100, 50, 0, 100, 0, 20 );
  linnpxvsxmym->SetMinimum(1);

  TProfile linncolvsxm( "linncolvsxm",
			"LIN cluster size vs xmod;x track mod 100 [#mum];LIN <cluster size> [columns]",
			50, 0, 100, 0, 20 );
  linncolvsxm.SetMinimum(1);

  TProfile linncolvsym( "linncolvsym",
			"LIN cluster size vs ymod;y track mod 100 [#mum];LIN <cluster size> [columns]",
			50, 0, 100, 0, 20 );
  linncolvsym.SetMinimum(1);

  TProfile2D * linnrowvsxmym = new
    TProfile2D( "linnrowvsxmym",
		"LIN cluster size vs xmod ymod;x track mod 100 [#mum];y track mod 100 [#mum];LIN <cluster size> [rows]",
		50, 0, 100, 50, 0, 100, 0, 20 );
  linnrowvsxmym->SetMinimum(1);

  TProfile linnrowvsxm( "linnrowvsxm",
			"LIN cluster size vs xmod;x track mod 100 [#mum];LIN <cluster size> [rows]",
			50, 0, 100, 0, 20 );
  linnrowvsxm.SetMinimum(1);

  TProfile linnrowvsxm5( "linnrowvsxm5",
			 "LIN cluster size vs xmod;x track mod 100 [#mum];LIN <cluster size> [rows]",
			 50, 0, 50, 0, 20 );
  linnrowvsxm5.SetMinimum(1);

  TProfile linnrowvsym( "linnrowvsym",
			"Lin cluster size vs ymod;y track mod 100 [#mum];Lin <cluster size> [rows]",
			50, 0, 100, 0, 20 );
  linnrowvsym.SetMinimum(1);

  TProfile2D * linqxvsxmym = new
    TProfile2D( "linqxvsxmym",
		"LIN cluster signal vs xmod ymod;x track mod 100 [#mum];y track mod 100 [#mum];LIN <cluster signal> [ToT]",
		50, 0, 100, 50, 0, 100, 0, qxmax );
  TProfile linqxvsxm( "linqxvsxm",
		     "LIN cluster signal vs xmod;x track mod 100 [#mum];LIN <cluster signal> [ToT]",
		     50, 0, 100, 0, qxmax );
  TProfile linqxvsxm2( "linqxvsxm2",
		     "LIN cluster signal vs xmod2;x track mod 200 [#mum];LIN <cluster signal> [ToT]",
		     100, 0, 200, 0, qxmax );
  TProfile linqxvsxm5( "linqxvsxm5",
		      "LIN cluster signal vs xmod;x track mod 50 [#mum];LIN <cluster signal> [ToT]",
		      50, 0, 50, 0, qxmax );
  TProfile linqxvsym( "linqxvsym",
		     "LIN cluster signal vs ymod;y track mod 100 [#mum];LIN <cluster signal> [ToT]",
		     50, 0, 100, 0, qxmax );

  TH2I * trixyclkHisto = new
    TH2I( "trixyclk", "linked triplets x-y;x [mm];y [mm];linked triplets",
	  240, -12, 12, 120, -6, 6 );
  TH1I trixclkHisto( "trixclk",
		    "linked track at DUT x;track x at DUT [mm];linked tracks",
		    240, -12, 12 );
  TH1I triyclkHisto( "triyclk",
		    "linked track at DUT y;track y at DUT [mm];linked tracks",
		    120, -6, 6 );

  TH2I * dutlkxmymHisto = new
    TH2I( "dutlkxmym",
	  "linked tracks at DUT xmod ymod;x track mod 100 [#mum];y track mod 100 [#mum];linked tracks",
	  50, 0, 100, 50, 0, 100 );

  int nbc = 400;
  int nbr = 192;
  if( !fifty ) {
    nbc = 200; // sensor
    nbr = 384;
  }

  TH1I dutlkcolHisto( "dutlkcol",
		      "DUT linked col;DUT linked col;linked DUT clusters",
		      nbc, 0, nbc );
  TH1I dutlkrowHisto( "dutlkrow",
		      "DUT linked row;DUT linked row;linked DUT clusters",
		      nbr, 0, nbr );

  TH1I dutpxcolHisto( "dutpxcol",
		      "DUT pixel column;DUT pixel column;linked pixels",
		      nbc, 0, nbc );
  TH1I dutpxrowHisto( "dutpxrow",
		      "DUT pixel row;DUT pixel row;linked pixels",
		      nbr, 0, nbr );
  TH1I dutpxcol1Histo( "dutpxcol1",
		       "DUT 1st pixel col;DUT pixel col;linked 1st pixels",
		       nbc, 0, nbc );
  TH1I dutpxrow1Histo( "dutpxrow1",
		       "DUT 1st pixel row;DUT pixel row;linked 1st pixels",
		       nbr, 0, nbr );
  TH1I dutpxcol11Histo( "dutpxcol11",
			"DUT 1-row 1st pixel col;DUT pixel col;linked 1-row 1st pixels",
			nbc, 0, nbc );
  TH1I dutpxrow11Histo( "dutpxrow11",
			"DUT 1-col 1st pixel row;DUT pixel row;linked 1-col 1st pixels",
			nbr, 0, nbr );
  TH1I dutpxrow111Histo( "dutpxrow111",
			 "DUT 1-BC 1-col 1st pixel row;DUT pixel row;linked 1-BC 1-col 1st pixels",
			 nbr, 0, nbr );
  TH1I dutpxrow1112Histo( "dutpxrow1112",
			  "DUT 1-BC 1-col 1st pixel row;DUT pixel row;linked 1-BC 1-col 2-row 1st pixels",
			  nbr, 0, nbr );
  TH1I dutpxcol111eveHisto( "dutpxcol111eve",
			    "DUT 1-BC 1-col 1st row even column;DUT pixel column;linked 1-BC 1-col 1st row even pixels",
			    nbc, 0, nbc );
  TH1I dutpxcol111oddHisto( "dutpxcol111odd",
			    "DUT 1-BC 1-col 1st row odd column;DUT pixel column;linked 1-BC 1-col 1st row odd pixels",
			    nbc, 0, nbc );

  TH1I dutpxqHisto( "dutpxq",
		     "DUT pixel signal;DUT pixel signal [ToT];linked pixels",
		     qnbins, 0, qvalmax );
  TH1I dutpxqbeamHisto( "dutpxqbeam",
		     "DUT beam pixel signal;DUT pixel signal [ToT];linked beam pixels",
		     qnbins, 0, qvalmax );
  TH1I dutpxq1Histo( "dutpxq1",
		     "DUT 1-col pixel signal;DUT pixel signal [ToT];linked 1-col pixels",
		     qnbins, 0, qvalmax );
  TH1I dutpxq2Histo( "dutpxq2",
		     "DUT 2-col pixel signal;DUT pixel signal [ToT];linked 2-col pixels",
		     qnbins, 0, qvalmax );

  TH1I dutpxqeveHisto( "dutpxqeve",
		       "DUT even pixel signal;DUT pixel signal [ToT];linked even pixels",
		       qnbins, 0, qvalmax );
  TH1I dutpxqoddHisto( "dutpxqodd",
		       "DUT odd pixel signal;DUT pixel signal [ToT];linked odd pixels",
		       qnbins, 0, qvalmax );

  TH1I dutpxdcolHisto( "dutpxdcol",
		       "DUT pixel dcol;DUT pixel dcol;linked subsequent pixels",
		       11, -5.5, 5.5 );
  TH1I dutpxdrowHisto( "dutpxdrow",
		       "DUT pixel drow;DUT pixel drow;linked subsequent pixels",
		       11, -5.5, 5.5 );
  TH1I dutpxdfrmHisto( "dutpxdfrm",
		       "DUT pixel dfrm;DUT pixel dfrm;linked subsequent pixels",
		       11, -5.5, 5.5 );
  TH2I dutpxqqHisto( "dutpxqq",
		     "DUT pixel correlation;DUT pixel signal [ToT];DUT pixel signal [ToT];linked pixels",
		     qnbins, 0, qvalmax, qnbins, 0, qvalmax );
  TH2I dutpxqq0Histo( "dutpxqq0",
		      "DUT even-odd row correlation;DUT pixel signal [ToT];DUT pixel signal [ToT];linked even-odd row pixels",
		     qnbins, 0, qvalmax, qnbins, 0, qvalmax );
  TH2I dutpxqq1Histo( "dutpxqq1",
		      "DUT odd-even row correlation;DUT pixel signal [ToT];DUT pixel signal [ToT];linked odd-even row pixels",
		     qnbins, 0, qvalmax, qnbins, 0, qvalmax );
  TH1I dutpxq1eveHisto( "dutpxq1eve",
			"DUT 1-col even row 1st pixel signal;DUT pixel signal [ToT];linked 1-col even row 1st pixels",
			qnbins, 0, qvalmax );
  TH1I dutpxq1oddHisto( "dutpxq1odd",
			"DUT 1-col odd row 1st pixel signal;DUT pixel signal [ToT];linked 1-col odd row 1st pixels",
			qnbins, 0, qvalmax );
  TH1I dutpxq2eveHisto( "dutpxq2eve",
			"DUT 1-col even row 2nd pixel signal;DUT pixel signal [ToT];linked 1-col even row 2nd pixels",
			qnbins, 0, qvalmax );
  TH1I dutpxq2oddHisto( "dutpxq2odd",
			"DUT 1-col odd row 2nd pixel signal;DUT pixel signal [ToT];linked 1-col odd row 2nd pixels",
			qnbins, 0, qvalmax );

  TH1I linpxqmaxHisto( "linpxqmax",
		     "Lin max pixel signal;Lin max pixel signal [ToT];Lin linked pixels",
			qnbins, 0, qvalmax );
  TH1I linpxqmax2Histo( "linpxqmax2",
			"Lin max pixel signal 2-pix clusters;Lin max pixel signal [ToT];Lin linked pixels in 2-pix clusters",
			qnbins, 0, qvalmax );
  TH1I linpxq2ndHisto( "linpxq2nd",
		     "Lin 2nd pixel signal;Lin 2nd pixel signal [ToT];Lin linked pixels",
			qnbins, 0, qvalmax );
  TProfile lintwvsq( "lintwvsq",
		     "Lin time walk vs q;Lin 2nd pixel signal [ToT];<time walk> [#DeltaBC]",
		     qnbins,  0, qvalmax );
  TProfile lintwvsq15( "lintwvsq15",
		     "Lin time walk vs q;Lin 2nd pixel signal [ToT];<time walk> [#DeltaBC]",
		     qnbins,  0, qvalmax );

  TH1I dutcolszHisto( "dutcolsz", "DUT column size;DUT column size [rows];columns inside linked clusters",
		      20, 0.5, 20.5 );
  TProfile dutcolszvsc( "dutcolszvsc", "DUT column size vs column;cluster column;DUT <column size> [rows]",
			40, 0.5, 40.5 );
  TH1I dutcolqHisto( "dutcolq", "DUT column signal;DUT column signal [ToT];columns inside linked clusters",
		     80, 0.5, 80.5 );
  TProfile dutcolqvsc( "dutcolqvsc", "DUT column signal vs column;cluster column;DUT <column signal> [ToT]",
		       40, 0.5, 40.5 );

  TH1I dutpdyHisto( "dutpdy",
		   "DUT row - track dy;DUT row - track #Deltay [mm];DUT rowq [ToT]",
		   100, -100, 100 );

  TProfile dutrowqvsdy( "dutrowqvsdy",
			"DUT row signal vs #Deltay;track-row #Deltay [#mum];DUT <row signal> [ToT]",
			100, -62.5, 62.5 );
  TProfile dutrowqvsdy0( "dutrowqvsdy0",
			 "DUT row signal vs #Deltay even;even track-row #Deltay [#mum];DUT <row signal> [ToT]",
			 100, -62.5, 62.5 );
  TProfile dutrowqvsdy1( "dutrowqvsdy1",
			 "DUT row signal vs #Deltay odd;odd track-row #Deltay [#mum];DUT <row signal> [ToT]",
			 100, -62.5, 62.5 );

  TProfile dutrowvsdy( "dutrowvsdy",
		       "DUT row vs #Deltay;track-row #Deltay [#mum];DUT <row>",
		       100, -62.5, 62.5 );
  TProfile dutrowvsdy0( "dutrowvsdy0",
			 "DUT row vs #Deltay even;even track-row #Deltay [#mum];DUT <row>",
			 100, -62.5, 62.5 );
  TProfile dutrowvsdy1( "dutrowvsdy1",
			 "DUT row vs #Deltay odd;odd track-row #Deltay [#mum];DUT <row>",
			 100, -62.5, 62.5 );

  TProfile2D * dutrowpqmap = new
    TProfile2D( "dutrowpqmap",
	  "row covariance;x track at DUT [mm];y track at DUT [mm];LIN < ToT_{row-1} #upoint ToT_{row} >",
		nbc, 0, nbc, nbr, 0, nbr );
  TProfile2D * dutrowppmap = new
    TProfile2D( "dutrowppmap",
	  "row variance;x track at DUT [mm];y track at DUT [mm];LIN < ToT_{row-1} #upoint ToT_{row-1} >",
	  nbc, 0, nbc, nbr, 0, nbr );
  TProfile2D * dutrowqqmap = new
    TProfile2D( "dutrowqqmap",
	  "row variance;x track at DUT [mm];y track at DUT [mm];LIN < ToT_{row} #upoint ToT_{row} >",
	  nbc, 0, nbc, nbr, 0, nbr );

  TProfile2D * dutcolpqmap = new
    TProfile2D( "dutcolpqmap",
	  "col covariance;x track at DUT [mm];y track at DUT [mm];LIN < ToT_{col-1} #upoint ToT_{col} >",
	  nbc, 0, nbc, nbr, 0, nbr );
  TProfile2D * dutcolppmap = new
    TProfile2D( "dutcolppmap",
	  "col variance;x track at DUT [mm];y track at DUT [mm];LIN < ToT_{col-1} #upoint ToT_{col-1} >",
	  nbc, 0, nbc, nbr, 0, nbr );
  TProfile2D * dutcolqqmap = new
    TProfile2D( "dutcolqqmap",
	  "col variance;x track at DUT [mm];y track at DUT [mm];LIN < ToT_{col} #upoint ToT_{col} >",
	  nbc, 0, nbc, nbr, 0, nbr );

  TH1I linqseedHisto( "linqseed",
		      "LIN seed row signal;LIN seed row signal [ToT];linked LIN clusters",
		      qnbins, 0, qvalmax );
  TH1I linqpairHisto( "linqpair",
		      "LIN pair row signal;LIN pair row signal [ToT];linked LIN clusters",
		      qnbins, 0, qvalmax );

  TProfile linqseedvsym( "linqseedvsym",
		      "LIN seed row signal vs ymod;y track mod 50 [#mum];1-col LIN <seed row signal> [ToT]",
		      50, 0, 50, -1, 99 );
  TProfile linqpairvsym( "linqpairvsym",
		      "LIN pair row signal vs ymod;y track mod 50 [#mum];1-col LIN <pair row signal> [ToT]",
		      50, 0, 50, -1, 99 );
  TH2I linqymHisto( "linqym",
		    "LIN cluster signal vs ymod;y track mod 50 [#mum];LIN <row signal> [ToT]",
		    25, 0, 50, 31, -15.5, 15.5 );

  TProfile linpqvsym( "linpqvsym",
		      "LIN signal covariance;y track mod 50 [#mum];LIN < ToT_{i-1} #upoint ToT_{i} >",
		      50, 0, 50, -1, 999 );
  TProfile linppvsym( "linppvsym",
		      "LIN signal variance;y track mod 50 [#mum];LIN < ToT_{i-1} #upoint ToT_{i-i} >",
		      50, 0, 50, -1, 999 );
  TProfile linqqvsym( "linqqvsym",
		      "LIN signal variance;y track mod 50 [#mum];LIN < ToT_{i} #upoint ToT_{i} >",
		      50, 0, 50, -1, 999 );

  TH1I linrowminHisto( "linrowmin",
			"LIN 1st pixel row;LIN 1st pixel row;linked 1-col clusters",
			nbr, 0, nbr );
  TH1I linrowmaxHisto( "linrowmax",
			"LIN lst pixel row;LIN lst pixel row;linked 1-col clusters",
			nbr, 0, nbr );
  TH1I linrowmin01Histo( "linrowmin01",
			"LIN 1st pixel row%2;LIN 1st pixel row mod 2;linked 1-col clusters",
			2, -0.5, 1.5 );
  TH1I linrowmax01Histo( "linrowmax01",
			"LIN lst pixel row%2;LIN lst pixel row mod 2;linked 1-col clusters",
			2, -0.5, 1.5 );
  TH1I linrowmin2Histo( "linrowmin2",
			"LIN 1st pixel row;LIN 1st pixel row;linked 1-col 2-row clusters",
			nbr, 0, nbr );
  TH1I linrowmax2Histo( "linrowmax2",
			"LIN lst pixel row;LIN lst pixel row;linked 1-col 2-row clusters",
			nbr, 0, nbr );

  TH1I dutetaHisto( "duteta",
		    "DUT eta;DUT eta;linked 1-col pixel pairs",
		    100, -1, 1 );
  TH1I duteta2Histo( "duteta2",
		    "DUT 2-row eta;DUT eta;linked 1-col 2-row pixel pairs",
		     100, -1, 1 );

  TH1I dutpxbclkHisto( "dutpxbclk", "DUT linked pixel BC;DUT pixel BC;linked pixels", 32, -0.5, 31.5 );
  TH1I dutpxbclk1pHisto( "dutpxbclk1p", "DUT linked 1-pix pixel BC;DUT pixel BC;linked pixels",
			 32, -0.5, 31.5 );
  TH1I dutpxbclk1cHisto( "dutpxbclk1c", "DUT linked 1-col pixel BC;DUT pixel BC;linked pixels",
			 32, -0.5, 31.5 );
  TH1I dutpxbclk1rHisto( "dutpxbclk1r", "DUT linked 1-row pixel BC;DUT pixel BC;linked pixels",
			 32, -0.5, 31.5 );

  TProfile dutpxbcvsq( "dutpxbcvsq",
		      "DUT pixel BC vs q;[ixel charge [ToT];<pixel BC>", 16, 0, 16 );
  TProfile dutpxbcvsx( "dutpxbcvsx",
		      "DUT pixel BC vs x;column;<pixel BC>", nbc, 0, nbc );
  TProfile dutpxbcvsy( "dutpxbcvsy",
		      "DUT pixel BC vs y;row;<pixel BC>", nbr, 0, nbr );
  TProfile2D * dutpxbcvsxy = new
    TProfile2D( "dutpxbcvsxy",
		"DUT pixel BC map;column;row;<pixel BC>", nbc, 0, nbc, nbr, 0, nbr );
  TProfile dutpxbcvst5( "dutpxbcvst5",
			"DUT pixel BC vs time;time [s];<pixel BC>", 1100, 0, 66000 );

  TH1I dutpxqbcHisto[32];
  for( int frm = 0; frm < 32; ++frm )
    dutpxqbcHisto[frm] =
      TH1I( Form( "dutpxqbc%i", frm ),
	    Form( "DUT pixel signal in BC %i;DUT pixel signal in BC %i [ToT];linked pixels", frm, frm ),
	    qnbins, 0, qvalmax );

  TH1I dutnlkHisto( "dutnlk", "DUT links per track;DUT links;in-time fiducial tracks",
		 11, -0.5, 10.5 );
  TH1I dutnlkshrHisto( "dutnlkshr", "DUT links per track, short cluster;DUT links;short cluster tracks",
		    11, -0.5, 10.5 );
  TH1I dutdxminHisto( "dutdxmin",
		      "nearest dx;nearest #Deltax [mm];efficient tracks",
		      400, -2, 2 );
  TH1I dutdyminHisto( "dutdymin",
		      "nearest dy;nearest #Deltay [mm];efficient tracks",
		      200, -1, 1 );
  TH1I dutqnearHisto( "dutqnear",
		      "DUT nearest linked cluster;DUT cluster signal [ToT];nearest linked DUT clusters",
		      80, 0, 80 );
  TH1I dutpxqnearHisto( "dutpxqnear",
			"DUT pixel signal on track;DUT pixel signal [ToT];DUT pixels on track",
			16, -0.5, 15.5 );
  TH1I dutpxbcnearHisto( "dutpxbcnear",
			 "DUT pixel BC on track;DUT pixel BC;DUT pixels on track",
			 32, -0.5, 31.5 );

  TH2I * dutxylkHisto = new
    TH2I( "dutxylk", "linked tracks at DUT;x [mm];y [mm];linked tracks",
	  480, -12, 12, 240, -6, 6 );

  TProfile effvsdmin( "effvsdmin",
		      "DUT efficiency vs triplet isolation at MOD;triplet isolation [mm];efficiency",
		      80, 0, 4, -1, 2 );

  TH2I * sixxylkHisto = new
    TH2I( "sixxylk",
	  "MOD-linked sixplets at z DUT;x [mm];y [mm];MOD-linked sixplets",
	  240, -12, 12, 120, -6, 6 );
  TH2I * sixxyeffHisto = new
    TH2I( "sixxyeff",
	  "MOD-linked sixplets with DUT;x [mm];y [mm];DUT-MOD-linked sixplets",
	  240, -12, 12, 120, -6, 6 );

  int nbx = 440;
  int nby = 240;
  if( rot90 ) { // HLL 25x100
    nbx = 880;
    nby = 120;
  }
  if( !fifty ) { // 100x25 straight
    nbx = 220;
    nby = 480;
  }
 
  float xcell_size_um = 50.0;
  float ycell_size_um = 50.0;
  if( !fifty )
  {
     xcell_size_um = 100.0;
     ycell_size_um = 25.0;
  }

  CellRelatedHistos hcell1(1,xcell_size_um,ycell_size_um);
  CellRelatedHistos hcell2(2,xcell_size_um,ycell_size_um);
  CellRelatedHistos hcell4(4,xcell_size_um,ycell_size_um);

  TH2I * dutpxxyHisto = new
    TH2I( "dutpxxy",
	  "DUT hits vs x;x [mm];y [mm];pixels",
	  nbx, -11, 11, nby, -6, 6 ); // bin = pix
  TProfile2D * effvsxy = new
    TProfile2D( "effvsxy",
		"DUT efficiency vs x;x track at DUT [mm];y track at DUT [mm];efficiency",
		nbx, -11, 11, nby, -6, 6, -1, 2 ); // bin = pix
  TProfile2D * effvsxydead = new
    TProfile2D( "effvsxydead",
		"DUT dead efficiency vs x;x track at DUT [mm];y track at DUT [mm];dead efficiency",
		nbx, -11, 11, nby, -6, 6, -1, 2 ); // bin = pix
  TProfile2D * effvsxylive = new
    TProfile2D( "effvsxylive",
		"DUT alive efficiency vs x;x track at DUT [mm];y track at DUT [mm];alive efficiency",
		nbx, -11, 11, nby, -6, 6, -1, 2 ); // bin = pix

  TProfile2D * effvsdxdydead = new
    TProfile2D( "effvsdxdydead",
		"efficiency around dead pixels;#Deltax track to dead pixel [mm];#Deltay track to dead pixel [mm];efficiency",
		60, -0.15, 0.15, 90, -4.5*0.025, 4.5*0.025, -1, 2 ); // for 100x25 straight

  TProfile effvsx( "effvsx",
		   "DUT efficiency vs x;x track at DUT [mm];efficiency",
		   nbx, -11, 11, -1, 2 ); // bin = col
  TProfile effvsx2( "effvsx2",
		   "DUT efficiency vs x;x track at DUT [mm];efficiency",
		   2*nbx, -11, 11, -1, 2 );
  TProfile effvsx4( "effvsx4",
		   "DUT efficiency vs x;x track at DUT [mm];efficiency",
		   4*nbx, -11, 11, -1, 2 );
  TProfile effvsy( "effvsy",
		   "DUT efficiency vs y;y track at DUT [mm];efficiency",
		   nby, -6, 6, -1, 2 ); // bin = row
  TProfile effvsydead( "effvsydead",
		       "DUT dead efficiency vs y;y track at DUT [mm];dead efficiency",
		       nby, -6, 6, -1, 2 ); // bin = row
  TProfile effvsylive( "effvsylive",
		       "DUT alive efficiency vs y;y track at DUT [mm];alive efficiency",
		       nby, -6, 6, -1, 2 ); // bin = row
  TProfile effvsr( "effvsr",
		   "DUT efficiency vs Rbeam;R track from beam [mm];efficiency",
		   100, 0, 10 );

  TProfile effvsev1( "effvsev1",
		     "DUT efficiency vs events;events;efficiency / 1k",
		     700, 0, 700*1000, -1, 2 );
  TProfile effvsev2( "effvsev2",
		     "DUT efficiency vs events;events;efficiency / 10k",
		     5000, 0, 50*1000*1000, -1, 2 );

  TProfile effvst1( "effvst1",
		    "DUT efficiency vs time;time [s];<efficiency> / s",
		    300, 0, 300, -1, 2 );
  TProfile effvst2( "effvst2",
		    "DUT efficiency vs time;time [s];<efficiency> / 5 s",
		    200, 0, 1000, -1, 2 );
  TProfile effvst3( "effvst3",
		    "DUT efficiency vs time;time [s];<efficiency> / 10 s",
		    600, 0, 6000, -1, 2 );
  TProfile effvst5( "effvst5",
		    "DUT efficiency vs time;time [s];<efficiency> / min",
		    1100, 0, 66000, -1, 2 );

  TH1I dutpdminHisto( "dutpdmin",
		      "pixel - Telescope min dxy;pixel - triplet min #Deltaxy [mm];tracks",
		      200, -1, 1 );

  TProfile effvsdxy( "effvsdxy",
		     "DUT efficiency vs dxy;xy search radius [mm];efficiency",
		     100, 0, 1, -1, 2 );

  TProfile effvsntri( "effvsntri",
		      "DUT efficiency vs triplets;triplets;efficiency",
		      20, 0.5, 20.5, -1, 2 );

  TProfile2D * effvsxmym = new
    TProfile2D( "effvsxmym",
		"DUT efficiency vs xmod ymod;x track mod 100 [#mum];y track mod 100 [#mum];efficiency",
		50, 0, 100, 50, 0, 100, -1, 2 );
  TProfile2D * effvsxmymdead = new
    TProfile2D( "effvsxmymdead",
		"DUT dead efficiency vs xmod ymod;x track mod 100 [#mum];y track mod 100 [#mum];dead efficiency",
		50, 0, 100, 50, 0, 100, -1, 2 );
  TProfile2D * effvsxmymlive = new
    TProfile2D( "effvsxmymlive",
		"DUT alive efficiency vs xmod ymod;x track mod 100 [#mum];y track mod 100 [#mum];alive efficiency",
		50, 0, 100, 50, 0, 100, -1, 2 );
  TProfile2D * effvsxmym2live = new
    TProfile2D( "effvsxmym2live",
		"DUT alive efficiency vs xmod ymod;x track mod 200 [#mum];y track mod 200 [#mum];alive efficiency",
		100, 0, 200, 100, 0, 200, -1, 2 );
  TProfile effvsxm( "effvsxm",
		    "DUT efficiency vs xmod;x track mod 100 [#mum];efficiency",
		    50, 0, 100, -1, 2 );
  TProfile effvsxm2( "effvsxm2",
		    "DUT efficiency vs xmod;x track mod 200 [#mum];efficiency",
		    100, 0, 200, -1, 2 );
  TProfile effvsxm4( "effvsxm4",
		    "DUT efficiency vs xmod;x track mod 400 [#mum];efficiency",
		    200, 0, 400, -1, 2 );
  TProfile effvsym( "effvsym",
		    "DUT efficiency vs ymod;y track mod 100 [#mum];efficiency",
		    50, 0, 100, -1, 2 );

  TProfile effvstx( "effvstx",
		    "DUT efficiency vs track slope x;x track slope [rad];efficiency",
		    100, -0.005, 0.005, -1, 2 );
  TProfile effvsty( "effvsty",
		    "DUT efficiency vs track slope y;y track slope [rad];efficiency",
		    100, -0.005, 0.005, -1, 2 );

  TProfile dutlkvst1( "dutlkvst1",
		      "track-DUT links vs time;time [s];tracks with DUT links / s",
		      300, 0, 300, -0.5, 1.5 );
  TProfile dutlkvst3( "dutlkvst3",
		      "track-DUT links vs time;time [s];tracks with DUT links / 10s",
		      300, 0, 3000, -0.5, 1.5 );
  TProfile dutlkvst5( "dutlkvst5",
		      "track-DUT links vs time;time [s];tracks with DUT links / min",
		      1100, 0, 66000, -0.5, 1.5 );

  TH1I ntrilkHisto( "ntrilk", "track - DUT links;track - DUT links;tracks",
		    11, -0.5, 10.5 );

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // event loop:

  cout << endl;

  FileReader * reader;
  if(      run <    100 )
    reader = new FileReader( runnum.c_str(), "data/run0000$2R$X" );
  else if( run <   1000 )
    reader = new FileReader( runnum.c_str(), "data/run000$3R$X" );
  else if( run <  10000 )
    reader = new FileReader( runnum.c_str(), "data/run00$4R$X" );
  else if( run < 100000 )
    reader = new FileReader( runnum.c_str(), "data/run0$5R$X" );
  else
    reader = new FileReader( runnum.c_str(), "data/run$6R$X" );

  DetectorEvent evt = reader->GetDetectorEvent();
  uint64_t evTLU0 = evt.GetTimestamp(); // 384 MHz = 2.6 ns
  if( evt.IsBORE() ) {
    eudaq::PluginManager::Initialize(evt);
    cout << "BORE TLU " << evTLU0 << endl;
    reader->NextEvent();
  }
  uint64_t prevTLU = evTLU0;

  int iev = 0;

  if( fev ) cout << "EU skip " << fev << endl;
  while( iev < fev ) {
    reader->NextEvent(); // fast forward
    ++iev;
  }

  int nevA = 0;
  int nevB = 0;
  int nmodlk = 0;

  int ntrck = 0;
  int ngood = 0;

  std::map<int,int> pxdutmap;

  do {

    evt = reader->GetDetectorEvent();

    uint64_t evTLU = evt.GetTimestamp(); // 384 MHz = 2.6 ns

    double evsec = (evTLU - evTLU0) / fTLU;
    t1Histo.Fill( evsec );
    t2Histo.Fill( evsec );
    t3Histo.Fill( evsec );
    t4Histo.Fill( evsec );
    t5Histo.Fill( evsec );
    t6Histo.Fill( evsec/3600 );

    double evdt = (evTLU - prevTLU) / fTLU;
    dtusHisto.Fill( evdt * 1E6 ); // [us]
    dtmsHisto.Fill( evdt * 1E3 ); // [ms]

    dt373Histo.Fill( (evTLU - prevTLU)%373 );
    dt374Histo.Fill( (evTLU - prevTLU)%374 );
    dt375Histo.Fill( (evTLU - prevTLU)%375 ); // best
    dt376Histo.Fill( (evTLU - prevTLU)%376 );
    dt377Histo.Fill( (evTLU - prevTLU)%377 );
    dt375vsdt.Fill( evdt*1E6, (evTLU - prevTLU)%375 ); // linear

    prevTLU = evTLU;

    bool ldbg = 0;

    if( iev <  0 )
      ldbg = 1;

    if( lev < 100 )
      ldbg = 1;

    if( iev < 10 || ldbg )
      cout << "scope53m  " << run << "." << iev << "  taken " << evsec << endl;
    else if( iev < 100 && iev%10 == 0 )
      cout << "scope53m  " << run << "." << iev << "  taken " << evsec << endl;
    else if( iev < 1000 && iev%100 == 0 )
      cout << "scope53m  " << run << "." << iev << "  taken " << evsec << endl;
    else if( iev%1000 == 0 ) {
      cout << "scope53m  " << run << "." << iev
	   << "  taken " << evsec
	   << "  modlk " << nmodlk
	   << "  eff " << ngood*1E2/max(1,ntrck)
	   << endl;
      ngood = 0;
      ntrck = 0;
      nmodlk = 0;
    }

    StandardEvent sevt = eudaq::PluginManager::ConvertToStandard(evt);

    if( ldbg ) cout << "planes " << sevt.NumPlanes() << endl;

    vector < cluster > cl[9];

    for( size_t iplane = 0; iplane < sevt.NumPlanes(); ++iplane ) {

      const eudaq::StandardPlane &plane = sevt.GetPlane(iplane);

      if(  ldbg )
	cout
	  << "  " << iplane
	  << ": plane " << plane.ID() // 1
	  << " " << plane.Type() // NI or BDAQ53
	  << " " << plane.Sensor() // MIMOSA26 or RD53A
	  << " frames " << plane.NumFrames() // 2 for NI or 32 for BDAQ53
	  << " pivot " << plane.PivotPixel()
	  << " total " << plane.TotalPixels()
	  << " hits " << plane.HitPixels()
	  ;

      int ipl = plane.ID(); // 0 = DUT, 1..6 = Mimosa
      // RD53A plane from converter, EUDAQ assigned id in [30 ,40)
      if(ipl >= 30 and ipl < 40)
      {
         ipl = 0;
      }

      if( ipl < 0 || ipl > 6 ) {
	cout << "event " << iev << " wrong plane number " << ipl << endl;
	continue;
      }

      hpivot[ipl].Fill( plane.PivotPixel() );
      hnpx[ipl].Fill( plane.HitPixels() ); // before masking
      hnframes[ipl].Fill( plane.NumFrames() ); // 32

      vector <pixel> pb; // for clustering

      // loop over frames, then pixels per frame

      for( unsigned frm = 0; frm < plane.NumFrames(); ++frm )

	for( size_t ipix = 0; ipix < plane.HitPixels( frm ); ++ipix ) {

	  if( ldbg )
	    cout << ": " << plane.GetX(ipix,frm)
		 << "." << plane.GetY(ipix,frm)
		 << "." << plane.GetPixel(ipix,frm) << " ";

	  int ix = plane.GetX(ipix,frm); // column
	  int iy = plane.GetY(ipix,frm); // row
	  int tot = plane.GetPixel(ipix,frm); // ToT 0..15

	  if( ipl == iDUT ) {
	    dutpxbcHisto.Fill( frm ); // before hot pixel masking
	    hrow[8].Fill( iy ); // before masking
	    hmap[8]->Fill( ix, iy ); // before masking
	    dutpxq8Histo.Fill( tot+1 );
	  }

	  // skip hot pixels:

	  int ipx = ix*ny[ipl] + iy;
	  if( hotset[ipl].count(ipx) ) continue;

	  if( ipl == iDUT ) { // after hot pixel masking

	    dutpxbcmHisto.Fill( frm ); // peaks at 7 and 9
	    if( tot+1 == 1 )
	      dutpxbc01Histo.Fill( frm );
	    else if( tot+1 == 15 )
	      dutpxbc15Histo.Fill( frm );
	    else
	      dutpxbc08Histo.Fill( frm );

	    dutpxbcmap->Fill( ix, iy, frm );
	    if( frm == 7 )
	      dutpxmapbc7->Fill( ix, iy );
	    else if ( frm == 9 )
	      dutpxmapbc9->Fill( ix, iy );

	    dutpxbcmap8->Fill( ix%8, iy%8, frm );
	    dutpxmap8->Fill( ix%8, iy%8 );
	    if( frm == 10 )
	      dutpxmap8bc10->Fill( ix%8, iy%8 );

	    dutpxbcvst.Fill( iev, frm ); // flat

	  }
	  
          // skip out-of-frame DUT pixels (the BC-ID) 
          // BUT Need a mechanism to avoid lost of efficiency!!
          //if( ipl == iDUT && ( frm < cuts.bcidmin || frm > cuts.bcidmax ) ) {
          //  continue;
          //}


	  pixel px;
	  px.col = ix; // ROC col
	  px.row = iy; // row
          if( ipl == iDUT )
          {
             int pcol = ix;
             int prow = iy;
	     if( !fifty ) 
	     {  // 100x25 from ROC to sensor:
	        pcol = ix/2; // 100 um

	        if( ix%2 == 1 ) 
                {
                    prow = 2*iy + 1; //
                }
	        else
                {
                    prow = 2*iy + 0;
                }
             }
             int thechan = pcol*nbr+prow;
             if( use_dut_calibration )
             {
               if(calibration_curves.find(thechan) == calibration_curves.end())
               { 
                  continue;
               }
               px.tot = calibration_curves[thechan](tot);
             }
          }
          else
          {
              px.tot = tot;
          }
	  px.frm = frm;
	  px.pivot = plane.GetPivot(ipix,frm);

	  if( ipl == iDUT ) 
          { 
	    // for hot pixel (c++11 already take cares of initialization)
	    if(create_duthotfile)
 	    {
                // Store ipx (see L2992 how is defined)
                ++pxdutmap[ipx];
            }

	    dutpxcol0Histo.Fill( ix + 0.5 );

	    px.tot += 1; // shift from zero: for COG clustering

	    if( ix < 128 )
	      cout << "pixel in Sync section: ignored" << endl;
	    else if( ix < 264 )
	      linpxqHisto.Fill( px.tot );
	    else
	      cout << "pixel in Diff section: ignored" << endl;

	    dutpxqvsx.Fill( ix, px.tot );
	    dutpxqvsxy->Fill( ix, iy, px.tot );

	    linpxbcHisto.Fill( frm );

	    if( px.tot < threshold ) continue; // offline threshold

	    dutpxcol9Histo.Fill( ix + 0.5 );

	    if( !fifty ) 
	    {  // 100x25 from ROC to sensor:
	       px.col = ix/2; // 100 um

	       if( ix%2 == 1 ) 
               {
                 px.row = 2*iy + 1; //
               }
	       else
               {
		 px.row = 2*iy + 0; // see ed53 for shallow angle
               }
	       if( chip0 == 182 || chip0 == 211 || chip0 == 512 ) 
               { // HLL
		if( ix%2 )
		  px.row = 2*iy + 1;
		else
		  px.row = 2*iy + 0;
	       }
	    }
	  } // DUT

	  pb.push_back(px);

	  hcol[ipl].Fill( ix ); // ROC
	  hrow[ipl].Fill( iy );
	  hmap[ipl]->Fill( ix, iy ); // after thr

	} // pix

      if( ldbg ) cout << endl;

      hnpxmsk[ipl].Fill( pb.size() ); // after masking
      if( ipl == iDUT )
	dutnpxvsev.Fill( iev, pb.size() );

      // clustering:

      if( ipl == iDUT )
	cl[ipl] = getClusq( pb, clflag );
      else
	cl[ipl] = getClusn( pb );

      if( ldbg ) cout << "    clusters " << cl[ipl].size() << endl;

      hncl[ipl].Fill( cl[ipl].size() );

      for( vector<cluster>::iterator c = cl[ipl].begin(); c != cl[ipl].end(); ++c ) {

	hsiz[ipl].Fill( c->size );
	hncol[ipl].Fill( c->ncol );
	hnrow[ipl].Fill( c->nrow );

	// cluster isolation:

	vector<cluster>::iterator d = c; // upper diagonal
	++d;
	for( ; d != cl[ipl].end(); ++d ) {
	  double dx = d->col - c->col;
	  double dy = d->row - c->row;
	  double dxy = sqrt( dx*dx + dy*dy );
	  if( dxy < c->mindxy ) c->mindxy = dxy;
	  if( dxy < d->mindxy ) d->mindxy = dxy;
	}

      } // cl

      for( vector<cluster>::iterator c = cl[ipl].begin(); c != cl[ipl].end(); ++c )
	hdxy[ipl].Fill( c->mindxy );

    } // eudaq planes

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // MOD:

    if( iev >= fev && modrun &&
	Astream.good() && ! Astream.eof() &&
	Bstream.good() && ! Bstream.eof()
	) {

      // one line = one trigger from one TBM channel

      string sl;
      getline( Astream, sl );

      istringstream Aev( sl ); // tokenize string

      //687
      //688 265 74 79 266 73 74 266 74 66

      int trg;
      Aev >> trg;

      ++nevA;

      int ierr = 0;

      bool ldb = 0;

      vector <pixel> pb; // for clustering

      while( ! Aev.eof()
	     // && Aev.good() && Aev.tellg() > 0 && Aev.tellg() < (int) sl.size()
	     ) {

	int xm;
	Aev >> xm; // 0..415

	if( Aev.eof() ) {
	  cout << "A truncated at line " << nevA << " col " << xm << endl;
	  cout << "Aev " << sl << endl;
	  ierr = 1;
	  break;
	}
	if( ldb ) cout << " : " << xm << flush;

	int ym;
	Aev >> ym; // 0..159

	if( Aev.eof() ) {
	  cout << "A truncated at line " << nevA << " row " << ym << endl;
	  cout << "Aev " << sl << endl;
	  ierr = 1;
	  break;
	}
	if( ldb ) cout << "." << ym << flush;

	int adc;
	Aev >> adc;
	if( ldb ) cout << "." << adc << flush;
	if( ldb ) cout << " (" << Aev.tellg() << ") " << flush;

	if( xm < 0 || xm > 415 || ym < 0 || ym > 159 || adc < 0 || adc > 255 ) {
	  cout << "data error at line " << nevA << " ev " << trg << endl;
	  ierr = 1;
	  break;
	}

	int roc = xm / 52; // 0..7
	int col = xm % 52; // 0..51
	//int row = ym;

	// leave space for big pixels:

	int ix = 1 + xm + 2*roc; // 1..52 per ROC
	int iy = ym;
	if( ym > 79 ) iy += 2; // 0..79, 82..161

	// fill pixel block for clustering:

	pixel px;
	px.col = ix; // col
	px.row = iy; // row
	px.tot = adc;
	pb.push_back(px);

	// double big pixels:
	// 0+1
	// 2..51
	// 52+53

	if( col == 0 ) {
	  px.col = ix-1; // double
	  px.row = iy;
	  pb[pb.size()-1].tot *= 0.5;
	  px.tot = 0.5*adc;
	  pb.push_back(px);
	}

	if( col == 51 ) {
	  px.col = ix+1; // double
	  px.row = iy;
	  pb[pb.size()-1].tot *= 0.5;
	  px.tot = 0.5*adc;
	  pb.push_back(px);
	}

	if( ym == 79 ) {
	  px.col = ix;
	  px.row = 80; // double
	  pb[pb.size()-1].tot *= 0.5;
	  px.tot = 0.5*adc;
	  pb.push_back(px);
	}

	if( ym == 80 ) {
	  px.col = ix;
	  px.row = 81; // double
	  pb[pb.size()-1].tot *= 0.5;
	  px.tot = 0.5*adc;
	  pb.push_back(px);
	}

      } // Aev

      if( ldb ) cout << endl;

      // B:

      getline( Bstream, sl );

      istringstream Bev( sl ); // tokenize string

      // 97875 358  81 47
      // 98144 373 126 47 372 126 46

      Bev >> trg;

      ++nevB; // lines

      if( ldb ) cout << trg << " B";

      //if( Bev.eof() ) { cout << "B empty" << endl; continue; }

      while( ! Bev.eof()
	     //&& Bev.good() && Bev.tellg() > 0 && Bev.tellg() < (int) sl.size()
	     ) {

	int xm;
	Bev >> xm; // 0..415

	if( Bev.eof() ) {
	  cout << "B truncated at line " << nevB << " col " << xm << endl;
	  cout << "Bev " << sl << endl;
	  ierr = 1;
	  break;
	}
	if( ldb ) cout << " : " << xm << flush;

	int ym;
	Bev >> ym; // 0..159

	if( Bev.eof() ) {
	  cout << "B truncated at line " << nevB << " row " << ym << endl;
	  cout << "Bev " << sl << endl;
	  ierr = 1;
	  break;
	}
	if( ldb ) cout << "." << ym << flush;

	int adc;
	Bev >> adc;
	if( ldb ) cout << "." << adc << flush;
	if( ldb ) cout << " (" << Bev.tellg() << ") " << flush;

	if( xm < 0 || xm > 415 || ym < 0 || ym > 159 || adc < 0 || adc > 255 ) {
	  cout << "data error at line " << nevB << " ev " << trg << endl;
	  ierr = 1;
	  break;
	}

	int roc = xm / 52; // 0..7
	int col = xm % 52; // 0..51
	//int row = ym;

	// leave space for big pixels:

	int ix = 1 + xm + 2*roc; // 1..52 per ROC
	int iy = ym;
	if( ym > 79 ) iy += 2; // 0..79, 82..161

	// fill pixel block for clustering:

	pixel px;
	px.col = ix; // col
	px.row = iy; // row
	px.tot = adc;
	pb.push_back(px);

	// double big pixels:
	// 0+1
	// 2..51
	// 52+53

	if( col == 0 ) {
	  px.col = ix-1; // double
	  px.row = iy;
	  pb[pb.size()-1].tot *= 0.5;
	  px.tot = 0.5*adc;
	  pb.push_back(px);
	}

	if( col == 51 ) {
	  px.col = ix+1; // double
	  px.row = iy;
	  pb[pb.size()-1].tot *= 0.5;
	  px.tot = 0.5*adc;
	  pb.push_back(px);
	}

	if( ym == 79 ) {
	  px.col = ix;
	  px.row = 80; // double
	  pb[pb.size()-1].tot *= 0.5;
	  px.tot = 0.5*adc;
	  pb.push_back(px);
	}

	if( ym == 80 ) {
	  px.col = ix;
	  px.row = 81; // double
	  pb[pb.size()-1].tot *= 0.5;
	  px.tot = 0.5*adc;
	  pb.push_back(px);
	}

      } // Bev

      if( ldb ) cout << endl;

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // clustering:

      hnpx[iMOD].Fill( pb.size() );

      if( !ierr )
	cl[iMOD] = getClusq( pb );

      hncl[iMOD].Fill( cl[iMOD].size() );

      for( vector<cluster>::iterator c = cl[iMOD].begin(); c != cl[iMOD].end(); ++c ) {

	hsiz[iMOD].Fill( c->size );
	hncol[iMOD].Fill( c->ncol );
	hnrow[iMOD].Fill( c->nrow );

	for( vector<pixel>::iterator px = c->vpix.begin(); px != c->vpix.end(); ++px ) {
	  hcol[iMOD].Fill( px->col );
	  hrow[iMOD].Fill( px->row );
	  hmap[iMOD]->Fill( px->col, px->row );
	}

      } // cl

    } // mod

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // DUT align vs event:

    if( run == 33333 ) { // dutdxvst5->Fit("pol9")

      double p0 =  8.99112e-05;
      double p1 =  5.22731e-07;
      double p2 = -1.76833e-10;
      double p3 =  3.12274e-14;
      double p4 = -3.10393e-18;
      double p5 =  1.86842e-22;
      double p6 = -6.75456e-27;
      double p7 =   1.4222e-31;
      double p8 = -1.60167e-36;
      double p9 =  7.44028e-42;
      DUTalignx = DUTalignx0 + p0 + ( p1 + ( p2 + ( p3 + ( p4 + ( p5 + ( p6 + ( p7 + ( p8 + p9 * evsec ) * evsec ) * evsec ) * evsec ) * evsec ) * evsec ) * evsec ) * evsec ) * evsec;

      p0 =   2.8632e-05;
      p1 = -5.09685e-07;
      p2 =  1.23226e-10;
      p3 = -3.65271e-15;
      p4 = -1.16106e-18;
      p5 =  1.31182e-22;
      p6 = -6.06594e-27;
      p7 =  1.43968e-31;
      p8 =  -1.7298e-36;
      p9 =  8.35097e-42;
      DUTaligny = DUTaligny0 + p0 + ( p1 + ( p2 + ( p3 + ( p4 + ( p5 + ( p6 + ( p7 + ( p8 + p9 * evsec ) * evsec ) * evsec ) * evsec ) * evsec ) * evsec ) * evsec ) * evsec ) * evsec;

    }

    if( run == 33373 ) { // dutdxvst5->Fit("pol9")

      double p0 =   0.00460585;
      double p1 =  6.49668e-07;
      double p2 = -9.02168e-11;
      double p3 = -6.19962e-15;
      double p4 =  1.58327e-18;
      double p5 =   -1.069e-22;
      double p6 =  3.58945e-27;
      double p7 = -6.56699e-32;
      double p8 =  6.26488e-37;
      double p9 =   -2.444e-42;
      DUTalignx = DUTalignx0 + p0 + ( p1 + ( p2 + ( p3 + ( p4 + ( p5 + ( p6 + ( p7 + ( p8 + p9 * evsec ) * evsec ) * evsec ) * evsec ) * evsec ) * evsec ) * evsec ) * evsec ) * evsec;

      p0 =  2.98704e-05;
      p1 = -8.07885e-07;
      p2 =  1.76479e-10;
      p3 = -1.63139e-14;
      p4 =   7.7256e-19;
      p5 = -2.01131e-23;
      p6 =  2.87059e-28;
      p7 =  -1.9848e-33;
      p8 =  3.16939e-39;
      p9 =   1.9522e-44;
      DUTaligny = DUTaligny0 + p0 + ( p1 + ( p2 + ( p3 + ( p4 + ( p5 + ( p6 + ( p7 + ( p8 + p9 * evsec ) * evsec ) * evsec ) * evsec ) * evsec ) * evsec ) * evsec ) * evsec ) * evsec;

    }

    if( run == 34135 ) { // dutdxvst5->Fit("pol9")

      double p0 =   0.00112574;
      double p1 = -1.06476e-06;
      double p2 =  5.26923e-10;
      double p3 = -4.78558e-14;
      double p4 =  1.43135e-18;
      double p5 =  1.94292e-23;
      double p6 = -2.33677e-27;
      double p7 =  6.10438e-32;
      double p8 = -7.07503e-37;
      double p9 =  3.16034e-42;
      DUTalignx = DUTalignx0 + p0 + ( p1 + ( p2 + ( p3 + ( p4 + ( p5 + ( p6 + ( p7 + ( p8 + p9 * evsec ) * evsec ) * evsec ) * evsec ) * evsec ) * evsec ) * evsec ) * evsec ) * evsec;

      p0 =  -0.00109749;
      p1 =  7.39726e-07;
      p2 = -3.71976e-10;
      p3 =  4.52165e-14;
      p4 = -2.57825e-18;
      p5 =  7.82913e-23;
      p6 = -1.24621e-27;
      p7 =  8.05672e-33;
      p8 =  2.02113e-38;
      p9 = -3.72545e-43;
      DUTaligny = DUTaligny0 + p0 + ( p1 + ( p2 + ( p3 + ( p4 + ( p5 + ( p6 + ( p7 + ( p8 + p9 * evsec ) * evsec ) * evsec ) * evsec ) * evsec ) * evsec ) * evsec ) * evsec ) * evsec;

    }

    if( run == 35119 ) { // dutdxvst5->Fit("pol9")

      double p0 = -0.000233727;
      double p1 =  9.71291e-07;
      double p2 = -1.02976e-10;
      double p3 =  3.67446e-14;
      double p4 = -6.69844e-18;
      double p5 =  6.18927e-22;
      double p6 = -3.15456e-26;
      double p7 =  8.99555e-31;
      double p8 = -1.34561e-35;
      double p9 =  8.22931e-41;
      DUTalignx = DUTalignx0 + p0 + ( p1 + ( p2 + ( p3 + ( p4 + ( p5 + ( p6 + ( p7 + ( p8 + p9 * evsec ) * evsec ) * evsec ) * evsec ) * evsec ) * evsec ) * evsec ) * evsec ) * evsec;

      p0 =   0.00019829;
      p1 = -5.29521e-07;
      p2 =  6.42272e-11;
      p3 = -1.59508e-14;
      p4 =  3.09556e-18;
      p5 = -3.29979e-22;
      p6 =  1.90697e-26;
      p7 = -6.00906e-31;
      p8 =  9.73253e-36;
      p9 = -6.34917e-41;
      DUTaligny = DUTaligny0 + p0 + ( p1 + ( p2 + ( p3 + ( p4 + ( p5 + ( p6 + ( p7 + ( p8 + p9 * evsec ) * evsec ) * evsec ) * evsec ) * evsec ) * evsec ) * evsec ) * evsec ) * evsec;

    }

    if( run == 35590 ) { // dutdxvst5->Fit("pol9")

      double p0 =   -0.0065022;
      double p1 =  5.39455e-07;
      double p2 = -1.35292e-11;
      double p3 =  1.49578e-14;
      double p4 = -1.75745e-18;
      double p5 =  9.13022e-23;
      double p6 = -2.57646e-27;
      double p7 =  4.09995e-32;
      double p8 = -3.46008e-37;
      double p9 =  1.20437e-42;
      DUTalignx = DUTalignx0 + p0 + ( p1 + ( p2 + ( p3 + ( p4 + ( p5 + ( p6 + ( p7 + ( p8 + p9 * evsec ) * evsec ) * evsec ) * evsec ) * evsec ) * evsec ) * evsec ) * evsec ) * evsec;

      p0 =  -0.00111423;
      p1 = -2.47965e-07;
      p2 = -1.27336e-11;
      p3 = -2.08114e-15;
      p4 =  3.72269e-19;
      p5 =  -2.1894e-23;
      p6 =  6.55495e-28;
      p7 =  -1.0745e-32;
      p8 =  9.17607e-38;
      p9 = -3.19272e-43;
      DUTaligny = DUTaligny0 + p0 + ( p1 + ( p2 + ( p3 + ( p4 + ( p5 + ( p6 + ( p7 + ( p8 + p9 * evsec ) * evsec ) * evsec ) * evsec ) * evsec ) * evsec ) * evsec ) * evsec ) * evsec;

    }

    if( run == 35670 ) { // dutdxvst5->Fit("pol9")

      double p0 =  0.000361424;
      double p1 = -2.02731e-06;
      double p2 =  7.27207e-10;
      double p3 = -3.04223e-13;
      double p4 =   6.4932e-17;
      double p5 = -7.13696e-21;
      double p6 =  4.33966e-25;
      double p7 = -1.47014e-29;
      double p8 =  2.57816e-34;
      double p9 = -1.80007e-39;
      DUTalignx = DUTalignx0 + p0 + ( p1 + ( p2 + ( p3 + ( p4 + ( p5 + ( p6 + ( p7 + ( p8 + p9 * evsec ) * evsec ) * evsec ) * evsec ) * evsec ) * evsec ) * evsec ) * evsec ) * evsec;

      p0 = -0.000258553;
      p1 =  1.22008e-06;
      p2 = -3.45351e-10;
      p3 =  1.41691e-13;
      p4 = -3.55133e-17;
      p5 =  4.60807e-21;
      p6 = -3.29672e-25;
      p7 =  1.32183e-29;
      p8 = -2.79026e-34;
      p9 =  2.41744e-39;
      DUTaligny = DUTaligny0 + p0 + ( p1 + ( p2 + ( p3 + ( p4 + ( p5 + ( p6 + ( p7 + ( p8 + p9 * evsec ) * evsec ) * evsec ) * evsec ) * evsec ) * evsec ) * evsec ) * evsec ) * evsec;

    }

    if( run == 38095 ) { // dutdxvst2->Fit("pol9")

      double p0 =  -0.00145845;
      double p1 = -6.12614e-05;
      double p2 =  2.07746e-06;
      double p3 = -2.33604e-08;
      double p4 =  1.32488e-10;
      double p5 = -4.28366e-13;
      double p6 =  8.23075e-16;
      double p7 = -9.28281e-19;
      double p8 =  5.66396e-22;
      double p9 = -1.44019e-25;
      DUTalignx = DUTalignx0 + p0 + ( p1 + ( p2 + ( p3 + ( p4 + ( p5 + ( p6 + ( p7 + ( p8 + p9 * evsec ) * evsec ) * evsec ) * evsec ) * evsec ) * evsec ) * evsec ) * evsec ) * evsec;

    }

    if( run == 38096 ) { // dutdxvst2->Fit("pol9")

      double p0 =    -0.002503;
      double p1 = -4.63296e-05;
      double p2 =  1.38021e-06;
      double p3 = -1.40307e-08;
      double p4 =  7.59069e-11;
      double p5 = -2.40585e-13;
      double p6 =  4.60265e-16;
      double p7 = -5.22692e-19;
      double p8 =  3.24234e-22;
      double p9 = -8.45737e-26;
      DUTalignx = DUTalignx0 + p0 + ( p1 + ( p2 + ( p3 + ( p4 + ( p5 + ( p6 + ( p7 + ( p8 + p9 * evsec ) * evsec ) * evsec ) * evsec ) * evsec ) * evsec ) * evsec ) * evsec ) * evsec;

    }

    if( run == 38097 ) { // dutdxvst2->Fit("pol9")

      double p0 =  -0.00257178;
      double p1 =  7.23449e-05;
      double p2 = -8.85536e-07;
      double p3 =  7.00554e-09;
      double p4 = -3.58695e-11;
      double p5 =  1.23641e-13;
      double p6 = -2.79799e-16;
      double p7 =   3.8982e-19;
      double p8 = -2.99069e-22;
      double p9 =  9.59883e-26;
      DUTalignx = DUTalignx0 + p0 + ( p1 + ( p2 + ( p3 + ( p4 + ( p5 + ( p6 + ( p7 + ( p8 + p9 * evsec ) * evsec ) * evsec ) * evsec ) * evsec ) * evsec ) * evsec ) * evsec ) * evsec;

    }

    if( run == 38105 ) { // dutdxvst2->Fit("pol9")

      double p0 =   0.00557853;
      double p1 = -7.22542e-05;
      double p2 =  1.15178e-06;
      double p3 = -1.32007e-08;
      double p4 =  8.64756e-11;
      double p5 = -3.33374e-13;
      double p6 =  7.74128e-16;
      double p7 =  -1.0657e-18;
      double p8 =  8.01204e-22;
      double p9 = -2.53382e-25;
      DUTalignx = DUTalignx0 + p0 + ( p1 + ( p2 + ( p3 + ( p4 + ( p5 + ( p6 + ( p7 + ( p8 + p9 * evsec ) * evsec ) * evsec ) * evsec ) * evsec ) * evsec ) * evsec ) * evsec ) * evsec;

    }

    if( run == 38106 ) { // dutdxvst2->Fit("pol9")

      double p0 =  0.000339581;
      double p1 =  8.78933e-05;
      double p2 = -2.27439e-06;
      double p3 =  2.23216e-08;
      double p4 = -1.15775e-10;
      double p5 =  3.45126e-13;
      double p6 = -6.04512e-16;
      double p7 =  6.05458e-19;
      double p8 =  -3.1325e-22;
      double p9 =  6.19965e-26;
      DUTalignx = DUTalignx0 + p0 + ( p1 + ( p2 + ( p3 + ( p4 + ( p5 + ( p6 + ( p7 + ( p8 + p9 * evsec ) * evsec ) * evsec ) * evsec ) * evsec ) * evsec ) * evsec ) * evsec ) * evsec;

    }

    if( run == 38107 ) { // dutdxvst2->Fit("pol9")

      double p0 =   0.00202886;
      double p1 =  4.13588e-05;
      double p2 = -1.60627e-06;
      double p3 =  1.98847e-08;
      double p4 = -1.23116e-10;
      double p5 =  4.28409e-13;
      double p6 = -8.75501e-16;
      double p7 =   1.0428e-18;
      double p8 = -6.70107e-22;
      double p9 =  1.79548e-25;
      DUTalignx = DUTalignx0 + p0 + ( p1 + ( p2 + ( p3 + ( p4 + ( p5 + ( p6 + ( p7 + ( p8 + p9 * evsec ) * evsec ) * evsec ) * evsec ) * evsec ) * evsec ) * evsec ) * evsec ) * evsec;

    }

    if( run == 38108 ) { // dutdxvst2->Fit("pol9")

      double p0 =   0.00432767;
      double p1 = -2.93722e-05;
      double p2 = -7.34686e-07;
      double p3 =  1.51137e-08;
      double p4 = -1.16217e-10;
      double p5 =   4.7268e-13;
      double p6 =  -1.1088e-15;
      double p7 =  1.50735e-18;
      double p8 = -1.10427e-21;
      double p9 =  3.37592e-25;
      DUTalignx = DUTalignx0 + p0 + ( p1 + ( p2 + ( p3 + ( p4 + ( p5 + ( p6 + ( p7 + ( p8 + p9 * evsec ) * evsec ) * evsec ) * evsec ) * evsec ) * evsec ) * evsec ) * evsec ) * evsec;

    }

    if( run == 38109 ) { // dutdxvst2->Fit("pol9")

      double p0 =   0.00400614;
      double p1 =  0.000185384;
      double p2 = -5.13357e-06;
      double p3 =  6.24512e-08;
      double p4 = -4.13974e-10;
      double p5 =  1.58577e-12;
      double p6 = -3.60821e-15;
      double p7 =  4.81435e-18;
      double p8 = -3.47961e-21;
      double p9 =  1.05183e-24;
      DUTalignx = DUTalignx0 + p0 + ( p1 + ( p2 + ( p3 + ( p4 + ( p5 + ( p6 + ( p7 + ( p8 + p9 * evsec ) * evsec ) * evsec ) * evsec ) * evsec ) * evsec ) * evsec ) * evsec ) * evsec;

    }

    if( run == 38111 ) { // dutdxvst2->Fit("pol9")

      double p0 =   0.00202815;
      double p1 =  2.99177e-05;
      double p2 = -9.52802e-07;
      double p3 =  8.64039e-09;
      double p4 = -4.19191e-11;
      double p5 =  1.07239e-13;
      double p6 = -1.23459e-16;
      double p7 =   8.0915e-22;
      double p8 =  1.20339e-22;
      double p9 = -7.35116e-26;
      DUTalignx = DUTalignx0 + p0 + ( p1 + ( p2 + ( p3 + ( p4 + ( p5 + ( p6 + ( p7 + ( p8 + p9 * evsec ) * evsec ) * evsec ) * evsec ) * evsec ) * evsec ) * evsec ) * evsec ) * evsec;

    }

    if( run == 38112 ) { // dutdxvst2->Fit("pol9")

      double p0 =   0.00260416;
      double p1 = -0.000136501;
      double p2 =   2.2008e-06;
      double p3 = -1.98835e-08;
      double p4 =  1.02248e-10;
      double p5 = -3.13804e-13;
      double p6 =  5.86944e-16;
      double p7 = -6.58118e-19;
      double p8 =  4.08061e-22;
      double p9 = -1.07997e-25;
      DUTalignx = DUTalignx0 + p0 + ( p1 + ( p2 + ( p3 + ( p4 + ( p5 + ( p6 + ( p7 + ( p8 + p9 * evsec ) * evsec ) * evsec ) * evsec ) * evsec ) * evsec ) * evsec ) * evsec ) * evsec;

    }

    if( run == 38117 ) { // dutdxvst2->Fit("pol9")

      double p0 =   0.00210596;
      double p1 = -0.000131263;
      double p2 =  3.78765e-06;
      double p3 =  -5.2112e-08;
      double p4 =   3.5484e-10;
      double p5 =  -1.3473e-12;
      double p6 =  3.00298e-15;
      double p7 = -3.90808e-18;
      double p8 =   2.7495e-21;
      double p9 = -8.08009e-25;
      DUTalignx = DUTalignx0 + p0 + ( p1 + ( p2 + ( p3 + ( p4 + ( p5 + ( p6 + ( p7 + ( p8 + p9 * evsec ) * evsec ) * evsec ) * evsec ) * evsec ) * evsec ) * evsec ) * evsec ) * evsec;

    }

    if( run == 38121 ) { // dutdxvst2->Fit("pol9")

      double p0 =   -0.0018171;
      double p1 =  6.92937e-06;
      double p2 =  7.64688e-07;
      double p3 = -1.17948e-08;
      double p4 =    8.575e-11;
      double p5 = -3.52165e-13;
      double p6 =  8.58283e-16;
      double p7 = -1.22876e-18;
      double p8 =  9.53647e-22;
      double p9 = -3.09386e-25;
      DUTalignx = DUTalignx0 + p0 + ( p1 + ( p2 + ( p3 + ( p4 + ( p5 + ( p6 + ( p7 + ( p8 + p9 * evsec ) * evsec ) * evsec ) * evsec ) * evsec ) * evsec ) * evsec ) * evsec ) * evsec;

      p0 =  -0.00148988;
      p1 =  0.000183153;
      p2 = -4.01971e-06;
      p3 =  3.72566e-08;
      p4 = -1.84841e-10;
      p5 =  5.28721e-13;
      p6 = -8.80165e-16;
      p7 =  8.09544e-19;
      p8 = -3.48618e-22;
      p9 =   3.8567e-26;
      DUTaligny = DUTaligny0 + p0 + ( p1 + ( p2 + ( p3 + ( p4 + ( p5 + ( p6 + ( p7 + ( p8 + p9 * evsec ) * evsec ) * evsec ) * evsec ) * evsec ) * evsec ) * evsec ) * evsec ) * evsec;

    }

    if( run == 38446 ) { // dutdxvst5->Fit("pol9")

      double p0 =  0.000672108;
      double p1 = -4.43748e-06;
      double p2 =  1.58843e-09;
      double p3 = -3.63893e-13;
      double p4 =  5.97706e-17;
      double p5 = -7.05132e-21;
      double p6 =  5.61676e-25;
      double p7 =  -2.7949e-29;
      double p8 =  7.75382e-34;
      double p9 = -9.10656e-39;
      DUTalignx = DUTalignx0 + p0 + ( p1 + ( p2 + ( p3 + ( p4 + ( p5 + ( p6 + ( p7 + ( p8 + p9 * evsec ) * evsec ) * evsec ) * evsec ) * evsec ) * evsec ) * evsec ) * evsec ) * evsec;

      p0 = -0.000411821;
      p1 =  1.66084e-06;
      p2 = -7.27193e-10;
      p3 =   2.1882e-13;
      p4 = -4.84624e-17;
      p5 =  7.15116e-21;
      p6 = -6.52794e-25;
      p7 =  3.50708e-29;
      p8 = -1.01324e-33;
      p9 =   1.2123e-38;
      DUTaligny = DUTaligny0 + p0 + ( p1 + ( p2 + ( p3 + ( p4 + ( p5 + ( p6 + ( p7 + ( p8 + p9 * evsec ) * evsec ) * evsec ) * evsec ) * evsec ) * evsec ) * evsec ) * evsec ) * evsec;

    }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // RD53A cross talk correction:

    if( ! fifty ) {

      double cx = 0.00; // cross talk
      //double cx = 0.04; // cross talk
      //double cx = 0.08; // cross talk
      //double cx = 0.12; // cross talk

      for( vector<cluster>::iterator c = cl[iDUT].begin(); c != cl[iDUT].end(); ++c ) {

	// a cluster is a fuzzy object: use a map of maps

	map < int, map < int, double > > crq; // col row tot, sorted

	for( vector<pixel>::iterator px = c->vpix.begin(); px != c->vpix.end(); ++px )
	  crq[px->col][px->row] = px->tot;

	for( auto cc = crq.begin(); cc != crq.end(); ++cc ) { // columns

	  // cc->second is the sorted list of rows in this column

	  if( cc->second.size() == 1 ) continue;

	  for( auto rr = cc->second.begin(); rr != cc->second.end(); ++rr ) { // rows

	    if( rr->first%2 ) continue; // we want even rows

	    auto r1 = rr;
	    ++r1; // next row
	    if( r1 == cc->second.end() ) break;

	    if( r1->first == rr->first + 1 ) { // cluster may have hole in rows

	      double q0 = rr->second;
	      double q1 = r1->second;

	      // unfold = invert cross talk matrix:

	      double det = 1-cx*cx; // Determinante

	      rr->second = ( q0 - cx*q1 ) / det; // overwrite!
	      r1->second = ( q1 - cx*q0 ) / det; // U1+U2 = q12

	    } // r1

	  } // rr

	} // cc

	double sumq = 0;
	double sumqr = 0;

	for( vector<pixel>::iterator px = c->vpix.begin(); px != c->vpix.end(); ++px ) {
	  px->tot = crq[px->col][px->row]; // overwrite
	  sumq += px->tot;
	  sumqr += px->tot * px->row;
	}

	c->row = sumqr / sumq; // overwrite !

      } // c

    } // 100x25

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // DUT:

    for( vector<cluster>::iterator c = cl[iDUT].begin(); c != cl[iDUT].end(); ++c ) {

      double ccol = c->col;
      double crow = c->row;

      double dutx = ( ccol + 0.5 - nx[iDUT]/2 ) * ptchx[iDUT]; // mm
      double duty = ( crow + 0.5 - ny[iDUT]/2 ) * ptchy[iDUT]; // mm

      // Mod:

      for( vector<cluster>::iterator c = cl[iMOD].begin(); c != cl[iMOD].end(); ++c ) {

	double ccol = c->col;
	double crow = c->row;
	double modx = ( ccol + 0.5 - nx[iMOD]/2 ) * ptchx[iMOD]; // -33..33 mm
	double mody = ( crow + 0.5 - ny[iMOD]/2 ) * ptchy[iMOD]; // -8..8 mm

	dutmodxxHisto->Fill( dutx, modx );
	dutmodyyHisto->Fill( duty, mody );

      } // Mod

      // cluster-cluster:

      int minx = 999;
      int maxx = 0;
      int miny = 999;
      int maxy = 0;
      for( vector<pixel>::iterator p = c->vpix.begin(); p != c->vpix.end(); ++p ) {
	if( p->col > maxx ) maxx = p->col;
	if( p->col < minx ) minx = p->col;
	if( p->row > maxy ) maxy = p->row;
	if( p->row < miny ) miny = p->row;
      }

      int mindx = 999;
      int mindy = 999;

      vector<cluster>::iterator c2 = c;
      ++c2;

      for( ; c2 != cl[iDUT].end(); ++c2 ) {

	for( vector<pixel>::iterator p = c2->vpix.begin(); p != c2->vpix.end(); ++p ) {
	  int dx = p->col - maxx;
	  if( p->col < minx ) dx = minx - p->col;
	  if( dx < mindx && dx > 0 ) mindx = dx;
	  int dy = p->row - maxy;
	  if( p->row < miny ) dy = miny - p->row;
	  if( dy < mindy && dy > 0 ) mindy = dy;
	}

      } // c2

      if( mindx < mindy )
	dutdpxHisto.Fill( mindx );
      else
	dutdpxHisto.Fill( mindy );

    } // DUT

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // make triplets 1+3-2:

    vector <triplet> triplets;

    //double triCut = 0.1; // [mm]
    double triCut = 0.05; // [mm] like tele

    for( vector<cluster>::iterator cA = cl[1].begin(); cA != cl[1].end(); ++cA ) {

      double xA = cA->col*ptchx[1] - alignx[1];
      double yA = cA->row*ptchy[1] - aligny[1];
      double xmid = xA - midx[1];
      double ymid = yA - midy[1];
      xA = xmid - ymid*rotx[1];
      yA = ymid + xmid*roty[1];

      double zA = zz[1] + alignz[1];
      double zC = zz[3] + alignz[3];
      double zB = zz[2] + alignz[2];

      for( vector<cluster>::iterator cC = cl[3].begin(); cC != cl[3].end(); ++cC ) {

	double xC = cC->col*ptchx[3] - alignx[3];
	double yC = cC->row*ptchy[3] - aligny[3];
	double xmid = xC - midx[3];
	double ymid = yC - midy[3];
	xC = xmid - ymid*rotx[3];
	yC = ymid + xmid*roty[3];

	double dx2 = xC - xA;
	double dy2 = yC - yA;
	double dzCA = zC - zA;
	hdx13.Fill( dx2 );
	hdy13.Fill( dy2 );

	if( fabs( dx2 ) > 0.005*f * dzCA ) continue; // angle cut *f?
	if( fabs( dy2 ) > 0.005*f * dzCA ) continue; // angle cut

	double avx = 0.5 * ( xA + xC ); // mid
	double avy = 0.5 * ( yA + yC );
	double avz = 0.5 * ( zA + zC ); // mid z

	double slpx = ( xC - xA ) / dzCA; // slope x
	double slpy = ( yC - yA ) / dzCA; // slope y

	// middle plane B = 2:

	for( vector<cluster>::iterator cB = cl[2].begin(); cB != cl[2].end(); ++cB ) {

	  double xB = cB->col*ptchx[2] - alignx[2];
	  double yB = cB->row*ptchy[2] - aligny[2];
	  double xmid = xB - midx[2];
	  double ymid = yB - midy[2];
	  xB = xmid - ymid*rotx[2];
	  yB = ymid + xmid*roty[2];

	  // interpolate track to B:

	  double dz = zB - avz;
	  double xm = avx + slpx * dz; // triplet at mid
	  double ym = avy + slpy * dz;

	  double dxm = xB - xm;
	  double dym = yB - ym;
	  htridx.Fill( dxm );
	  htridy.Fill( dym );

	  if( fabs( dym ) < 0.05 ) {

	    htridxc.Fill( dxm );
	    tridxvsx.Fill( xB, dxm );
	    trimadxvsx.Fill( xB, fabs(dxm)*1e3 );
	    tridxvsy.Fill( yB, dxm );
	    tridxvstx.Fill( slpx, dxm );
	    tridxvst3.Fill( evsec, dxm );
	    tridxvst5.Fill( evsec, dxm );

	  } // dy

	  if( fabs( dxm ) < 0.05 ) {
	    htridyc.Fill( dym );
	    tridyvsx.Fill( xB, dym );
	    tridyvsty.Fill( slpy, dym );
	    tridyvst3.Fill( evsec, dym );
	    tridyvst5.Fill( evsec, dym );
	  }

	  // telescope triplet cuts:

	  if( fabs(dxm) > triCut ) continue;
	  if( fabs(dym) > triCut ) continue;

	  triplet tri;
	  tri.xm = avx;
	  tri.ym = avy;
	  tri.zm = avz;
	  tri.sx = slpx;
	  tri.sy = slpy;
	  tri.lk = 0;
	  tri.ttdmin = 99.9; // isolation [mm]
	  tri.iA = distance( cl[1].begin(), cA );
	  tri.iB = distance( cl[2].begin(), cB );
	  tri.iC = distance( cl[3].begin(), cC );

	  triplets.push_back(tri);

	  trixHisto.Fill( avx );
	  triyHisto.Fill( avy );
	  trixyHisto->Fill( avx, avy );
	  tritxHisto.Fill( slpx*1e3 );
	  tritxwvsx.Fill( avx, fabs(slpx)*1e3 );
	  trityHisto.Fill( slpy*1e3 );

	} // cl B

      } // cl C

    } // cl A

    ntriHisto.Fill( triplets.size() );

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // make driplets 4+6-5:

    vector <triplet> driplets;

    //double driCut = 0.1; // [mm]
    double driCut = 0.05; // [mm] like tele

    for( vector<cluster>::iterator cA = cl[4].begin(); cA != cl[4].end(); ++cA ) {

      double xA = cA->col*ptchx[4] - alignx[4];
      double yA = cA->row*ptchy[4] - aligny[4];
      double zA = zz[4] + alignz[4];
      double xmid = xA - midx[4];
      double ymid = yA - midy[4];
      xA = xmid - ymid*rotx[4];
      yA = ymid + xmid*roty[4];

      double zC = zz[6] + alignz[6];
      double zB = zz[5] + alignz[5];

      for( vector<cluster>::iterator cC = cl[6].begin(); cC != cl[6].end(); ++cC ) {

	double xC = cC->col*ptchx[6] - alignx[6];
	double yC = cC->row*ptchy[6] - aligny[6];
	double xmid = xC - midx[6];
	double ymid = yC - midy[6];
	xC = xmid - ymid*rotx[6];
	yC = ymid + xmid*roty[6];

	double dx2 = xC - xA;
	double dy2 = yC - yA;
	double dzCA = zC - zA; // from 4 to 6 in z
	hdx46.Fill( dx2 );
	hdy46.Fill( dy2 );

	if( fabs( dx2 ) > 0.005 * dzCA ) continue; // angle cut *f?
	if( fabs( dy2 ) > 0.005 * dzCA ) continue; // angle cut

	double avx = 0.5 * ( xA + xC ); // mid
	double avy = 0.5 * ( yA + yC );
	double avz = 0.5 * ( zA + zC ); // mid z

	double slpx = ( xC - xA ) / dzCA; // slope x
	double slpy = ( yC - yA ) / dzCA; // slope y

	// middle plane B = 5:

	for( vector<cluster>::iterator cB = cl[5].begin(); cB != cl[5].end(); ++cB ) {

	  double xB = cB->col*ptchx[5] - alignx[5];
	  double yB = cB->row*ptchy[5] - aligny[5];
	  double xmid = xB - midx[5];
	  double ymid = yB - midy[5];
	  xB = xmid - ymid*rotx[5];
	  yB = ymid + xmid*roty[5];

	  // interpolate track to B:

	  double dz = zB - avz;
	  double xm = avx + slpx * dz; // driplet at m
	  double ym = avy + slpy * dz;

	  double dxm = xB - xm;
	  double dym = yB - ym;
	  hdridx.Fill( dxm );
	  hdridy.Fill( dym );

	  if( fabs( dym ) < 0.05 ) {
	    hdridxc.Fill( dxm );
	    dridxvsy.Fill( yB, dxm );
	    dridxvstx.Fill( slpx, dxm );
	    dridxvst3.Fill( evsec, dxm );
	    dridxvst5.Fill( evsec, dxm );
	  }

	  if( fabs( dxm ) < 0.05 ) {
	    hdridyc.Fill( dym );
	    dridyvsx.Fill( xB, dym );
	    dridyvsty.Fill( slpy, dym );
	    dridyvst3.Fill( evsec, dym );
	    dridyvst5.Fill( evsec, dym );
	  }

	  // telescope driplet cuts:

	  if( fabs(dxm) > driCut ) continue;
	  if( fabs(dym) > driCut ) continue;

	  triplet dri;

	  dri.xm = avx;
	  dri.ym = avy;
	  dri.zm = avz;
	  dri.sx = slpx;
	  dri.sy = slpy;
	  dri.lk = 0;
	  dri.ttdmin = 99.9; // isolation [mm]
	  dri.iA = distance( cl[4].begin(), cA );
	  dri.iB = distance( cl[5].begin(), cB );
	  dri.iC = distance( cl[6].begin(), cC );

	  driplets.push_back(dri);

	  drixHisto.Fill( avx );
	  driyHisto.Fill( avy );
	  drixyHisto->Fill( avx, avy );
	  dritxHisto.Fill( slpx );
	  drityHisto.Fill( slpy );

	} // cl B

      } // cl C

    } // cl A

    ndriHisto.Fill( driplets.size() );

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // which Mod hits have a track:

    if( ldbmod )
      cout << iev << endl;

    for( vector<cluster>::iterator c = cl[iMOD].begin(); c != cl[iMOD].end(); ++c ) {

      double ccol = c->col;
      double crow = c->row;

      modcolHisto.Fill( ccol );
      modrowHisto.Fill( crow );

      double modx = ( ccol + 0.5 - nx[iMOD]/2 ) * ptchx[iMOD]; // -33..33 mm
      double mody = ( crow + 0.5 - ny[iMOD]/2 ) * ptchy[iMOD]; // -8..8 mm

      modxHisto.Fill( modx );
      modyHisto.Fill( mody );

      if( ldbmod )
	cout << "         " << modx << "  " << mody << endl;

      double mindy = 99;
      double mindx = 99;

      // triplets:

      for( unsigned int iA = 0; iA < triplets.size(); ++iA ) { // iA = upstream

	double xmA = triplets[iA].xm;
	double ymA = triplets[iA].ym;
	double zmA = triplets[iA].zm;
	double sxA = triplets[iA].sx;
	double syA = triplets[iA].sy;

	// intersect inclined track with tilted MOD plane:

	double zB = MODz - zmA; // z MOD from mid of triplet
	double zd = (Nzm*zB - Nym*ymA - Nxm*xmA) / (Nxm*sxA + Nym*syA + Nzm); // from zmB
	double yd = ymA + syA * zd;
	double xd = xmA + sxA * zd;

	double dzd = zd + zmA - MODz; // from MOD z0 [-8,8] mm

	// transform into MOD system: (passive).
	// large rotations don't commute: careful with order

	double x1m = com*xd - som*dzd; // turn o
	double y1m = yd;
	double z1m = som*xd + com*dzd;

	double x2m = x1m;
	double y2m = cam*y1m + sam*z1m; // tilt a

	double x3m = cfm*x2m + sfm*y2m; // rot
	double y3m =-sfm*x2m + cfm*y2m;

	double x4m =-x3m + MODalignx; // shift to mid
	double y4m = y3m + MODaligny; // invert y, shift to mid

	double moddx = modx - x4m;
	double moddy = mody - y4m;
	if( fabs(moddx) < fabs(mindx) )
	  mindx = moddx;
	if( fabs(moddy) < fabs(mindy) )
	  mindy = moddy;

      } // tri

      modmindxHisto.Fill( mindx );
      modmindyHisto.Fill( mindy );

      if( fabs( mindx ) < 0.15 && fabs( mindy ) < 0.1 ) {
	modxlkHisto.Fill( modx );
	modylkHisto.Fill( mody );
      }

    } // Mod

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // triplets vs MOD and DUT:

    int nmdm = 0;
    int ntrimod = 0;

    double xcutMOD = 0.15;
    double ycutMOD = 0.10;

    int nmtd = 0;
    int ntrilk = 0;
    int nsix = 0;

    double xcutDUT = 0.150; // 100 um
    double ycutDUT = 0.100; //  25 um
    if( rot90 ) {
      xcutDUT = 0.100; //  25 um
      ycutDUT = 0.150;
    }
    if( fifty ) {
      xcutDUT = 0.100;
      ycutDUT = 0.100;
    }
    if( fabs(DUTturn) > 33 ) // shallow
      xcutDUT = 0.300;

    for( unsigned int iA = 0; iA < triplets.size(); ++iA ) { // iA = upstream

      double xmA = triplets[iA].xm;
      double ymA = triplets[iA].ym;
      double zmA = triplets[iA].zm;
      double sxA = triplets[iA].sx;
      double syA = triplets[iA].sy;

      double zB = MODz - zmA; // z MOD from mid of triplet
      double xB = xmA + sxA * zB; // triplet impact point on MOD
      double yB = ymA + syA * zB;

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // tri vs tri: isolation at MOD

      double ttdminmod = 99.9;

      for( unsigned int jj = 0; jj < triplets.size(); ++jj ) {

	if( jj == iA ) continue;

	double xmj = triplets[jj].xm;
	double ymj = triplets[jj].ym;
	double sxj = triplets[jj].sx;
	double syj = triplets[jj].sy;

	double dz = MODz - triplets[jj].zm;
	double xj = xmj + sxj * dz; // triplet impact point on MOD
	double yj = ymj + syj * dz;

	double dx = xB - xj;
	double dy = yB - yj;
	double dd = sqrt( dx*dx + dy*dy );
	if( dd < ttdminmod )
	  ttdminmod = dd;

      } // jj

      ttdminmod1Histo.Fill( ttdminmod );
      ttdminmod2Histo.Fill( ttdminmod );
      triplets[iA].ttdmin = ttdminmod;

      // intersect inclined track with tilted MOD plane:

      double zd = (Nzm*zB - Nym*ymA - Nxm*xmA) / (Nxm*sxA + Nym*syA + Nzm); // from zmB
      double yd = ymA + syA * zd;
      double xd = xmA + sxA * zd;

      trixmHisto.Fill( xd );
      triymHisto.Fill( yd );

      double dzd = zd + zmA - MODz; // from MOD z0 [-8,8] mm

      // transform into MOD system: (passive).
      // large rotations don't commute: careful with order

      double x1m = com*xd - som*dzd; // turn o
      double y1m = yd;
      double z1m = som*xd + com*dzd;

      double x2m = x1m;
      double y2m = cam*y1m + sam*z1m; // tilt a

      double x3m = cfm*x2m + sfm*y2m; // rot
      double y3m =-sfm*x2m + cfm*y2m;

      double x4m =-x3m + MODalignx; // shift to mid
      double y4m = y3m + MODaligny; // invert y, shift to mid

      double xmodm = fmod( 36.000 + x4m, 0.3 ); // [0,0.3] mm, 2 pixel wide
      double ymodm = fmod(  9.000 + y4m, 0.2 ); // [0,0.2] mm

      if( ldbmod )
	cout << "                                   " << x4m << "  " << y4m;

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // triplet vs MOD clusters:

      bool ltrimod = 0;

      for( vector<cluster>::iterator c = cl[iMOD].begin(); c != cl[iMOD].end(); ++c ) {

	double ccol = c->col;
	double crow = c->row;
	double modx = ( ccol + 0.5 - nx[iMOD]/2 ) * ptchx[iMOD]; // -33..33 mm
	double mody = ( crow + 0.5 - ny[iMOD]/2 ) * ptchy[iMOD]; // -8..8 mm
	double q = c->signal;
	double q0 = q*normm;

	// residuals for pre-alignment:

	modsxaHisto.Fill( modx + x3m );
	moddxaHisto.Fill( modx - x3m );

	modsyaHisto.Fill( mody + y3m );
	moddyaHisto.Fill( mody - y3m );

	double moddx = modx - x4m;
	double moddy = mody - y4m;

	moddxHisto.Fill( moddx );
	moddyHisto.Fill( moddy );

	if( fabs( moddx ) < xcutMOD ) {

	  moddycHisto.Fill( moddy );
	  moddyvsx.Fill( -x3m, moddy ); // for rot
	  moddyvsy.Fill( y2m, moddy ); // for tilt
	  modmadyvsy.Fill( y2m, fabs(moddy) );
	  moddyvsty.Fill( syA, moddy );
	  moddyvst5.Fill( evsec, moddy );

	}

	if( fabs( moddy ) < ycutMOD ) {

	  moddxcHisto.Fill( moddx );
	  moddxvsx.Fill( -x1m, moddx ); // for turn
	  moddxvsy.Fill( y3m, moddx ); // for rot
	  modmadxvsy.Fill( y3m, fabs(moddx) );
	  moddxvstx.Fill( sxA, moddx );
	  moddxvst5.Fill( evsec, moddx );

	}

	if( fabs( moddx ) < xcutMOD &&
	    fabs( moddy ) < ycutMOD ) {

	  modnpxHisto.Fill( c->size );
	  modqHisto.Fill( q );
	  modq0Histo.Fill( q0 );
	  modnpxvsxmym->Fill( xmodm*1E3, ymodm*1E3, c->size );

	  modcollkHisto.Fill( ccol );
	  modrowlkHisto.Fill( crow );
	  trixmlkHisto.Fill( xd ); // telescope coordinates
	  triymlkHisto.Fill( yd );

	  //if( crow > 80 ) cout << "B link " << crow << " ev " << iev << endl;
	  //if( crow < 80 ) cout << "A link " << crow << " ev " << iev << endl;
	  //if( iev > 5100200 && crow < 80 ) cout << "link " << crow << " ev " << iev << endl; // A

	  triplets[iA].lk = 1;
	  nmdm = 1; // we have a MOD-triplet match in this event
	  ++ntrimod;
	  ltrimod = 1;

	  if( ldbmod )
	    cout << " match";

	  for( vector<pixel>::iterator px = c->vpix.begin(); px != c->vpix.end(); ++px ) {
	    tripxfrmlkHisto.Fill( px->frm ); // 0 or 1
	    tripxpivlkHisto.Fill( px->pivot ); // 0 or 1, cannot reduce pile-up
	  }

	} // MOD link x and y

      } // MOD

      if( ldbmod )
	cout << endl;

      if( ltrimod ) {
	trixmodHisto.Fill( xd );
	triymodHisto.Fill( yd );
      }

      double zA = DUTz - zmA; // z DUT from mid of triplet
      //double xA = xmA + sxA * zA; // triplet impact point on DUT
      //double yA = ymA + syA * zA;

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // intersect inclined track with turned DUT plane:

      double zc = (Nz*zA - Ny*ymA - Nx*xmA) / (Nx*sxA + Ny*syA + Nz); // from zmA
      double yc = ymA + syA * zc;
      double xc = xmA + sxA * zc;

      trixcHisto.Fill( xc );
      triycHisto.Fill( yc );
      trixycHisto->Fill( xc, yc );
      if( ltrimod ) {
	trixcmodHisto.Fill( xc );
	triycmodHisto.Fill( yc );
	trixycmodHisto->Fill( xc, yc );
      }

      double dzc = zc + zmA - DUTz; // from DUT z0 [-8,8] mm

      // transform into DUT system: (passive).
      // large rotations don't commute: careful with order

      double x1 = co*xc - so*dzc; // turn o
      double y1 = yc;
      double z1 = so*xc + co*dzc;

      double x2 = x1;
      double y2 = ca*y1 + sa*z1; // tilt a

      double x3 = cf*x2 + sf*y2; // rot
      double y3 =-sf*x2 + cf*y2;

      double x4 = x3 + DUTalignx; // shift to mid
      double y4 = y3 + DUTaligny; // shift to mid

      double xmod = fmod( 9.000 + x4, 0.100 ); // [0,0.100] mm
      double ymod = fmod( 9.000 + y4, 0.100 ); // [0,0.100] mm
      double xmod2 = fmod( 9.000 + x4, 0.200 ); // [0,0.200] mm
      double ymod2 = fmod( 9.000 + y4, 0.200 ); // [0,0.200] mm
      double xmod5 = fmod( 9.000 + x4, 0.050 ); // [0,0.050] mm
      double xmod25 = fmod( 9.000 + x4, 0.025 ); // [0,0.025] mm
      double ymod5 = fmod( 9.000 + y4, 0.050 ); // [0,0.050] mm
      double xmod4 = fmod( 9.000 + x4, 0.400 ); // [0,0.400] mm RD53 8x8 core

      double dxbeam = x4 - xbeam;
      double dybeam = y4 - ybeam;
      double drbeam = sqrt( dxbeam*dxbeam + dybeam*dybeam );

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // tri vs tri: isolation at DUT

      double ttdmin = 99.9;

      for( unsigned int jj = 0; jj < triplets.size(); ++jj ) {

	if( jj == iA ) continue;

	double xmj = triplets[jj].xm;
	double ymj = triplets[jj].ym;
	double sxj = triplets[jj].sx;
	double syj = triplets[jj].sy;

	double dz = zc + zmA - triplets[jj].zm;
	double xj = xmj + sxj * dz; // triplet impact point on DUT
	double yj = ymj + syj * dz;

	double dx = xc - xj;
	double dy = yc - yj;
	double dd = sqrt( dx*dx + dy*dy );
	if( dd < ttdmin )
	  ttdmin = dd;

      } // jj

      ttdmin1Histo.Fill( ttdmin );
      ttdmin2Histo.Fill( ttdmin );
      triplets[iA].ttdmin = ttdmin;

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // match triplet and driplet:

      double sixcut = 0.1; // [mm]

      for( unsigned int jB = 0; jB < driplets.size(); ++jB ) { // j = B = downstream

	double xmB = driplets[jB].xm;
	double ymB = driplets[jB].ym;
	double zmB = driplets[jB].zm;
	double sxB = driplets[jB].sx;
	double syB = driplets[jB].sy;

	// driplet at z DUT:
	/*
	  double zB = DUTz - zmB; // z from mid of triplet to mid
	  double xB = xmB + sxB * zB; // triplet at mid
	  double yB = ymB + syB * zB;
	*/
	// driplet at DUT:

	double zd = zc + zmA - zmB; // z from mid of driplet to DUT intersect
	double xd = xmB + sxB * zd; // driplet at DUT
	double yd = ymB + syB * zd;

	// driplet - triplet:

	//double dx = xB - xA; // at DUT
	//double dy = yB - yA;
	double dx = xd - xc; // at DUT intersect
	double dy = yd - yc;
	//double dxy = sqrt( dx*dx + dy*dy );
	double dtx = sxB - sxA;
	double dty = syB - syA;
	double dtxy = sqrt( dtx*dtx + dty*dty );

	hsixdx.Fill( dx ); // for align fit
	hsixdy.Fill( dy ); // for align fit

	if( fabs(dy) < sixcut ) {

	  hsixdxc.Fill( dx );

	  sixdxvsx.Fill( x4, dx );
	  sixmadxvsx.Fill( x4, fabs(dx) );
	  sixdxvsy.Fill( yc, dx );
	  sixdxvstx.Fill( sxA, dx );
	  sixdxvsdtx.Fill( dtx, dx );
	  sixdxvst3.Fill( evsec, dx );
	  sixdxvst5.Fill( evsec, dx );
	  sixmadxvsy.Fill( y4, fabs(dx) );
	  sixmadxvstx.Fill( sxA, fabs(dx) );
	  sixmadxvsdtx.Fill( dtx, fabs(dx) ); // U-shape

	} // dy

	if( fabs(dx) < sixcut ) {

	  hsixdyc.Fill( dy );

	  sixdyvsx.Fill( x4, dy );
	  sixmadyvsx.Fill( x4, fabs(dy) );
	  sixdyvsy.Fill( y4, dy );
	  sixdyvsty.Fill( syA, dy );
	  sixdyvsdty.Fill( dty, dy );
	  sixdyvst3.Fill( evsec, dy );
	  sixdyvst5.Fill( evsec, dy );
	  sixmadyvsy.Fill( y4, fabs(dy) );
	  sixmadyvsty.Fill( syA, fabs(dy) );
	  sixmadyvsdty.Fill( dty, fabs(dy) ); // U-shape

	} // dx

	// match:

	if( fabs(dx) < sixcut && fabs(dy) < sixcut ) {

	  ++nsix;

	  // average driplet and triplet at DUT:

	  double xa = 0.5 * ( xd + xc );
	  double ya = 0.5 * ( yd + yc );

	  sixxyHisto->Fill( xa, ya );

	  // compare slopes:

	  hsixdtx.Fill( dtx );
	  hsixdty.Fill( dty ); // width: 0.3 mrad
	  sixdtvsxav.Fill( xa, dtxy );
	  sixdtvsyav.Fill( ya, dtxy );
	  sixdtvsxyav->Fill( xa, ya, dtxy );
	  sixdtvsx.Fill( x4, dtxy );
	  sixdtvsxy->Fill( x4, y4, dtxy );
	  sixdtvsxm.Fill( xmod*1E3, dtxy );
	  sixdtvsym.Fill( ymod*1E3, dtxy );
	  sixdtvsxmym->Fill( xmod*1E3, ymod*1E3, dtxy ); // gStyle->SetPalette(55) 105 107

	  // transform into DUT system: (passive)

	  double dzc = zc + zmA - DUTz; // from DUT z0 [-8,8] mm

	  double x5 = co*xa - so*dzc; // turn o
	  double y5 = ya;
	  double z5 = so*xa + co*dzc;

	  double x6 = x5;
	  double y6 = ca*y5 + sa*z5; // tilt a

	  double x7 = cf*x6 + sf*y6; // rot
	  double y7 =-sf*x6 + cf*y6;

	  double x8 = x7 + DUTalignx; // shift to mid
	  double y8 = y7 + DUTaligny;

	  // update xy with six-plane average if no cooling box material:

	  if( chip0 == 182 || chip0 == 211 || // fresh: no box
	      chip0 == 501 || chip0 == 504 ||
	      chip0 == 520 || chip0 == 524 || chip0 == 529 ||
	      chip0 == 531 || chip0 == 543 || chip0 == 550 ||
	      chip0 == 534 || chip0 == 535 || // Zh card test Nov 2019
	      chip0 == 331 || chip0 == 731 || // 3D 
	      chip0 == 719 || // 3D
	      chip0 == 577 || chip0 == 578 || // Zh card
	      chip0 == 599 || chip0 == 605 || chip0 == 606 ||
	      chip0 == 364 || chip0 == 3811 || // FBK
	      chip0 >= 60100 ) 
          {
	    x4 = x8;
	    y4 = y8;

	    xmod = fmod( 9.000 + x4, 0.100 ); // [0,0.100] mm
	    ymod = fmod( 9.000 + y4, 0.100 ); // [0,0.100] mm

	    xmod2 = fmod( 9.000 + x4, 0.200 ); // [0,0.200] mm
	    ymod2 = fmod( 9.000 + y4, 0.200 ); // [0,0.200] mm

	    xmod5 = fmod( 9.000 + x4, 0.050 ); // [0,0.050] mm
	    xmod25 = fmod( 9.000 + x4, 0.025 ); // [0,0.025] mm
	    ymod5 = fmod( 9.000 + y4, 0.050 ); // [0,0.050] mm

	    xmod4 = fmod( 9.000 + x4, 0.400 ); // [0,0.400] mm RD53 8x8 core

	    dxbeam = x4 - xbeam;
	    dybeam = y4 - ybeam;
	    drbeam = sqrt( dxbeam*dxbeam + dybeam*dybeam );
	  }

	} // six match

      } // driplets

      // Coordinates for the cell histograms
      hcell1.add_coordinates(x4*1E3, y4*1E3);
      hcell2.add_coordinates(x4*1E3, y4*1E3);
      hcell4.add_coordinates(x4*1E3, y4*1E3);

      dutxyHisto->Fill( x4, y4 );

      // from track x, y (at DUT) to sensor col, row:
      // for straight 50x50:

      int kcol = ( x4 / ptchx[iDUT] + 0.5*nx[iDUT] ); // straight 50x50
      int krow = ( y4 / ptchy[iDUT] + 0.5*ny[iDUT] );
      if( rot90 ) { // 25x100
	kcol = (-y4 / ptchx[iDUT] + 0.5*nx[iDUT] ); // 0..398 even
	krow = ( x4 / ptchy[iDUT] + 0.5*ny[iDUT] ); // 0..191 full
      }

      if( kcol < 0 ) kcol = 0;
      if( kcol >= nx[iDUT] ) kcol = nx[iDUT]-1;

      if( krow < 0 ) krow = 0;
      if( krow >= ny[iDUT] ) krow = ny[iDUT]-1;

      int lrow;
      if( krow%2 ) // odd
	lrow = krow-1; // for 100x25: row pairs
      else
	lrow = krow+1;

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // triplet vs DUT clusters:

      double pdmin = 19;
      vector<cluster>::iterator cnear;
      int nlk = 0;
      bool shrt = 0;
      double dxmin = 99;
      double dymin = 99;

      for( vector<cluster>::iterator c = cl[iDUT].begin(); c != cl[iDUT].end(); ++c ) {

	double ccol = c->col;
	double crow = c->row;

	double dutx = ( ccol + 0.5 - nx[iDUT]/2 ) * ptchx[iDUT]; // mm
	double duty = ( crow + 0.5 - ny[iDUT]/2 ) * ptchy[iDUT]; // mm

	if( rot90 ) {
	  dutx = ( crow + 0.5 - ny[iDUT]/2 ) * ptchy[iDUT]; // mm
	  duty = ( ccol + 0.5 - nx[iDUT]/2 ) * ptchx[iDUT]; // mm
	}

	dutxxHisto->Fill( x4, dutx );
	dutyyHisto->Fill( y4, duty );

	//if( ccol > 200 && ccol < 201 ) cout << "col " << ccol << "  " << dutx << endl;
	//if( crow > 95.5 && crow < 96.5 ) cout << "row " << crow << "  " << duty << endl;

	double Q = c->signal + 0.5; // bin center correction
	double Q0 = Q*norm;
	double Qx = exp( -Q0 / qwid ); // Moyal weighting

	int npx = c->size;

	// residuals for pre-alignment:

	dutdxaHisto.Fill( dutx - x3 );
	dutsxaHisto.Fill(-dutx - x3 );
 	dutdyaHisto.Fill( duty - y3 );
 	dutsyaHisto.Fill(-duty - y3 );

	double dutdx = dutx - x4;
	double dutdy = duty - y4;
	if( rot90 )
	  dutdy = -duty - y4; // HLL, HPK, FBK
	if( rot90 && chip0 > 790100 )
	  dutdx = -dutx - x4; // FBK Nov 2019

	dutdxHisto.Fill( dutdx );
	dutdyHisto.Fill( dutdy );

	if( fabs(dutdx) < fabs(dxmin) ) dxmin = dutdx;
	if( fabs(dutdy) < fabs(dymin) ) dymin = dutdy;

	if( fabs( dutdy ) < ycutDUT ) {

	  dutdxcHisto.Fill( dutdx );
	  dutresxHisto.Fill( dutdx );

	  if( xmod25*1e3 < 6.25 || xmod25 > 18.75 )
	    dutdxccHisto.Fill( dutdx ); // crack
	  else
	    dutdxcmHisto.Fill( dutdx ); // mid

	  lindxcHisto.Fill( dutdx );

	  if(      c->nrow == 1 )
	    dutdxc1Histo.Fill( dutdx );
	  else if( c->nrow == 2 )
	    dutdxc2Histo.Fill( dutdx );

	  if( Q0 > qL && Q0 < qR )
	    dutdxcqHisto.Fill( dutdx );

	  dutdxvsx.Fill( x4, dutdx ); // for turn
	  dutdxvsy.Fill( y4, dutdx ); // for rot
	  dutdxvstx.Fill( sxA, dutdx );
	  if( x4 < 0 )
	    dutdxvstx0.Fill( sxA, dutdx );
	  else
	    dutdxvstx1.Fill( sxA, dutdx );
	  dutdxvsxmym->Fill(xmod*1E3, ymod*1E3, dutdx);
	  dutdxvsxm.Fill( xmod*1E3, dutdx );
	  dutdxvst2.Fill( evsec, dutdx );
	  dutdxvst5.Fill( evsec, dutdx );

	  dutmadxvsx.Fill( x4, fabs(dutdx) );
	  dutmadxvsxm.Fill( xmod*1E3, fabs(dutdx) );
	  dutmadxvsxm5.Fill( xmod5*1E3, fabs(dutdx) );
	  dutmadxvstx.Fill( sxA*1e3, fabs(dutdx) );
	  dutmadxvsq.Fill( Q0, fabs(dutdx) );

	} // cut y

	if( fabs( dutdx ) < xcutDUT ) {

	  dutdy0Histo.Fill( dutdy );
	  dutdycHisto.Fill( dutdy );
	  dutresyHisto.Fill( dutdy );

	  dutmadyvsq.Fill( Q0, fabs(dutdy) );
	  if( Q0 > qL && Q0 < qR )
	    dutdycqHisto.Fill( dutdy );

	  dutdyvsx.Fill( x4, dutdy ); // for rot
	  dutdyvsy.Fill( y4, dutdy ); // for tilt
	  dutdyvsty.Fill( syA, dutdy );
	  dutdyvsxm.Fill( xmod*1E3, dutdy );
	  dutdyvsym.Fill( ymod*1E3, dutdy );
	  dutdyvsxmym->Fill( xmod*1E3, ymod*1E3, dutdy );
	  dutdyvst2.Fill( evsec, dutdy );
	  dutdyvst5.Fill( evsec, dutdy );

	  dutmadyvsx.Fill( x4, fabs(dutdy) );
	  dutmadyvsy.Fill( y4, fabs(dutdy) );
	  dutmadyvstx.Fill( sxA, fabs(dutdy) );
	  dutmadyvsty.Fill( syA, fabs(dutdy) );
	  dutmadyvsxm.Fill( xmod*1E3, fabs(dutdy) );
	  dutmadyvsym.Fill( ymod*1E3, fabs(dutdy) );
	  dutmadyvst.Fill( evsec, fabs(dutdy) );

	} // cut x

	// cut x and y: cluster on track

	if( fabs( dutdx ) < xcutDUT &&
	    fabs( dutdy ) < ycutDUT ) {

	  dutncolHisto.Fill( c->ncol );
	  if( y4 < 4.7 && y4 > -4.7 ) { // fiducial
	    dutncolvsx.Fill( x4, c->ncol );
	    if( x4 <  8  && x4 > -1.5 ) { // fiducial
	      dutncolfHisto.Fill( c->ncol );
	      if( c->ncol < 9 ) {
		trixy4shrHisto->Fill( x4, y4 );
		trix4shrHisto.Fill( x4 );
		triy4shrHisto.Fill( y4 );
		shrt = 1;
	      }
	    }
	  }

	  minbcHisto.Fill( c->minbc );
	  maxbcHisto.Fill( c->maxbc );
	  nbcHisto.Fill( c->maxbc - c->minbc + 1 );

	  linqHisto.Fill( Q );
	  linq0Histo.Fill( Q0 );
	  linqxvsx.Fill( x4, Qx );
	  linqxvsy.Fill( y4, Qx );
	  linqxvsxy->Fill( x4, y4, Qx );
	  linqxvsr.Fill( drbeam, Qx );

	  linnpxHisto.Fill( npx );
	  linncolHisto.Fill( c->ncol );
	  linnrowHisto.Fill( c->nrow );

	  if( c->ncol == 1 ) {
	    linnrow1Histo.Fill( c->nrow );
	    vector<pixel>::iterator px = c->vpix.begin(); // 1st pixel
	    if( px->row%2 )
	      linnrow1oddHisto.Fill( c->nrow );
	    else
	      linnrow1eveHisto.Fill( c->nrow ); // twice more, similar shape
	  }
	  
          if( c->minbc == c->maxbc ) {
	    linnpxfHisto.Fill( npx );
	    linncolfHisto.Fill( c->ncol );
	    linnrowfHisto.Fill( c->nrow );
	  }
	  else {
	    linnpxffHisto.Fill( npx );
	    linncolffHisto.Fill( c->ncol );
	    linnrowffHisto.Fill( c->nrow );
	  }

	  linnpxvsxmym->Fill( xmod*1E3, ymod*1E3, npx );
          hcell1.fill_clsize(npx);
          hcell2.fill_clsize(npx);
          hcell4.fill_clsize(npx);
	  linncolvsxm.Fill( xmod*1E3, c->ncol );
	  linncolvsym.Fill( ymod*1E3, c->ncol );

	  linnrowvsxmym->Fill( xmod*1E3, ymod*1E3, c->nrow );
	  linnrowvsxm.Fill( xmod*1E3, c->nrow );
	  linnrowvsxm5.Fill( xmod5*1E3, c->nrow );
	  if( !fifty && !rot90 && xmod > 0.010 && xmod < 0.090 )
	    linnrowvsym.Fill( ymod*1E3, c->nrow );
	  if( fifty )
	    linnrowvsym.Fill( ymod*1E3, c->nrow );

	  linqxvsxmym->Fill( xmod*1E3, ymod*1E3, Qx );
          hcell1.fill_charge(Qx);
          hcell2.fill_charge(Qx);
          hcell4.fill_charge(Qx);
	  linqxvsxm.Fill( xmod*1E3, Qx );
	  linqxvsxm2.Fill( xmod2*1E3, Qx );
	  linqxvsxm5.Fill( xmod5*1E3, Qx );
	  linqxvsym.Fill( ymod*1E3, Qx );

	  trixyclkHisto->Fill( xc, yc );
	  trixclkHisto.Fill( xc );
	  triyclkHisto.Fill( yc );
	  dutlkxmymHisto->Fill( xmod*1E3, ymod*1E3 );
	  dutlkcolHisto.Fill( ccol );
	  dutlkrowHisto.Fill( crow );

	  // pixels in the cluster:

	  int colmin = 999;
	  int colmax = -1;
	  int rowmin = 999;
	  int rowmax = -1;
	  vector <int> colsz(400); // initialized as zero
	  vector <int> colq(400); // initialized as zero
	  vector <int> rowq(384); // initialized as zero
	  int totmax = 0;
	  vector<pixel>::iterator pxmax;

	  bool ldb = 0;
	  if( ldb )
	    cout << iev << endl;

	  for( vector<pixel>::iterator px = c->vpix.begin(); px != c->vpix.end(); ++px ) {

	    int icol = px->col; // sensor
	    int irow = px->row;
	    int itot = px->tot;
	    int ifrm = px->frm;

	    if( ldb )
	      cout
		<< "  " << icol
		<< "  " << irow
		<< "  " << itot
		<< "  " << ifrm // 7-8 or 9-10
		<< endl;

	    dutpxcolHisto.Fill( icol + 0.5 );
	    dutpxrowHisto.Fill( irow + 0.5 );

	    if( icol < colmin ) colmin = icol;
	    if( icol > colmax ) colmax = icol;
	    if( irow < rowmin ) rowmin = irow;
	    if( irow > rowmax ) rowmax = irow;
	    if( itot > totmax ) {
	      totmax = itot;
	      pxmax = px;
	    }

	    ++colsz[icol];
	    colq[icol] += itot;
	    rowq[irow] += itot;

	    if( px == c->vpix.begin() ) { // 1st px in readout order: cross talk ?

	      dutpxrow1Histo.Fill( irow + 0.5 ); // even-odd pattern in 25x100
	      dutpxcol1Histo.Fill( icol + 0.5 ); // 8-pattern in 50x50

	      if( c->ncol == 1 ) {

		dutpxrow11Histo.Fill( irow + 0.5 ); // 70% even

		if( c->minbc == c->maxbc ) { // 1-frame clusters: pixels read out sequentially

		  dutpxrow111Histo.Fill( irow + 0.5 ); // 75% even

		  if( c->nrow == 2 )
		    dutpxrow1112Histo.Fill( irow + 0.5 ); // 100% even as it should be for 25x100

		  if( irow%2 )
		    dutpxcol111oddHisto.Fill( icol + 0.5 );
		  else
		    dutpxcol111eveHisto.Fill( icol + 0.5 ); // same shape

		} // frm

	      } // ncol

	      if( c->nrow == 1 )
		dutpxcol11Histo.Fill( icol + 0.5 );

	    } // 1st px

	    dutpxqHisto.Fill( itot );

	    if( drbeam < rbeam )
	      dutpxqbeamHisto.Fill( itot );

	    if( c->ncol == 1 )
	      dutpxq1Histo.Fill( itot );
	    else
	      dutpxq2Histo.Fill( itot ); // smaller

	    if( fifty ) {
	      if( icol%2 )
		dutpxqoddHisto.Fill( itot );
	      else
		dutpxqeveHisto.Fill( itot );
	    }
	    else {
	      if( irow%2 )
		dutpxqoddHisto.Fill( itot );
	      else
		dutpxqeveHisto.Fill( itot ); // similar
	    }

	    dutpxbclkHisto.Fill( ifrm );
	    if( c->size == 1 )
	      dutpxbclk1pHisto.Fill( ifrm );
	    if( c->ncol == 1 )
	      dutpxbclk1cHisto.Fill( ifrm );
	    if( c->nrow == 1 )
	      dutpxbclk1rHisto.Fill( ifrm );

	    dutpxqbcHisto[ifrm].Fill( itot ); // charge per BC
	    dutpxbcvsq.Fill( itot+0.5, ifrm );

	    dutpxbcvsx.Fill( icol+0.5, ifrm );
	    dutpxbcvsy.Fill( irow+0.5, ifrm );
	    dutpxbcvsxy->Fill( icol+0.5, irow+0.5, ifrm );

	    dutpxbcvst5.Fill( evsec, ifrm ); // long run 34135: not stable

	    // next px:

	    vector<pixel>::iterator nx = px;
	    ++nx;

	    if( nx != c->vpix.end() ) {

	      dutpxdcolHisto.Fill( nx->col - icol );
	      dutpxdrowHisto.Fill( nx->row - irow );
	      dutpxdfrmHisto.Fill( nx->frm - ifrm ); // mostly zero or +1
	      dutpxqqHisto.Fill( itot, nx->tot );

	      if( fifty ) {

		if( nx->row == irow && // same row
		    ( nx->col == icol + 1 || nx->col + 1 == icol ) ) { // adjacent cols

		  double eta = ( nx->tot - itot ) / double( nx->tot + itot );
		  if( nx->col < icol ) eta *= -1; // symmetrize
		  dutetaHisto.Fill( eta );
		  if( c->ncol == 2 )
		    duteta2Histo.Fill( eta ); // spiky

		} // topo

	      } // 50

	      else {

		if( nx->col == icol && // same column
		    ( nx->row == irow + 1 || nx->row + 1 == irow ) ) { // adjacent rows

		  double eta = ( nx->tot - itot ) / double( nx->tot + itot );
		  if( nx->row < irow ) eta *= -1; // symmetrize
		  dutetaHisto.Fill( eta );
		  if( c->nrow == 2 )
		    duteta2Histo.Fill( eta ); // spiky

		  if( nx->row%2 ) { // odd
		    dutpxqq0Histo.Fill( itot, nx->tot );
		    dutpxq1oddHisto.Fill( itot );
		    dutpxq2oddHisto.Fill( nx->tot );
		  }
		  else {
		    dutpxqq1Histo.Fill( itot, nx->tot );
		    dutpxq1eveHisto.Fill( itot );
		    dutpxq2eveHisto.Fill( nx->tot );
		  }

		} // topo

	      } // 25

	    } // nx

	  } // loop px

	  linpxqmaxHisto.Fill( totmax );
	  if( c->size > 1 )
	    linpxqmax2Histo.Fill( totmax );
	  int maxbc = pxmax->frm;

	  for( vector<pixel>::iterator px = c->vpix.begin(); px != c->vpix.end(); ++px ) {
	    if( px == pxmax ) continue;
	    linpxq2ndHisto.Fill( px->tot );
	    lintwvsq.Fill( px->tot, px->frm - maxbc ); // timewalk
	    if( pxmax->tot == 15 ) // overflow = fast
	      lintwvsq15.Fill( px->tot, px->frm - maxbc ); // timewalk
	  }

	  for( int icol = colmin+1; icol < colmax; ++icol ) {
	    dutcolszHisto.Fill( colsz[icol] );
	    dutcolszvsc.Fill( icol-colmin, colsz[icol] );
	    dutcolqHisto.Fill( colq[icol] );
	    dutcolqvsc.Fill( icol-colmin, colq[icol] );
	  }

	  // row-charge vs distance to track:

	  for( int irow = krow-2; irow <= krow+2; ++irow ) {

	    if( irow <   0 ) continue;
	    if( irow > 383 ) break;

	    double py = ( irow + 0.5 - ny[iDUT]/2 ) * ptchy[iDUT]; // mm
	    double pdy = py - y4;
	    dutpdyHisto.Fill( pdy*1E3, rowq.at(irow) );

	    if( krow%2 )
	      dutrowqvsdy1.Fill( pdy*1E3, rowq.at(irow) );
	    else
	      dutrowqvsdy0.Fill( pdy*1E3, rowq.at(irow) );

	    dutrowvsdy.Fill( pdy*1E3, (rowq.at(irow)>0?1:0) );
	    if( krow%2 )
	      dutrowvsdy1.Fill( pdy*1E3, (rowq.at(irow)>0?1:0) );
	    else
	      dutrowvsdy0.Fill( pdy*1E3, (rowq.at(irow)>0?1:0) );

	  }

	  // correlation maps:

	  if( krow > 0 ) { // p = k-1, q = k
	    dutrowpqmap->Fill( kcol + 0.5, krow + 0.5, rowq.at(krow-1) * rowq.at(krow) );
	    dutrowppmap->Fill( kcol + 0.5, krow + 0.5, rowq.at(krow-1) * rowq.at(krow-1) );
	    dutrowqqmap->Fill( kcol + 0.5, krow + 0.5, rowq.at(krow) * rowq.at(krow) );
	  }
	  if( kcol > 0 ) { // p = k-1, q = k
	    dutcolpqmap->Fill( kcol + 0.5, krow + 0.5, colq.at(kcol-1) * colq.at(kcol) );
	    dutcolppmap->Fill( kcol + 0.5, krow + 0.5, colq.at(kcol-1) * colq.at(kcol-1) );
	    dutcolqqmap->Fill( kcol + 0.5, krow + 0.5, colq.at(kcol) * colq.at(kcol) );
	  }

	  if( c->ncol == 1 ) { // for 100x25

	    // seed row and pair:

	    linqseedHisto.Fill( rowq.at(krow) );
	    linqpairHisto.Fill( rowq.at(lrow) );

	    linqseedvsym.Fill( ymod5*1E3, rowq.at(krow) );
	    linqpairvsym.Fill( ymod5*1E3, rowq.at(lrow) );
	    linqymHisto.Fill( ymod5*1E3, -rowq.at(lrow) );
	    linqymHisto.Fill( ymod5*1E3,  rowq.at(krow) );

	    // correlation profile: p = l, q = k

	    linpqvsym.Fill( ymod5*1E3, rowq.at(lrow) * rowq.at(krow) );
	    linppvsym.Fill( ymod5*1E3, rowq.at(lrow) * rowq.at(lrow) );
	    linqqvsym.Fill( ymod5*1E3, rowq.at(krow) * rowq.at(krow) );

	    linrowminHisto.Fill( rowmin + 0.5 );
	    linrowmaxHisto.Fill( rowmax + 0.5 );
	    linrowmin01Histo.Fill( rowmin%2 );
	    linrowmax01Histo.Fill( rowmax%2 );

	    if( c->nrow == 2 ) {
	      linrowmin2Histo.Fill( rowmin + 0.5 ); // mostly even
	      linrowmax2Histo.Fill( rowmax + 0.5 );
	    }

	  } // 1-col

	  triplets[iA].lk = 1;
	  nmtd = 1; // we have a DUT-triplet match in this event
	  ++ntrilk;
	  ++nlk;

	} // DUT link x and y

	// for eff: nearest pixel

	if( c->minbc >= minBC && c->maxbc <= maxBC ) { // good timing

	  for( unsigned ipx = 0; ipx < c->vpix.size(); ++ipx ) {

	    double px = ( c->vpix[ipx].col + 0.5 - nx[iDUT]/2 ) * ptchx[iDUT]; // mm
	    double py = ( c->vpix[ipx].row + 0.5 - ny[iDUT]/2 ) * ptchy[iDUT]; // mm
	    if( rot90 ) { // HPK
	      px = ( c->vpix[ipx].row + 0.5 - ny[iDUT]/2 ) * ptchy[iDUT]; // mm
	      py =-( c->vpix[ipx].col + 0.5 - nx[iDUT]/2 ) * ptchx[iDUT]; // mm
	    }
	    if( rot90 && chip0 > 790100 ) { // FBK Nov 2019
	      px =-( c->vpix[ipx].row + 0.5 - ny[iDUT]/2 ) * ptchy[iDUT]; // mm
	      py = ( c->vpix[ipx].col + 0.5 - nx[iDUT]/2 ) * ptchx[iDUT]; // mm
	    }
	    dutpxxyHisto->Fill( px, py );
	    double pdx = px - x4; // triplet extrapol
	    double pdy = py - y4;
	    double pdxy = sqrt( pdx*pdx + pdy*pdy );
	    if( pdxy < pdmin ) {
	      pdmin = pdxy; // [mm]
	      cnear = c;
	    }

	  } // pix

	} // BC

      } // loop DUT clusters

      int nm[99] = {0};
      for( int iw = 1; iw < 99; ++iw )
	if( pdmin < iw*0.010 ) // 10 um bins
	  nm[iw] = 1; // eff

      if( nm[49] )
	dutxylkHisto->Fill( x4, y4 ); // tracks with cluster link

      if( nm[20]
	  && y4 < 4.7 && y4 > -4.7
	  && x4 <  8  && x4 > -1.5 // Lin fiducial
	  ) {
	dutnlkHisto.Fill( nlk );
	dutdxminHisto.Fill( dxmin );
	dutdyminHisto.Fill( dymin );
      }
      if( shrt )
	dutnlkshrHisto.Fill( nlk ); // nlk = 2 => split clusters

      if( pdmin < 0.5 ) {

	dutqnearHisto.Fill( cnear->signal + 0.5 ); // nearest cluster

	for( unsigned ipx = 0; ipx < cnear->vpix.size(); ++ipx ) {
	  dutpxqnearHisto.Fill( cnear->vpix[ipx].tot );
	  dutpxbcnearHisto.Fill( cnear->vpix[ipx].frm );
	}

      } // good link

	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // DUT efficiency vs isolated MOD-linked fiducial tracks:

      if( ltrimod ) {

	bool fidx = 1;
	bool fidy = 1;

	// fiducial region: Lin section, straight (vertical) mount:

	if( y4 >  4.7 ) fidy = 0; // [mm] track on DUT
	if( y4 < -4.7 ) fidy = 0;

	if( x4 >  3.1 ) fidx = 0;
	if( x4 < -3.5 ) fidx = 0;

	if( run >= 37621 && run <= 37635 ) { // FBK Nov 2019
	  if( x4 >  3.1 ) fidx = 0; // Lin from 136
	  if( x4 < -3.1 ) fidx = 0;
	}

	if( rot90 ) { // rot90 Lin
	  fidx = 1;
	  if( x4 >  4.7 ) fidx = 0;
	  if( x4 < -4.7 ) fidx = 0;
	  fidy = 1;
	  if( y4 >  3.5 ) fidy = 0;
	  if( y4 < -3.1 ) fidy = 0;
	}

	// from track x, y (at DUT) to sensor col, row:

	int kcol = ( x4 / ptchx[iDUT] + 0.5*nx[iDUT] ); // straight 50x50
	int krow = ( y4 / ptchy[iDUT] + 0.5*ny[iDUT] );
	if( rot90 ) { // 25x100
	  kcol = (-y4 / ptchx[iDUT] + 0.5*nx[iDUT] ); // 0..398 even
	  krow = ( x4 / ptchy[iDUT] + 0.5*ny[iDUT] ); // 0..191 full
	}
	if( rot90 && chip0 > 790100 ) {
	  kcol = ( y4 / ptchx[iDUT] + 0.5*nx[iDUT] ); // 0..398 even
	  krow = (-x4 / ptchy[iDUT] + 0.5*ny[iDUT] ); // 0..191 full
	}

	//if( x4 > 0 && x4 < 0.05 ) cout << "x " << x4 << "  " << kcol << endl;
	//if( y4 > 0 && y4 < 0.05 ) cout << "y " << y4 << "  " << krow << endl;

	int kpx = kcol*384 + krow; // sens pix
	bool dead = 0;
	if( deadset.count(kpx) )
	  dead = 1;

	if( fidx && fidy && !dead )
	  effvsdmin.Fill( ttdminmod, nm[49] ); // triplet isolation at MOD

	if( ttdminmod > 0.6 ) {

	  sixxylkHisto->Fill( xc, yc );
	  if( nm[49] ) sixxyeffHisto->Fill( xc, yc );

	  effvsxy->Fill( x4, y4, nm[49] ); // map

	  if( dead )
	    effvsxydead->Fill( x4, y4, nm[49] ); // map
	  else
	    effvsxylive->Fill( x4, y4, nm[49] ); // map

	  if( fidy ) {
	    effvsx.Fill( x4, nm[49] );
	    effvsx2.Fill( x4, nm[49] );
	    effvsx4.Fill( x4, nm[49] );
	  }

	  if( fidx ) {
	    effvsy.Fill( y4, nm[49] );
	    if( dead )
	      effvsydead.Fill( y4, nm[49] ); // map
	    else
	      effvsylive.Fill( y4, nm[49] ); // map
	  }

	  if( fidx && fidy ) {

	    effvsr.Fill( drbeam, nm[49] );
	    effvsxmym->Fill( xmod*1E3, ymod*1E3, nm[49] );
	    hcell1.fill_eff(nm[49]);
	    hcell2.fill_eff(nm[49]);
	    hcell4.fill_eff(nm[49]);

	    // efficiency around dead pixels:

	    for( int idead : deadset ) { // C++11 range-based for loop over all dead px
	      int dcol = idead/384;
	      int drow = idead%384; // sensor-based
	      double px = ( dcol + 0.5 - nx[iDUT]/2 ) * ptchx[iDUT]; // mm
	      double py = ( drow + 0.5 - ny[iDUT]/2 ) * ptchy[iDUT]; // mm
	      double pdx = px - x4; // triplet extrapol
	      double pdy = py - y4;
	      if( rot90 ) // HPK
		pdy = -py - y4;
	      if( rot90 && chip0 > 790100 ) // FBK
		pdx = -px - x4;
	      effvsdxdydead->Fill( pdx, pdy, nm[49] ); // map
	    }

	    if( dead )
	      effvsxmymdead->Fill( xmod*1E3, ymod*1E3, nm[49] );

	    else {

	      effvsev1.Fill( iev, nm[49] );
	      effvsev2.Fill( iev, nm[49] );
	      effvst1.Fill( evsec, nm[49] );
	      effvst2.Fill( evsec, nm[49] );
	      effvst3.Fill( evsec, nm[49] );
	      effvst5.Fill( evsec, nm[49] );
	      ++ntrck;
	      ngood += nm[49];

	      dutpdminHisto.Fill( pdmin );

	      for( int iw = 1; iw < 99; ++iw )
		effvsdxy.Fill( iw*0.010-0.001, nm[iw] );

	      effvsntri.Fill( triplets.size(), nm[49] ); // flat
	      effvsxmymlive->Fill( xmod*1E3, ymod*1E3, nm[49] );
	      effvsxmym2live->Fill( xmod2*1E3, ymod2*1E3, nm[49] );
	      effvsxm.Fill( xmod*1E3, nm[49] ); // bias dot
	      effvsxm2.Fill( xmod2*1E3, nm[49] ); // bias dot
	      effvsxm4.Fill( xmod4*1E3, nm[49] ); // bias dot
	      effvsym.Fill( ymod*1E3, nm[49] ); // bias dot
	      effvstx.Fill( sxA, nm[49] );
	      effvsty.Fill( syA, nm[49] );

	    } // live

	  } // fid

	} // iso triplet

      } // tri-mod link

    } // loop triplets

    if( nmdm )
      ++nmodlk;

    modlkvst1.Fill( evsec, nmdm ); // MOD yield vs time
    modlkvst3.Fill( evsec, nmdm );
    modlkvst5.Fill( evsec, nmdm );
    modlkvsev.Fill( iev, nmdm );
    modlkvsev9.Fill( iev, nmdm );
    modlkvsev1.Fill( iev, nmdm );
    modlkvsev2.Fill( iev, nmdm );
    ntrimodHisto.Fill( ntrimod );

    dutlkvst1.Fill( evsec, nmtd ); // DUT yield vs time
    dutlkvst3.Fill( evsec, nmtd );
    dutlkvst5.Fill( evsec, nmtd );
    ntrilkHisto.Fill( ntrilk ); // DUT links

    ++iev;

  } while( reader->NextEvent() && iev < lev );

  delete reader;

  cout << "done after " << iev << " events" << endl;
  histoFile.Write();
  //histoFile->Close();

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // MOD alignment:

  if( ldbmod == 0 && moddxaHisto.GetEntries() > 9999 ) {

    double newMODalignx = MODalignx;
    double newMODaligny = MODaligny;

    if( moddxaHisto.GetMaximum() > modsxaHisto.GetMaximum() ) {
      cout << endl << moddxaHisto.GetTitle()
	   << " bin " << moddxaHisto.GetBinWidth(1)
	   << endl;
      TF1 * fgp0 = new TF1( "fgp0", "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]", -1, 1 );
      double xpk = moddxaHisto.GetBinCenter( moddxaHisto.GetMaximumBin() );
      fgp0->SetParameter( 0, moddxaHisto.GetMaximum() ); // amplitude
      fgp0->SetParameter( 1, xpk );
      fgp0->SetParameter( 2, moddxaHisto.GetBinWidth(1) ); // sigma
      fgp0->SetParameter( 3, moddxaHisto.GetBinContent( moddxaHisto.FindBin(xpk-1) ) ); // BG
      moddxaHisto.Fit( "fgp0", "q", "", xpk-1, xpk+1 );
      cout << "Fit Gauss + BG:"
	   << endl << "  A " << fgp0->GetParameter(0)
	   << endl << "mid " << fgp0->GetParameter(1)
	   << endl << "sig " << fgp0->GetParameter(2)
	   << endl << " BG " << fgp0->GetParameter(3)
	   << endl;
      newMODalignx = fgp0->GetParameter(1);
      delete fgp0;
    }
    else {
      cout << endl << modsxaHisto.GetTitle()
	   << " bin " << modsxaHisto.GetBinWidth(1)
	   << endl;
      TF1 * fgp0 = new TF1( "fgp0", "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]", -1, 1 );
      double xpk = modsxaHisto.GetBinCenter( modsxaHisto.GetMaximumBin() );
      fgp0->SetParameter( 0, modsxaHisto.GetMaximum() ); // amplitude
      fgp0->SetParameter( 1, xpk );
      fgp0->SetParameter( 2, modsxaHisto.GetBinWidth(1) ); // sigma
      fgp0->SetParameter( 3, modsxaHisto.GetBinContent( modsxaHisto.FindBin(xpk-1) ) ); // BG
      modsxaHisto.Fit( "fgp0", "q", "", xpk-1, xpk+1  );
      cout << "Fit Gauss + BG:"
	   << endl << "  A " << fgp0->GetParameter(0)
	   << endl << "mid " << fgp0->GetParameter(1)
	   << endl << "sig " << fgp0->GetParameter(2)
	   << endl << " BG " << fgp0->GetParameter(3)
	   << endl;
      newMODalignx = fgp0->GetParameter(1);
      delete fgp0;
    }

    if( moddyaHisto.GetMaximum() > modsyaHisto.GetMaximum() ) {
      cout << endl << moddyaHisto.GetTitle()
	   << " bin " << moddyaHisto.GetBinWidth(1)
	   << endl;
      TF1 * fgp0 = new TF1( "fgp0", "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]", -1, 1 );
      double xpk = moddyaHisto.GetBinCenter( moddyaHisto.GetMaximumBin() );
      fgp0->SetParameter( 0, moddyaHisto.GetMaximum() ); // amplitude
      fgp0->SetParameter( 1, xpk );
      fgp0->SetParameter( 2, moddyaHisto.GetBinWidth(1) ); // sigma
      fgp0->SetParameter( 3, moddyaHisto.GetBinContent( moddyaHisto.FindBin(xpk-1) ) ); // BG
      moddyaHisto.Fit( "fgp0", "q", "", xpk-1, xpk+1 );
      cout << "Fit Gauss + BG:"
	   << endl << "  A " << fgp0->GetParameter(0)
	   << endl << "mid " << fgp0->GetParameter(1)
	   << endl << "sig " << fgp0->GetParameter(2)
	   << endl << " BG " << fgp0->GetParameter(3)
	   << endl;
      newMODaligny = fgp0->GetParameter(1);
      delete fgp0;
    }
    else {
      cout << endl << modsyaHisto.GetTitle()
	   << " bin " << modsyaHisto.GetBinWidth(1)
	   << endl;
      TF1 * fgp0 = new TF1( "fgp0", "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]", -1, 1 );
      double xpk = modsyaHisto.GetBinCenter( modsyaHisto.GetMaximumBin() );
      fgp0->SetParameter( 0, modsyaHisto.GetMaximum() ); // amplitude
      fgp0->SetParameter( 1, xpk );
      fgp0->SetParameter( 2, modsyaHisto.GetBinWidth(1) ); // sigma
      fgp0->SetParameter( 3, modsyaHisto.GetBinContent( modsyaHisto.FindBin(xpk-1) ) ); // BG
      modsyaHisto.Fit( "fgp0", "q", "", xpk-1, xpk+1 );
      cout << "Fit Gauss + BG:"
	   << endl << "  A " << fgp0->GetParameter(0)
	   << endl << "mid " << fgp0->GetParameter(1)
	   << endl << "sig " << fgp0->GetParameter(2)
	   << endl << " BG " << fgp0->GetParameter(3)
	   << endl;
      newMODaligny = fgp0->GetParameter(1);
      delete fgp0;
    }

    cout << endl << "coarse MODalign x changed by " << newMODalignx - MODalignx << " mm" << endl;

    cout << endl << "coarse MODalign y changed by " << newMODaligny - MODaligny << " mm" << endl;

    // finer alignment:

    if( MODaligniteration > 0 && fabs( newMODalignx - MODalignx ) < 0.1 ) {

      cout << endl << moddxcHisto.GetTitle()
	   << " bin " << moddxcHisto.GetBinWidth(1)
	   << endl;
      TF1 * fgp0 = new TF1( "fgp0", "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]", -1, 1 );
      fgp0->SetParameter( 0, moddxcHisto.GetMaximum() ); // amplitude
      fgp0->SetParameter( 1, moddxcHisto.GetBinCenter( moddxcHisto.GetMaximumBin() ) );
      fgp0->SetParameter( 2, 8*moddxcHisto.GetBinWidth(1) ); // sigma
      fgp0->SetParameter( 3, moddxcHisto.GetBinContent(1) ); // BG
      moddxcHisto.Fit( "fgp0", "q" );
      cout << "Fit Gauss + BG:"
	   << endl << "  A " << fgp0->GetParameter(0)
	   << endl << "mid " << fgp0->GetParameter(1)
	   << endl << "sig " << fgp0->GetParameter(2)
	   << endl << " BG " << fgp0->GetParameter(3)
	   << endl;
      newMODalignx = MODalignx + fgp0->GetParameter(1);
      delete fgp0;

      // dxvsx -> turn:

      if( fabs(som) > 0.01 &&
	  moddxvsx.GetEntries() > 999
	  ) {

	double x0 = -midx[iMOD]+0.2; // fit range
	for( int ix = 1; ix < moddxvsx.GetNbinsX(); ++ix ) {
	  if( moddxvsx.GetBinEntries( ix ) > 11 ) {
	    x0 = moddxvsx.GetBinLowEdge(ix) + 2*moddxvsx.GetBinWidth(ix);
	    break;
	  }
	}

	double x9 = midx[iMOD]-0.2; // [mm] full range
	for( int ix = moddxvsx.GetNbinsX(); ix > 0; --ix ) {
	  if( moddxvsx.GetBinEntries( ix ) > 11 ) {
	    x9 = moddxvsx.GetBinLowEdge(ix)-moddxvsx.GetBinWidth(ix);
	    break;
	  }
	}

	moddxvsx.Fit( "pol1", "q", "", x0, x9 );

	TF1 * fdxvsx = moddxvsx.GetFunction( "pol1" );
	cout << endl << moddxvsx.GetTitle()
	     << ": slope " << fdxvsx->GetParameter(1)
	     << ", extra turn " << fdxvsx->GetParameter(1)/wt/som
	     << " deg"
	     << endl;
	MODturn += fdxvsx->GetParameter(1)/wt/som; // [deg] min 0.6 deg
	//delete fdxvsx;
      }

    } // finer x

    if( MODaligniteration > 0 && fabs( newMODaligny - MODaligny ) < 0.1 ) {

      cout << endl << moddycHisto.GetTitle()
	   << " bin " << moddycHisto.GetBinWidth(1)
	   << endl;
      TF1 * fgp0 = new TF1( "fgp0", "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]", -1, 1 );
      fgp0->SetParameter( 0, moddycHisto.GetMaximum() ); // amplitude
      fgp0->SetParameter( 1, moddycHisto.GetBinCenter( moddycHisto.GetMaximumBin() ) );
      fgp0->SetParameter( 2, 5*moddycHisto.GetBinWidth(1) ); // sigma
      fgp0->SetParameter( 3, moddycHisto.GetBinContent(1) ); // BG
      moddycHisto.Fit( "fgp0", "q" );
      cout << "Fit Gauss + BG:"
	   << endl << "  A " << fgp0->GetParameter(0)
	   << endl << "mid " << fgp0->GetParameter(1)
	   << endl << "sig " << fgp0->GetParameter(2)
	   << endl << " BG " << fgp0->GetParameter(3)
	   << endl;
      newMODaligny = MODaligny + fgp0->GetParameter(1);
      delete fgp0;

      // dyvsx -> rot

      if( moddyvsx.GetEntries() > 999 ) {

	double x0 = -midx[iMOD]+0.2; // fit range
	for( int ix = 1; ix < moddyvsx.GetNbinsX(); ++ix ) {
	  if( moddyvsx.GetBinEntries( ix ) > 11 ) {
	    x0 = moddyvsx.GetBinLowEdge(ix);
	    break;
	  }
	}

	double x9 = midx[iMOD]-0.2;
	for( int ix = moddyvsx.GetNbinsX(); ix > 0; --ix ){
	  if( moddyvsx.GetBinEntries( ix ) > 11 ) {
	    x9 = moddyvsx.GetBinLowEdge(ix)+moddyvsx.GetBinWidth(ix);
	    break;
	  }
	}

	moddyvsx.Fit( "pol1", "q", "", x0, x9 );

	TF1 * fdyvsx = moddyvsx.GetFunction( "pol1" );
	cout << endl << moddyvsx.GetTitle()
	     << ": extra rot " << fdyvsx->GetParameter(1)*1E3 << " mrad" << endl;
	MODrot += fdyvsx->GetParameter(1);
	//delete fdyvsx;

      }

      // dyvsy -> tilt:

      if( fabs( sam ) > 0.01 &&
	  moddyvsy.GetEntries() > 999
	  ) {

	double x0 = -midy[iMOD]+0.2; // fit range
	for( int ix = 1; ix < moddyvsy.GetNbinsX(); ++ix ){
	  if( moddyvsy.GetBinEntries( ix ) > 11 ) {
	    x0 = moddyvsy.GetBinLowEdge(ix);
	    break;
	  }
	}

	double x9 = midy[iMOD]-0.2;
	for( int ix = moddyvsy.GetNbinsX(); ix > 0; --ix ){
	  if( moddyvsy.GetBinEntries( ix ) > 11 ) {
	    x9 = moddyvsy.GetBinLowEdge(ix)+moddyvsy.GetBinWidth(ix);
	    break;
	  }
	}

	moddyvsy.Fit( "pol1", "q", "", x0, x9 );

	TF1 * fdyvsy = moddyvsy.GetFunction( "pol1" );
	cout << endl << moddyvsy.GetTitle()
	     << ": slope " << fdyvsy->GetParameter(1)
	     << ", extra tilt " << fdyvsy->GetParameter(1)/wt/sam
	     << " deg"
	     << endl;

	MODtilt += fdyvsy->GetParameter(1)/wt/sam; // [deg] min 0.6 deg
	//delete fdyvsy;
      }

      // dyvsty -> dz:

      if( moddyvsty.GetEntries() > 999 ) {

	double x0 = -0.002;
	for( int ix = 1; ix < moddyvsty.GetNbinsX(); ++ix ){
	  if( moddyvsty.GetBinEntries( ix ) > 11 ) {
	    x0 = moddyvsty.GetBinLowEdge(ix);
	    break;
	  }
	}

	double x9 = 0.002;
	for( int ix = moddyvsty.GetNbinsX(); ix > 0; --ix ){
	  if( moddyvsty.GetBinEntries( ix ) > 11 ) {
	    x9 = moddyvsty.GetBinLowEdge(ix)+moddyvsty.GetBinWidth(ix);
	    break;
	  }
	}

	moddyvsty.Fit( "pol1", "q", "", x0, x9 );

	TF1 * fdyvsty = moddyvsty.GetFunction( "pol1" );
	cout << endl << moddyvsty.GetTitle()
	     << ": z shift " << fdyvsty->GetParameter(1)
	     << " mm"
	     << endl;
	MODz += fdyvsty->GetParameter(1);
	//delete fdyvsty;
      }

    } // finer y

    cout << endl << "new MOD alignment:" << endl
	 << "  alignx " << newMODalignx << endl
	 << "  aligny " << newMODaligny << endl
	 << "  rot    " << MODrot << endl
	 << "  tilt   " << MODtilt << endl
	 << "  turn   " << MODturn << endl
	 << "  dz     " << MODz - zz[1] << endl
      ;

  cout << endl
       << "MOD yield " << 100*modlkvst3.GetMean(2) << "%"
       << " from " << iev << " events"
       << endl;

    cout << "update MOD alignment file? (y/n)" << endl;
    string ans{"y"};
    string YES{"y"};

    //cin >> ans;

    if( ans == YES ) {

      // write new MOD alignment:

      ofstream MODalignFile( MODalignFileName.str() );

      MODalignFile << "# MOD alignment for run " << run << endl;
      ++MODaligniteration;
      MODalignFile << "iteration " << MODaligniteration << endl;
      MODalignFile << "alignx " << newMODalignx << endl;
      MODalignFile << "aligny " << newMODaligny << endl;
      MODalignFile << "rot " << MODrot << endl;
      MODalignFile << "tilt " << MODtilt << endl;
      MODalignFile << "turn " << MODturn << endl;
      MODalignFile << "dz " << MODz - zz[1] << endl;

      MODalignFile.close();

      cout << endl << "wrote MOD alignment iteration " << MODaligniteration
	   << " to " << MODalignFileName.str() << endl
	;

    }

  } // MOD

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // DUT alignment:

  if( DUTaligniteration == 0 ) {

    if( dutdxaHisto.GetEntries() > 999 ) {

      if( dutdxaHisto.GetMaximum() > dutsxaHisto.GetMaximum() ) {

	cout << endl << dutdxaHisto.GetTitle()
	     << " bin " << dutdxaHisto.GetBinWidth(1)
	     << " peak " << dutdxaHisto.GetMaximum()
	     << " at " << dutdxaHisto.GetBinCenter( dutdxaHisto.GetMaximumBin() )
	     << endl;
	DUTalignx0 = dutdxaHisto.GetBinCenter( dutdxaHisto.GetMaximumBin() );
      }
      else {
	cout << endl << dutsxaHisto.GetTitle()
	     << " bin " << dutsxaHisto.GetBinWidth(1)
	     << " peak " << dutsxaHisto.GetMaximum()
	     << " at " << dutsxaHisto.GetBinCenter( dutsxaHisto.GetMaximumBin() )
	     << endl;
	DUTalignx0 = dutsxaHisto.GetBinCenter( dutsxaHisto.GetMaximumBin() );
      }

      // y:

      if( dutdyaHisto.GetMaximum() > dutsyaHisto.GetMaximum() ) {
	cout << endl << dutdyaHisto.GetTitle()
	     << " bin " << dutdyaHisto.GetBinWidth(1)
	     << " peak " << dutdyaHisto.GetMaximum()
	     << " at " << dutdyaHisto.GetBinCenter( dutdyaHisto.GetMaximumBin() )
	     << endl;
	DUTaligny0 = dutdyaHisto.GetBinCenter( dutdyaHisto.GetMaximumBin() );
      }
      else {
	cout << endl << dutsyaHisto.GetTitle()
	     << " bin " << dutsyaHisto.GetBinWidth(1)
	     << " peak " << dutsyaHisto.GetMaximum()
	     << " at " << dutsyaHisto.GetBinCenter( dutsyaHisto.GetMaximumBin() )
	     << endl;
	DUTaligny0 = dutsyaHisto.GetBinCenter( dutsyaHisto.GetMaximumBin() );
      }

    } // stat
    else
      cout << "not enough for coarse alignment" << endl;

  } // iteration 0

  // finer alignment:

  if( DUTaligniteration > 0 ) {

    cout << endl;

    if( dutdxcHisto.GetEntries() > 999 ) {

      cout << "finer x " << dutdxcHisto.GetTitle()
	   << " bin " << dutdxcHisto.GetBinWidth(1)
	   << endl;
      TF1 * fgp0 = new TF1( "fgp0", "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]", -1, 1 );
      fgp0->SetParameter( 0, dutdxcHisto.GetMaximum() ); // amplitude
      fgp0->SetParameter( 1, dutdxcHisto.GetBinCenter( dutdxcHisto.GetMaximumBin() ) );
      fgp0->SetParameter( 2, 8*dutdxcHisto.GetBinWidth(1) ); // sigma
      fgp0->SetParameter( 3, dutdxcHisto.GetBinContent(1) ); // BG
      dutdxcHisto.Fit( "fgp0", "q" );
      cout << "Fit Gauss + BG:"
	   << endl << "  A " << fgp0->GetParameter(0)
	   << endl << "mid " << fgp0->GetParameter(1)
	   << endl << "sig " << fgp0->GetParameter(2)
	   << endl << " BG " << fgp0->GetParameter(3)
	   << endl;
      DUTalignx0 += fgp0->GetParameter(1);

    }
    else
      cout << "not enough for fine x alignment" << endl;

    // y:

    cout << endl;

    if( dutdycHisto.GetEntries() > 999 ) {

      cout << "finer y " << dutdycHisto.GetTitle()
	   << " bin " << dutdycHisto.GetBinWidth(1)
	   << endl;
      TF1 * fgp0 = new TF1( "fgp0", "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]", -1, 1 );
      fgp0->SetParameter( 0, dutdycHisto.GetMaximum() ); // amplitude
      fgp0->SetParameter( 1, dutdycHisto.GetBinCenter( dutdycHisto.GetMaximumBin() ) );
      fgp0->SetParameter( 2, 5*dutdycHisto.GetBinWidth(1) ); // sigma
      fgp0->SetParameter( 3, dutdycHisto.GetBinContent(1) ); // BG
      dutdycHisto.Fit( "fgp0", "q" );
      cout << "Fit Gauss + BG:"
	   << endl << "  A " << fgp0->GetParameter(0)
	   << endl << "mid " << fgp0->GetParameter(1)
	   << endl << "sig " << fgp0->GetParameter(2)
	   << endl << " BG " << fgp0->GetParameter(3)
	   << endl;
      DUTaligny0 += fgp0->GetParameter(1);

    }
    else
      cout << "not enough for fine y alignment" << endl;

  } // iteration > 0

  // angles from projections:

  if( DUTaligniteration > 1 ) {

    // dxvsy => rot

    if( rot90 || fifty ) {

      if( dutdxvsy.GetEntries() > 999 ) {

	double x0 = -midx[iDUT]+0.2; // fit range
	for( int ix = 1; ix < dutdxvsy.GetNbinsX(); ++ix ) {
	  if( dutdxvsy.GetBinEntries( ix ) > 11 ) {
	    x0 = dutdxvsy.GetBinLowEdge(ix);
	    break;
	  }
	}
	double x9 = midx[iDUT]-0.2;
	for( int ix = dutdxvsy.GetNbinsX(); ix > 0; --ix ){
	  if( dutdxvsy.GetBinEntries( ix ) > 11 ) {
	    x9 = dutdxvsy.GetBinLowEdge(ix)+dutdxvsy.GetBinWidth(ix);
	    break;
	  }
	}
	dutdxvsy.Fit( "pol1", "q", "", x0, x9 );
	TF1 * fdxvsy = dutdxvsy.GetFunction( "pol1" );
	cout << endl
	     << "fit " << dutdxvsy.GetTitle()
	     << " from " << x0
	     << " to " << x9
	     << ": extra rot " << fdxvsy->GetParameter(1)*1E3 << " mrad" << endl;
	DUTrot += fdxvsy->GetParameter(1);

      }
      else
	cout << "not enough for rot in dutdxvsy" << endl;
    }
    else { // !rot90 && !fifty

      // dyvsx => rot

      if( dutdyvsx.GetEntries() > 999 ) {

	double x0 = -midx[iDUT]+0.2; // fit range
	for( int ix = 1; ix < dutdyvsx.GetNbinsX(); ++ix ) {
	  if( dutdyvsx.GetBinEntries( ix ) > 11 ) {
	    x0 = dutdyvsx.GetBinLowEdge(ix);
	    break;
	  }
	}

	double x9 = midx[iDUT]-0.2;
	for( int ix = dutdyvsx.GetNbinsX(); ix > 0; --ix ){
	  if( dutdyvsx.GetBinEntries( ix ) > 11 ) {
	    x9 = dutdyvsx.GetBinLowEdge(ix)+dutdyvsx.GetBinWidth(ix);
	    break;
	  }
	}

	dutdyvsx.Fit( "pol1", "q", "", x0, x9 );

	TF1 * fdyvsx = dutdyvsx.GetFunction( "pol1" );
	cout << endl
	     << "fit " << dutdyvsx.GetTitle()
	     << " from " << x0
	     << " to " << x9
	     << ": extra rot " << fdyvsx->GetParameter(1)*1E3 << " mrad" << endl;
	DUTrot -= fdyvsx->GetParameter(1);

      }
      else
	cout << "not enough for rot in dutdyvsx" << endl;

    }

    // dxvsx => turn:

    if( dutdxvsx.GetEntries() > 999 && fabs( DUTturn ) > 0.3 ) { // [deg]

      double x0 = -midx[iDUT]+0.2; // fit range
      for( int ix = 1; ix < dutdxvsx.GetNbinsX(); ++ix ) {
	if( dutdxvsx.GetBinEntries( ix ) > 11 ) {
	  x0 = dutdxvsx.GetBinLowEdge(ix) + 2*dutdxvsx.GetBinWidth(ix);
	  break;
	}
      }
      double x9 = midx[iDUT]-0.2; // [mm] full range
      for( int ix = dutdxvsx.GetNbinsX(); ix > 0; --ix ) {
	if( dutdxvsx.GetBinEntries( ix ) > 11 ) {
	  x9 = dutdxvsx.GetBinLowEdge(ix)-dutdxvsx.GetBinWidth(ix);
	  break;
	}
      }
      dutdxvsx.Fit( "pol1", "q", "", x0, x9 );
      TF1 * fdxvsx = dutdxvsx.GetFunction( "pol1" );
      cout << endl
	   << "fit " << dutdxvsx.GetTitle()
	   << " from " << x0
	   << " to " << x9
	   << ": slope " << fdxvsx->GetParameter(1)
	   << ", extra turn " << fdxvsx->GetParameter(1)/wt/so
	   << " deg"
	   << endl;
      DUTturn += fdxvsx->GetParameter(1)/wt/so; // [deg] min 0.6 deg

    } // turn x
    else
      cout << "not enough for turn in dutdxvsx" << endl;

    // dyvsy => tilt:

    if( fabs( DUTtilt ) > 0.3 ) { // [deg]

      if( dutdyvsy.GetEntries() > 999 ) {

	double x0 = -midy[iDUT]+0.2; // fit range
	for( int ix = 1; ix < dutdyvsy.GetNbinsX(); ++ix ){
	  if( dutdyvsy.GetBinEntries( ix ) > 11 ) {
	    x0 = dutdyvsy.GetBinLowEdge(ix);
	    break;
	  }
	}

	double x9 = midy[iDUT]-0.2;
	for( int ix = dutdyvsy.GetNbinsX(); ix > 0; --ix ){
	  if( dutdyvsy.GetBinEntries( ix ) > 11 ) {
	    x9 = dutdyvsy.GetBinLowEdge(ix)+dutdyvsy.GetBinWidth(ix);
	    break;
	  }
	}

	dutdyvsy.Fit( "pol1", "q", "", x0, x9 );

	TF1 * fdyvsy = dutdyvsy.GetFunction( "pol1" );
	cout << endl
	     << "fit " << dutdyvsy.GetTitle()
	     << " from " << x0
	     << " to " << x9
	     << ": slope " << fdyvsy->GetParameter(1)
	     << ", extra tilt " << fdyvsy->GetParameter(1)/wt/sa
	     << " deg"
	     << endl;

	DUTtilt += fdyvsy->GetParameter(1)/wt/sa; // sa might be neg

      }
      else
	cout << "not enough for tilt in dutdyvsy" << endl;

    } // tilt

    if( ( rot90 || fifty ) ) {

      // dxvstx => dz:

      if( dutdxvstx.GetEntries() > 999 ) {

	double x0 = -0.002;
	for( int ix = 1; ix < dutdxvstx.GetNbinsX(); ++ix ){
	  if( dutdxvstx.GetBinEntries( ix ) > 11 ) {
	    x0 = dutdxvstx.GetBinLowEdge(ix);
	    break;
	  }
	}

	double x9 = 0.002;
	for( int ix = dutdxvstx.GetNbinsX(); ix > 0; --ix ){
	  if( dutdxvstx.GetBinEntries( ix ) > 11 ) {
	    x9 = dutdxvstx.GetBinLowEdge(ix)+dutdxvstx.GetBinWidth(ix);
	    break;
	  }
	}

	dutdxvstx.Fit( "pol1", "q", "", x0, x9 );

	TF1 * fdxvstx = dutdxvstx.GetFunction( "pol1" );
	cout << endl << dutdxvstx.GetTitle()
	     << ": z shift " << fdxvstx->GetParameter(1)
	     << " mm"
	     << endl;
	DUTz += fdxvstx->GetParameter(1);

      }
      else
	cout << "not enough for z in dutdxvstx" << endl;
    }
    else { // !rot90 && !fifty

      // dyvsty => dz:

      if( dutdyvsty.GetEntries() > 999 ) {

	double x0 = -0.002;
	for( int ix = 1; ix < dutdyvsty.GetNbinsX(); ++ix ){
	  if( dutdyvsty.GetBinEntries( ix ) > 11 ) {
	    x0 = dutdyvsty.GetBinLowEdge(ix);
	    break;
	  }
	}

	double x9 = 0.002;
	for( int ix = dutdyvsty.GetNbinsX(); ix > 0; --ix ){
	  if( dutdyvsty.GetBinEntries( ix ) > 11 ) {
	    x9 = dutdyvsty.GetBinLowEdge(ix)+dutdyvsty.GetBinWidth(ix);
	    break;
	  }
	}

	dutdyvsty.Fit( "pol1", "q", "", x0, x9 );

	TF1 * fdyvsty = dutdyvsty.GetFunction( "pol1" );
	cout << endl << dutdyvsty.GetTitle()
	     << ": z shift " << fdyvsty->GetParameter(1)
	     << " mm"
	     << endl;
	DUTz += fdyvsty->GetParameter(1);

      }
      else
	cout << "not enough for z in dutdyvsty" << endl;
    }

  } // iteration > 1

  cout << endl
       << "DUT efficiency " << 100*effvst5.GetMean(2) << "%"
       << " from " << effvst5.GetEntries() << " in-time tracks"
       << " from " << iev << " events"
       << endl;

  // write new DUT alignment:

  cout << endl
       << "DUT alignment iteration " << DUTaligniteration + 1 << endl
       << "  alignx " << DUTalignx0 << endl
       << "  aligny " << DUTaligny0 << endl
       << "  rot    " << DUTrot << endl
       << "  tilt   " << DUTtilt << endl
       << "  turn   " << DUTturn << endl
       << "  dz     " << DUTz - zz[3] << endl
    ;

  cout << "update DUT alignment file? (y/n)" << endl;
  string ans{"y"};
  string YES{"y"};

  if( ldbmod == 0 && fabs(DUTturn) > 38 )
    cin >> ans;

  if( ans == YES ) {

    ofstream DUTalignFile( DUTalignFileName.str() );

    DUTalignFile << "# DUT alignment for run " << run
		 << " using " << iev << " events" << endl;
    ++DUTaligniteration;
    DUTalignFile << "iteration " << DUTaligniteration << endl;
    DUTalignFile << "alignx " << DUTalignx0 << endl;
    DUTalignFile << "aligny " << DUTaligny0 << endl;
    DUTalignFile << "rot " << DUTrot << endl;
    DUTalignFile << "tilt " << DUTtilt << endl;
    DUTalignFile << "turn " << DUTturn << endl;
    DUTalignFile << "dz " << DUTz - zz[3] << endl;

    DUTalignFile.close();

    cout << " to " << DUTalignFileName.str() << endl;

  }
  else
    cout << "no" << endl;

  // --------------- 
  // hotpixels for DUT
  if(create_duthotfile)
  {
     std::cout << "DUT hot pixel list for run " << run << std::endl;
     std::ofstream DUThotFile( DUThotFileName.str() );

     DUThotFile << "# DUT  hot pixel list for run " << run
	  << " with " << iev << " events"
	  << std::endl;
     DUThotFile << std::endl;
     int nmax = 0;
     int ntot = 0;
     int nhot = 0;
     for(const auto & px: pxdutmap) 
     {
        ntot += px.second;
        if( px.second > nmax )
        { 
          nmax = px.second;
        }
        if( px.second > iev/128 ) 
        {
          ++nhot;
          int ix = px.first/ny[0];
          int iy = px.first%ny[0];
	  DUThotFile << "pix "
             << setw(4) << ix
             << setw(5) << iy
             << "  " << px.second
             << std::endl;
         }
      } 
      std::cout 
         << ": active " << pxdutmap.size()
         << ", sum " << ntot 
         << ", max " << nmax
         << ", hot " << nhot
         << std::endl;

      std::cout << "hot pixel list written to " << DUThotFileName.str() << std::endl;
      DUThotFile.close();
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // done

  cout << endl << histoFile.GetName() << endl;

  cout << endl;

  return 0;
}
