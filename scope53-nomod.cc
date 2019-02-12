
// Daniel Pitzl, DESY, Jun 2018
// telescope analysis with RD53A

// make scope53
// scope53 33095
// needs runs.dat
// needs align_33095.dat from tele

#include "eudaq/FileReader.hh"
#include "eudaq/PluginManager.hh"
#include "/home/pitzl/eudaq17/main/lib/plugins/BDAQ53ConverterPlugin.cc"

#include <TFile.h>
#include <TH1I.h> // counting
#include <TH1D.h> // weighted counts
#include <TH2I.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TF1.h>

#include <sstream> // stringstream
#include <fstream> // filestream
#include <set>
#include <cmath>

using namespace std;
using namespace eudaq;

struct pixel {
  int col;
  int row;
  int tot;
};

struct cluster {
  vector <pixel> vpix; // Armin Burgmeier: list
  int size;
  int ncol, nrow;
  double col, row;
  int signal;
};

struct triplet {
  double xm;
  double ym;
  double zm;
  double sx;
  double sy;
  bool lk;
  double ttdmin;
  vector <double> vx;
  vector <double> vy;
};

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

    for( vector<pixel>::iterator p = c.vpix.begin();  p != c.vpix.end();  ++p ) {

      int nn = max( 1, p->tot ); // neighbours
      sumnn += nn;
      c.col += p->col * nn;
      c.row += p->row * nn;

      if( p->col > maxx ) maxx = p->col;
      if( p->col < minx ) minx = p->col;
      if( p->row > maxy ) maxy = p->row;
      if( p->row < miny ) miny = p->row;

    }

    //cout << "(cluster with " << c.vpix.size() << " pixels)" << endl;

    if( sumnn > 0 ) {
      c.col /= sumnn;
      c.row /= sumnn;
    }
    else {
      c.col = (*c.vpix.begin()).col;
      c.row = (*c.vpix.begin()).row;
      cout << "GetClus: cluster with signal" << sumnn << endl;
    }

    c.signal = sumnn;

    c.ncol = maxx-minx+1;
    c.nrow = maxy-miny+1;

    vc.push_back(c); // add cluster to vector

    // look for a new seed = used pixel:

    while( ( ++seed < pb.size() ) && gone[seed] );

  } // while over seeds

  delete gone;

  return vc; // vector of clusters
}

//------------------------------------------------------------------------------
vector <cluster> getClusq( vector <pixel> pb, int fCluCut = 1 ) // 1 = no gap
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

    // added all I could. determine position and append it to the list of clusters:

    c.size = c.vpix.size();
    c.col = 0;
    c.row = 0;
    double sumQ = 0;
    int minx = 999;
    int maxx = 0;
    int miny = 999;
    int maxy = 0;

    for( vector<pixel>::iterator p = c.vpix.begin();  p != c.vpix.end();  ++p ) {

      double Qpix = p->tot;

      sumQ += Qpix;

      c.col += p->col*Qpix;
      c.row += p->row*Qpix;

      if( p->col > maxx ) maxx = p->col;
      if( p->col < minx ) minx = p->col;
      if( p->row > maxy ) maxy = p->row;
      if( p->row < miny ) miny = p->row;

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

    vc.push_back(c); // add cluster to vector

    // look for a new seed = used pixel:

    while( ( ++seed < pb.size() ) && gone[seed] );

  } // while over seeds

  delete gone;

  return vc; // vector of clusters
}

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

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // further arguments:

  int lev = 999222111; // last event

  for( int i = 1; i < argc; ++i ) {

    if( !strcmp( argv[i], "-l" ) )
      lev = atoi( argv[++i] ); // last event

  } // argc

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // runs.dat:

  cout << endl;

  string geoFileName( "geo.dat" );
  double pbeam = 4.8;
  double DUTtilt0 = 0.0;
  double DUTturn0 = 0.5; // small turn will not be aligned
  double DUTtilt = DUTtilt0; // [deg]
  double DUTturn = DUTturn0; // [deg]
  int chip0 = 501;
  int fifty = 0; // default is 100x25
  int rot90 = 0; // default is straight

  ifstream runsFile( "runs.dat" );

  if( runsFile.bad() || ! runsFile.is_open() ) {
    cout << "Error opening runs.dat" << endl;
    return 1;
  }
  // can there be instructions between if and else ? no

  else {

    cout << "read runs from runs.dat" << endl;

    string hash( "#" );
    string RUN( "run" );
    string GEO( "geo" );
    string GeV( "GeV" );
    string CHIP( "chip" );
    string TURN( "turn" );
    string FIFTY( "fifty" );
    string ROT90( "rot90" );
    bool found = 0;

    while( ! runsFile.eof() ) {

      string line;
      getline( runsFile, line );

      if( line.empty() ) continue;

      stringstream tokenizer( line );
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

      if( tag == TURN ) {
	tokenizer >> DUTturn0;
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

      if( tag == FIFTY ) {
	tokenizer >> fifty;
	continue;
      }

      if( tag == ROT90 ) {
	tokenizer >> rot90;
	continue;
      }

      // anything else on the line and in the file gets ignored

    } // while getline

    if( found )
      cout 
	<< "settings for run " << run << ":" << endl
	<< "  beam " << pbeam << " GeV" << endl
	<< "  geo file " << geoFileName << endl
	<< "  nominal DUT turn " << DUTturn0 << " deg" << endl
	<< "  DUT chip " << chip0 << endl
	<< "  fifty " << fifty << endl
	<< "  rot90 " << rot90 << endl
	<< endl;
    else {
      cout << "run " << run << " not found in runs.dat" << endl;
      return 1;
    }

  } // runsFile

  runsFile.close();

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

      stringstream tokenizer( line );
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

      stringstream tokenizer( line );
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

      stringstream tokenizer( line );
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

      stringstream tokenizer( line );
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

  if( DUTaligniteration <= 1 ) {
    DUTtilt = DUTtilt0;
    DUTturn = DUTturn0;
  }

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

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // (re-)create root file:

  ostringstream rootFileName; // output string stream

  rootFileName << "scopeRD" << run << ".root";

  TFile* histoFile = new TFile( rootFileName.str(  ).c_str(  ), "RECREATE" );

  // book histos:

  double f = 4.8/pbeam;

  TH1I t1Histo = TH1I( "t1", "event time;event time [s];events / 10 ms", 100, 0, 1 );
  TH1I t2Histo = TH1I( "t2", "event time;event time [s];events / s", 300, 0, 300 );
  TH1I t3Histo = TH1I( "t3", "event time;event time [s];events / 10 s", 150, 0, 1500 );
  TH1I t4Histo = TH1I( "t4", "event time;event time [s];events /10 s", 600, 0, 6000 );
  TH1I t5Histo = TH1I( "t5", "event time;event time [s];events / 100 s", 600, 0, 60000 );
  TH1I t6Histo = TH1I( "t6", "event time;event time [h];events / 3 min", 1000, 0, 50 );

  TH1I hdtus = TH1I( "dtus", "time between events;time between events [us];events", 100, 0, 1000 );
  TH1I hdtms = TH1I( "dtms", "time between events;time between events [ms];events", 100, 0, 1000 );

  TH1I dutnframesHisto = TH1I( "dutnframes", "DUT frames;DUT frames;events", 32, -0.5, 31.5 );

  TH1I hcol[9];
  TH1I hrow[9];
  TH1I hnpx[9];
  TH2I * hmap[9];

  TH1I hncl[9];
  TH1I hsiz[9];
  TH1I hncol[9];
  TH1I hnrow[9];

  int nbx = 400; // RD53
  int nby = 192;
  int mx = nbx;
  int my = nby;

  for( int ipl = 0; ipl <= 6; ++ipl ) {

    hcol[ipl] = TH1I( Form( "col%i", ipl ),
		      Form( "%i col;col;plane %i pixels", ipl, ipl ), 
		      nbx, 0, mx );
    hrow[ipl] = TH1I( Form( "row%i", ipl ),
		      Form( "%i row;row;plane %i pixels", ipl, ipl ),
		      nby, 0, my );
    hmap[ipl] = new TH2I( Form( "map%i", ipl ),
			  Form( "%i map;col;row;plane %i pixels", ipl, ipl ),
			  nbx, 0, mx, nby, 0, my );

    nbx = nx[1]/2; // Mimosa
    nby = ny[1]/2; // for next round
    mx = nx[1];
    my = ny[1];

    hnpx[ipl] = TH1I( Form( "npx%i", ipl ),
		      Form( "%i pixel per event;pixels;plane %i events", ipl, ipl ),
		      200, 0, 200 );

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

  } // planes

  // driplets:

  TH1I hdx46 = TH1I( "dx46", "4-6 dx;4-6 dx [mm];cluster pairs", 100, -f, f );
  TH1I hdy46 = TH1I( "dy46", "4-6 dy;4-6 dy [mm];cluster pairs", 100, -f, f );

  TH1I hdridx = TH1I( "dridx", "driplet dx;driplet dx [mm];driplets", 100, -0.1*f, 0.1*f );
  TH1I hdridy = TH1I( "dridy", "driplet dy;driplet dy [mm];driplets", 100, -0.1*f, 0.1*f );

  TH1I hdridxc = TH1I( "dridxc", "driplet dx;driplet dx [mm];driplets", 100, -0.1*f, 0.1*f );
  TH1I hdridyc = TH1I( "dridyc", "driplet dy;driplet dy [mm];driplets", 100, -0.1*f, 0.1*f );

  TProfile dridxvsy =
    TProfile( "dridxvsy",
	      "driplet dx vs y;driplet yB [mm];<driplets #Deltax> [mm]",
	      110, -5.5, 5.5, -0.05*f, 0.05*f );
  TProfile dridyvsx =
    TProfile( "dridyvsx",
	      "driplet dy vs x;driplet xB [mm];<driplets #Deltay> [mm]",
	      110, -11, 11, -0.05*f, 0.05*f );

  TProfile dridxvstx =
    TProfile( "dridxvstx",
	      "driplet dx vs slope x;driplet slope x [rad];<driplets #Deltax> [mm]",
	      60, -0.003, 0.003, -0.05*f, 0.05*f );
  TProfile dridyvsty =
    TProfile( "dridyvsty",
	      "driplet dy vs slope y;driplet slope y [rad];<driplets #Deltay> [mm]",
	      60, -0.003, 0.003, -0.05*f, 0.05*f );

  TH1I drixHisto = TH1I( "drix", "driplets x;x [mm];driplets",
			  240, -12, 12 );
  TH1I driyHisto = TH1I( "driy", "driplets x;y [mm];driplets",
			  120, -6, 6 );
  TH2I * drixyHisto = new
    TH2I( "drixy", "driplets x-y;x [mm];y [mm];driplets",
	  240, -12, 12, 120, -6, 6 );
  TH1I dritxHisto = TH1I( "dritx", "driplet slope x;slope x [rad];driplets",
			    100, -0.005*f, 0.005*f );
  TH1I drityHisto = TH1I( "drity", "driplet slope y;slope y [rad];driplets",
			    100, -0.005*f, 0.005*f );

  TH1I ndriHisto = TH1I( "ndri", "driplets;driplets;events", 51, -0.5, 50.5 );

  // triplets:

  TH1I hdx13 = TH1I( "dx13", "1-3 dx;1-3 dx [mm];cluster pairs", 100, -f, f );
  TH1I hdy13 = TH1I( "dy13", "1-3 dy;1-3 dy [mm];cluster pairs", 100, -f, f );

  TH1I htridx = TH1I( "tridx", "triplet dx;triplet dx [mm];triplets", 100, -0.1, 0.1 );
  TH1I htridy = TH1I( "tridy", "triplet dy;triplet dy [mm];triplets", 100, -0.1, 0.1 );

  TH1I htridxc = TH1I( "tridxc", "triplet dx;triplet dx [mm];triplets", 100, -0.05, 0.05 );
  TH1I htridyc = TH1I( "tridyc", "triplet dy;triplet dy [mm];triplets", 100, -0.05, 0.05 );

  TProfile tridxvsx =
    TProfile( "tridxvsx",
	      "triplet dx vs x;triplet x [mm];<triplet #Deltax> [mm]",
	      120, -12, 12, -0.05, 0.05 );
  TProfile tridxvsy =
    TProfile( "tridxvsy",
	      "triplet dx vs y;triplet y [mm];<triplet #Deltax> [mm]",
	      110, -5.5, 5.5, -0.05, 0.05 );
  TProfile tridxvstx =
    TProfile( "tridxvstx",
	      "triplet dx vs slope x;triplet slope x [rad];<triplet #Deltax> [mm]",
	      60, -0.003, 0.003, -0.05, 0.05 );
  TProfile tridxvst3 =
    TProfile( "tridxvst3",
	      "triplet dx vs time;time [s];<triplet #Deltax> [mm]",
	      300, 0, 6000, -0.05, 0.05 );
  TProfile tridxvst5 =
    TProfile( "tridxvst5",
	      "triplet dx vs time;time [s];<triplet #Deltax> [mm]",
	      1000, 0, 60000, -0.05, 0.05 );

  TProfile tridyvsx =
    TProfile( "tridyvsx",
	      "triplet dy vs x;triplet x [mm];<triplet #Deltay> [mm]",
	      110, -11, 11, -0.05, 0.05 );
  TProfile tridyvsty =
    TProfile( "tridyvsty",
	      "triplet dy vs slope y;triplet slope y [rad];<triplet #Deltay> [mm]",
	      60, -0.003, 0.003, -0.05, 0.05 );
  TProfile tridyvst3 =
    TProfile( "tridyvst3",
	      "triplet dy vs time;time [s];<triplet #Deltay> [mm]",
	      300, 0, 6000, -0.05, 0.05 );
  TProfile tridyvst5 =
    TProfile( "tridyvst5",
	      "triplet dy vs time;time [h];<triplet #Deltay> [mm]",
	      1000, 0, 60000, -0.05, 0.05 );

  TH1I trixHisto = TH1I( "trix", "triplets x;x [mm];triplets",
			  240, -12, 12 );
  TH1I triyHisto = TH1I( "triy", "triplets y;y [mm];triplets",
			  120, -6, 6 );
  TH2I * trixyHisto = new
    TH2I( "trixy", "triplets x-y;x [mm];y [mm];triplets",
	  240, -12, 12, 120, -6, 6 );
  TH1I tritxHisto = TH1I( "tritx", "triplet slope x;slope x [rad];triplets",
			    100, -0.005*f, 0.005*f );
  TH1I trityHisto = TH1I( "trity", "triplet slope y;slope y [rad];triplets",
			    100, -0.005*f, 0.005*f );

  TH1I ntriHisto = TH1I( "ntri", "triplets;triplets;events", 51, -0.5, 50.5 );

  TH1I trix1Histo = TH1I( "trix1", "lone triplets x;x [mm];lone triplets",
			  240, -12, 12 );
  TH1I triy1Histo = TH1I( "triy1", "lone triplets y;y [mm];lone triplets",
			  120, -6, 6 );
  TH2I * trixy1Histo = new
    TH2I( "trixy1", "lone triplets x-y;x [mm];y [mm];lone triplets",
	  240, -12, 12, 120, -6, 6 );

  // DUT clusters:

  TH1I dutpxqHisto =
    TH1I( "dutpxq",
	  "DUT pixel signal;DUT pixel signal [ToT];pixels",
	  16, -0.5, 15.5 );

  TH1I ttdmin1Histo = TH1I( "ttdmin1",
			    "telescope triplets isolation;triplet min #Delta_{xy} [mm];triplet pairs",
			    100, 0, 1 );
  TH1I ttdmin2Histo = TH1I( "ttdmin2",
			    "telescope triplets isolation;triplet min #Delta_{xy} [mm];triplet pairs",
			    150, 0, 15 );

  // dripets - triplets:

  TH1I hsixdx = TH1I( "sixdx", "six dx;dx [mm];triplet-driplet pairs", 200, -0.2*f, 0.2*f );
  TH1I hsixdy = TH1I( "sixdy", "six dy;dy [mm];triplet-driplet pairs", 200, -0.2*f, 0.2*f );
  TH1I hsixdxc = TH1I( "sixdxc", "six dx;dx [mm];triplet-driplet pairs", 200, -0.2*f, 0.2*f );
  TH1I hsixdyc = TH1I( "sixdyc", "six dy;dy [mm];triplet-driplet pairs", 200, -0.2*f, 0.2*f );

  TProfile sixdxvsx =
    TProfile( "sixdxvsx",
	      "six #Deltax vs x;x [mm];<driplet - triplet #Deltax [mm]",
	      220, -11, 11, -0.1, 0.1 );
  TProfile sixmadxvsx =
    TProfile( "sixmadxvsx",
	      "six MAD x vs x;x [mm];driplet - triplet MAD #Deltax [mm]",
	      220, -11, 11, 0, 0.1 );
  TProfile sixmadxvsy =
    TProfile( "sixmadxvsy",
	      "six MAD x vs y;y [mm];driplet - triplet MAD #Deltax [mm]",
	      110, -5.5, 5.5, 0, 0.1 );
  TProfile sixmadxvstx =
    TProfile( "sixmadxvstx",
	      "six MAD x vs x;triplet #theta_{x} [rad];driplet - triplet MAD #Deltax [mm]",
	      80, -0.002, 0.002, 0, 0.1 );
  TProfile sixmadxvsdtx =
    TProfile( "sixmadxvsdtx",
	      "six MAD x vs x;driplet-triplet #Delta#theta_{x} [rad];driplet - triplet MAD #Deltax [mm]",
	      80, -0.002, 0.002, 0, 0.1 );
  TProfile sixdxvsy =
    TProfile( "sixdxvsy",
	      "six #Deltax vs y;y [mm];<driplet - triplet #Deltax> [mm]",
	      100, -5, 5, -0.5, 0.5 );
  TProfile sixdxvstx =
    TProfile( "sixdxvstx",
	      "six #Deltax vs slope x;slope x [rad];<driplet - triplet #Deltax> [mm]",
	      100, -0.002, 0.002, -0.5, 0.5 );
  TProfile sixdxvsdtx =
    TProfile( "sixdxvsdtx",
	      "six #Deltax vs #Delta slope x;#Delta slope x [rad];<driplet - triplet #Deltax> [mm]",
	      100, -0.002, 0.002, -0.5, 0.5 );
  TProfile sixdxvst3 =
    TProfile( "sixdxvst3",
	      "sixplet dx vs time;time [s];<sixplet #Deltax> [mm]",
	      300, 0, 6000, -0.05, 0.05 );
  TProfile sixdxvst5 =
    TProfile( "sixdxvst5",
	      "sixplet dx vs time;time [s];<sixplet #Deltax> [mm]",
	      1000, 0, 60000, -0.05, 0.05 );

  TProfile sixdyvsx =
    TProfile( "sixdyvsx",
	      "six #Deltay vs x;x [mm];<driplet - triplet #Deltay> [mm]",
	      200, -10, 10, -0.5, 0.5 );
  TProfile sixdyvsy =
    TProfile( "sixdyvsy",
	      "six #Deltay vs y;y [mm];<driplet - triplet #Deltay [mm]",
	      110, -5.5, 5.5, -0.1, 0.1 );
  TProfile sixdyvsty =
    TProfile( "sixdyvsty",
	      "six #Deltay vs slope y;slope y [rad];<driplet - triplet #Deltay> [mm]",
	      100, -0.002, 0.002, -0.5, 0.5 );
  TProfile sixdyvsdty =
    TProfile( "sixdyvsdty",
	      "six #Deltay vs #Delta slope y;#Delta slope y [rad];<driplet - triplet #Deltay> [mm]",
	      100, -0.002, 0.002, -0.5, 0.5 );
  TProfile sixdyvst3 =
    TProfile( "sixdyvst3",
	      "sixplet dy vs time;time [s];<sixplet #Deltay> [mm]",
	      300, 0, 6000, -0.05, 0.05 );
  TProfile sixdyvst5 =
    TProfile( "sixdyvst5",
	      "sixplet dy vs time;time [s];<sixplet #Deltay> [mm]",
	      1000, 0, 60000, -0.05, 0.05 );
  TProfile sixmadyvsx =
    TProfile( "sixmadyvsx",
	      "six MAD y vs x;x [mm];driplet - triplet MAD #Deltay [mm]",
	      220, -11, 11, 0, 0.1 );
  TProfile sixmadyvsy =
    TProfile( "sixmadyvsy",
	      "six MAD y vs y;y [mm];driplet - triplet MAD #Deltay [mm]",
	      110, -5.5, 5.5, 0, 0.1 );
  TProfile sixmadyvsty =
    TProfile( "sixmadyvsty",
	      "six MAD y vs #theta_{y};triplet #theta_{y} [rad];driplet - triplet MAD #Deltay [mm]",
	      80, -0.002, 0.002, 0, 0.1 );
  TProfile sixmadyvsdty =
    TProfile( "sixmadyvsdty",
	      "six MAD y vs #Delta#theta_{y};driplet-triplet #Delta#theta_{y} [rad];driplet - triplet MAD #Deltay [mm]",
	      80, -0.002, 0.002, 0, 0.1 );

  TProfile sixdtvsxav =
    TProfile( "sixdtvsxav",
	      "driplet - triplet kink_{xy} vs x;x_{avg} [mm];<sqrt(#Delta#theta_{x}^{2}+#Delta#theta_{y}^{2})> [rad]",
	      240, -12, 12, 0, 0.1 );
  TProfile sixdtvsyav =
    TProfile( "sixdtvsyav",
	      "driplet - triplet kink_{xy} vs y;y_{avg} [mm];<sqrt(#Delta#theta_{x}^{2}+#Delta#theta_{y}^{2})> [rad]",
	      120, -6, 6, 0, 0.1 );
  TProfile2D * sixdtvsxyav = new
    TProfile2D( "sixdtvsxyav",
		"driplet - triplet kink_{xy} vs x-y;x_{avg} [mm];y_{avg} [mm];<sqrt(#Delta#theta_{x}^{2}+#Delta#theta_{y}^{2})> [rad]",
		120, -12, 12, 60, -6, 6, 0, 0.1 );

  TH1I hsixdtx =
    TH1I( "sixdtx",
	  "driplet slope x - triplet slope x;driplet slope x - triplet slope x;driplet-triplet pairs",
	  100, -0.005*f, 0.005*f );     
  TH1I hsixdty =
    TH1I( "sixdty",
	  "driplet slope y - triplet slope y;driplet slope y - triplet slope y;driplet-triplet pairs",
	  100, -0.005*f, 0.005*f );     

  TProfile sixdtvsx =
    TProfile( "sixdtvsx",
	      "driplet - triplet kink_{xy} vs x;x_{mid} [mm];<sqrt(#Delta#theta_{x}^{2}+#Delta#theta_{y}^{2})> [rad]",
	      110, -11, 11, 0, 0.1 );
  TProfile2D * sixdtvsxy = new
    TProfile2D( "sixdtvsxy",
		"driplet - triplet kink_{xy} vs x-y;x_{DUT} [mm];y_{DUT} [mm];<sqrt(#Delta#theta_{x}^{2}+#Delta#theta_{y}^{2})> [rad]",
		110, -11, 11, 55, -5.5, 5.5, 0, 0.1 );

  TProfile sixdtvsxm =
    TProfile( "sixdtvsxm",
	      "driplet - triplet kink_{xy} vs xmod;track x mod 0.1 [mm];<sqrt(#Delta#theta_{x}^{2}+#Delta#theta_{y}^{2})> [rad]",
	      50, 0, 0.1, 0, 0.1 );
  TProfile sixdtvsym =
    TProfile( "sixdtvsym",
	      "driplet - triplet kink_{xy} vs ymod;track y mod 0.1 [mm];<sqrt(#Delta#theta_{x}^{2}+#Delta#theta_{y}^{2})> [rad]",
	      50, 0, 0.1, 0, 0.1 );
  TProfile2D * sixdtvsxmym = new // bump bonds ?
    TProfile2D( "sixdtvsxmym",
	      "driplet - triplet kink_{xy} vs xmod ymod;track x mod 100 [#mum];track y mod 100 [#mum];<sqrt(#Delta#theta_{x}^{2}+#Delta#theta_{y}^{2})> [rad]",
		50, 0, 100, 50, 0, 100, 0, 0.1 );

  TH2I * sixxyHisto = new
    TH2I( "sixxy", "sixplets at z DUT;x [mm];y [mm];sixplets",
	  240, -12, 12, 120, -6, 6 );

  // DUT vs triplets:

  TH2I * dutxyHisto = new
    TH2I( "dutxy", "tracks at DUT;x [mm];y [mm];tracks",
	  240, -12, 12, 120, -6, 6 );
  TH2I * dutxy1Histo = new
    TH2I( "dutxy1", "lone tracks at DUT;x [mm];y [mm];lone tracks",
	  240, -12, 12, 120, -6, 6 );
  TH1I dutx1Histo = TH1I( "dutx1",
			  "lone track at DUT x;track x at DUT [mm];lone tracks",
			  240, -12, 12 );
  TH1I duty1Histo = TH1I( "duty1",
			  "lone track at DUT y;track y at DUT [mm];lone tracks",
			  120, -6, 6 );

  TH2I * dutxxHisto = new
    TH2I( "dutxx", "tracks vs DUT in x;track x [mm];DUT cluster x [mm];track-cluster pairs",
	  500, -5, 5, 500, -5, 5 );
  TH2I * dutyyHisto = new
    TH2I( "dutyy", "tracks vs DUT in y;track y [mm];DUT cluster y [mm];track-cluster pairs",
	  500, -5, 5, 500, -5, 5 );

  TH1I dutdxaHisto = TH1I( "dutdxa",
			   "DUT - track dx;DUT cluster - track #Deltax [mm];DUT clusters",
			   400, -8, 8 );
  TH1I dutsxaHisto = TH1I( "dutsxa",
			   "DUT + track sx;DUT cluster + track #Deltax [mm];DUT clusters",
			   400, -8, 8 );
  TH1I dutdyaHisto = TH1I( "dutdya",
			   "DUT - track dy;DUT cluster - track #Deltay [mm];DUT clusters",
			   500, -5, 5 );
  TH1I dutsyaHisto = TH1I( "dutsya",
			   "DUT + track sy;DUT cluster + track #Deltay [mm];DUT clusters",
			   500, -5, 5 );

  TH1I dutdxHisto = TH1I( "dutdx",
			   "DUT - track dx;DUT cluster - track #Deltax [mm];DUT clusters",
			   200, -0.2, 0.2 );
  TH1I dutdx0Histo = TH1I( "dutdx0",
			   "DUT - track dx;DUT cluster - track #Deltax [mm];DUT clusters",
			   500, -5, 5 ); // binning like dutxx
  TH1I dutdxcHisto = TH1I( "dutdxc",
			   "DUT - track dx;DUT cluster - track #Deltax [mm];DUT clusters",
			   200, -0.2, 0.2 );
  TH1I dutdxc1Histo = TH1I( "dutdxc1",
			   "DUT - track dx;DUT cluster - track #Deltax [mm];DUT 1-px clusters",
			   200, -0.2, 0.2 );
  TH1I dutdxc2Histo = TH1I( "dutdxc2",
			   "DUT - track dx;DUT cluster - track #Deltax [mm];DUT 2-px clusters",
			   200, -0.2, 0.2 );
  TH1I dutdxcqHisto = TH1I( "dutdxcq",
			    "DUT - track dx Landau peak;DUT cluster - track #Deltax [mm];Landau peak DUT clusters",
			    500, -0.2, 0.2 );
  TProfile dutdxvsx =
    TProfile( "dutdxvsx",
	      "DUT #Deltax vs x;x track [mm];<cluster - track #Deltax> [mm]",
	      200, -10, 10, -0.1, 0.1 );
  TProfile dutdxvsy =
    TProfile( "dutdxvsy",
	      "DUT #Deltax vs y;y track [mm];<cluster - track #Deltax> [mm]",
	      160, -8, 8, -0.1, 0.1 );
  TProfile dutdxvstx =
    TProfile( "dutdxvstx",
	      "DUT #Deltax vs #theta_{x};x track slope [rad];<cluster - track #Deltax> [mm]",
	      80, -0.002, 0.002, -0.1, 0.1 );
  TProfile dutdxvsxm =
    TProfile( "dutdxvsxm",
	      "DUT #Deltax vs xmod;x track mod 0.1 [mm];<cluster - track #Deltax> [mm]",
	     50, 0, 0.1, -0.1, 0.1 );
  TProfile dutdxvst5 =
    TProfile( "dutdxvst5",
	      "DUT #Deltax vs time;time [s];<DUT #Deltax> [mm]",
	      1000, 0, 60000, -0.1, 0.1 );
  TProfile dutdyvst5 =
    TProfile( "dutdyvst5",
	      "DUT #Deltay vs time;time [s];<DUT #Deltay> [mm]",
	      1000, 0, 60000, -0.1, 0.1 );

  TProfile dutmadxvsx =
    TProfile( "dutmadxvsx",
	      "DUT MAD(#Deltax) vs x;x track [mm];MAD(#Deltax) [mm]",
	      200, -10, 10, 0, 0.1 );
  TProfile dutmadxvsxm =
    TProfile( "dutmadxvsxm",
	      "DUT MAD(#Deltax) vs xmod;x track mod 0.1 [mm];MAD(#Deltax) [mm]",
	      50, 0, 0.1, 0, 0.1 );
  TProfile dutmadxvstx =
    TProfile( "dutmadxvstx",
	      "DUT MAD(#Deltax) vs #theta_{x};x track slope [rad];MAD(#Deltax) [mm]",
	      80, -0.002, 0.002, 0, 0.1 );
  TProfile dutmadxvsq =
    TProfile( "dutmadxvsq",
	      "DUT MAD(#Deltax) vs Q;cluster signal [ToT];MAD(#Deltax) [mm]",
	      80, 0, 80, 0, 0.1 );

  TH1I dutdyHisto =
    TH1I( "dutdy",
	  "DUT - track dy;DUT cluster - track #Deltay [mm];DUT clusters",
	  500, -0.5, 0.5 );
  TH1I dutdy0Histo =
    TH1I( "dutdy0",
	  "DUT - track dy;DUT cluster - track #Deltay [mm];DUT clusters",
	  500, -5, 5 ); // binning like dutyy
  TH1I dutdycHisto =
    TH1I( "dutdyc",
	  "DUT - track dy;DUT cluster - track #Deltay [mm];DUT clusters",
	  500, -0.5, 0.5 );
  TH1I dutdyqHisto =
    TH1I( "dutdyq",
	  "DUT - track dy Landau peak;DUT cluster - track #Deltay [mm];Landau peak DUT clusters",
	  500, -0.5, 0.5 );
  TH1I dutdyqtHisto =
    TH1I( "dutdyqt",
	  "DUT - track dy Landau peak;DUT cluster - track #Deltay [mm];Landau peak DUT clusters",
	  500, -0.5, 0.5 );
  TH1I dutdyq3tHisto =
    TH1I( "dutdyq3t",
	  "DUT - track dy Landau peak;DUT cluster - track #Deltay [mm];Landau peak DUT clusters",
	  500, -0.5, 0.5 );
  TProfile dutdyvsx =
    TProfile( "dutdyvsx",
	      "DUT #Deltay vs x;x track [mm];<cluster - track #Deltay> [mm]",
	      200, -10, 10, -0.1, 0.1 );
  TProfile dutdyvsy =
    TProfile( "dutdyvsy",
	      "DUT #Deltay vs y;y track [mm];<cluster - track #Deltay> [mm]",
	      160, -8, 8, -0.1, 0.1 );
  TProfile dutdyvsty =
    TProfile( "dutdyvsty",
	      "DUT #Deltay vs #theta_{y};y track slope [rad];<cluster - track #Deltay> [mm]",
	      80, -0.002, 0.002, -0.1, 0.1 );
  TProfile dutdyvsxm =
    TProfile( "dutdyvsxm",
	      "DUT #Deltay vs xmod;x track mod 0.1 [mm];<cluster - track #Deltay> [mm]",
	     50, 0, 0.1, -0.1, 0.1 );
  TProfile dutdyvsym =
    TProfile( "dutdyvsym",
	      "DUT #Deltay vs ymod;y track mod 0.1 [mm];<cluster - track #Deltay> [mm]",
	      50, 0, 0.1, -0.1, 0.1 );
  TProfile2D * dutdyvsxmym = new
    TProfile2D( "dutdyvsxmym",
		"DUT #Deltay vs xmod ymod;x track mod 0.1 [mm];y track mod 0.1 [mm];<#Deltay> [mm]",
		50, 0, 0.1, 50, 0, 0.1, -0.1, 0.1 );

  TProfile dutmadyvsx =
    TProfile( "dutmadyvsx",
	      "DUT MAD(#Deltay) vs x;x track [mm];MAD(#Deltay) [mm]",
	      200, -10, 10, 0, 0.1 );
  TProfile dutmadyvsy =
    TProfile( "dutmadyvsy",
	      "DUT MAD(#Deltay) vs y;y track [mm];MAD(#Deltay) [mm]",
	      160, -8, 8, 0, 0.1 );
  TProfile dutmadyvstx =
    TProfile( "dutmadyvstx",
	      "DUT MAD(#Deltay) vs #theta_{x};x track slope [rad];MAD(#Deltay) [mm]",
	      80, -0.002, 0.002, 0, 0.1 );
  TProfile dutmadyvsty =
    TProfile( "dutmadyvsty",
	      "DUT MAD(#Deltay) vs #theta_{y};y track slope [rad];MAD(#Deltay) [mm]",
	      80, -0.002, 0.002, 0, 0.1 );
  TProfile dutmadyvsxm =
    TProfile( "dutmadyvsxm",
	      "DUT MAD(#Deltay) vs xmod;x track mod 0.1 [mm];MAD(#Deltay) [mm]",
	      50, 0, 0.1, 0, 0.1 );
  TProfile dutmadyvsym =
    TProfile( "dutmadyvsym",
	      "DUT MAD(#Deltay) vs ymod;y track mod 0.1 [mm];MAD(#Deltay) [mm]",
	      50, 0, 0.1, 0, 0.1 );
  TProfile dutmadyvsq =
    TProfile( "dutmadyvsq",
	      "DUT MAD(#Deltay) vs Q;cluster signal [ToT];MAD(#Deltay) [mm]",
	      80, 0, 80, 0, 0.1 );
  TProfile dutmadyvst =
    TProfile( "dutmadyvst",
	      "DUT MAD(#Deltay) vs time;time [s];MAD(#Deltay) [mm]",
	      100, 0, 1000, 0, 0.1 );

  TH1I dutqHisto = TH1I( "dutq",
			 "DUT linked clusters;DUT cluster signal [ToT];linked DUT cluster",
			 80, 0, 80 );
  TH1I dutq0Histo = TH1I( "dutq0",
			 "DUT linked clusters;DUT normal cluster signal [ToT];linked DUT cluster",
			 80, 0, 80 );

  TH1I dutncolHisto =
    TH1I( "dutncol",
	  "DUT linked cluster cols;DUT cluster size [columns];linked DUT cluster",
	  20, 0.5, 20.5 );
  TH1I dutnrowHisto =
    TH1I( "dutnrow",
	  "DUT linked cluster rows;DUT cluster size [rows];linked DUT cluster",
	  20, 0.5, 20.5 );
  TH1I dutnpxHisto =
    TH1I( "dutnpx",
	  "DUT linked cluster size;DUT cluster size [pixels];linked DUT cluster",
	  20, 0.5, 20.5 );

  TProfile dutncolvsxm =
    TProfile( "dutncolvsxm",
	      "DUT cluster size vs xmod;x track mod 0.1 [mm];DUT <cluster size> [columns]",
	      50, 0, 0.1, 0, 20 );
  TProfile dutncolvsym =
    TProfile( "dutncolvsym",
	      "DUT cluster size vs ymod;y track mod 0.1 [mm];DUT <cluster size> [columns]",
	      50, 0, 0.1, 0, 20 );
  TProfile dutnrowvsxm =
    TProfile( "dutnrowvsxm",
	      "DUT cluster size vs xmod;x track mod 0.1 [mm];DUT <cluster size> [rows]",
	      50, 0, 0.1, 0, 20 );
  TProfile dutnrowvsym =
    TProfile( "dutnrowvsym",
	      "DUT cluster size vs ymod;y track mod 0.1 [mm];DUT <cluster size> [rows]",
	      50, 0, 0.1, 0, 20 );
  TProfile2D * dutnpxvsxmym = new
    TProfile2D( "dutnpxvsxmym",
		"DUT cluster size vs xmod ymod;x track mod 0.1 [mm];y track mod 0.1 [mm];DUT <cluster size> [pixels]",
		50, 0, 0.1, 50, 0, 0.1, 0, 20 );

  TH2I * trixylkHisto = new
    TH2I( "trixylk", "linked triplets x-y;x [mm];y [mm];linked triplets",
	  240, -12, 12, 120, -6, 6 );
  TH1I trixlkHisto = TH1I( "trixlk",
			   "linked track at DUT x;track x at DUT [mm];linked tracks",
			   240, -12, 12 );
  TH1I triylkHisto = TH1I( "triylk",
			   "linked track at DUT y;track y at DUT [mm];linked tracks",
			   120, -6, 6 );
  TH1I trix1lkHisto = TH1I( "trix1lk",
			   "linked lone track at DUT x;track x at DUT [mm];linked lone tracks",
			    240, -12, 12 );
  TH1I triy1lkHisto = TH1I( "triy1lk",
			    "linked lone track at DUT y;track y at DUT [mm];linked lone tracks",
			    120, -6, 6 );

  TH1I dutx1lkHisto = TH1I( "dutx1lk",
			    "linked lone track at DUT x;track x at DUT [mm];linked lone tracks",
			    240, -12, 12 );
  TH1I duty1lkHisto = TH1I( "duty1lk",
			    "linked lone track at DUT y;track y at DUT [mm];linked lone tracks",
			    120, -6, 6 );

  TH2I * dutxylkHisto = new
    TH2I( "dutxylk", "linked tracks at DUT;x [mm];y [mm];linked tracks",
	  480, -12, 12, 240, -6, 6 );
  TH2I * dutxy1lkHisto = new
    TH2I( "dutxy1lk", "linked lone tracks at DUT;x [mm];y [mm];linked lone tracks",
	  480, -12, 12, 240, -6, 6 );

  TH2I * dutlkxmymHisto = new
    TH2I( "dutlkxmym",
	  "linked tracks at DUT xmod ymod;x track mod 0.1 [mm];y track mod 0.1 [mm];linked tracks",
	  50, 0, 0.1, 50, 0, 0.1 );

  TH1I dutlkcolHisto = TH1I( "dutlkcol",
			     "DUT linked col;DUT linked col;linked DUT cluster",
			     216, 0, 432 );
  TH1I dutlkrowHisto = TH1I( "dutlkrow",
			     "DUT linked row;DUT linked row;linked DUT cluster",
			     192, 0, 192 );
  TH1I dutpxqlkHisto =
    TH1I( "dutpxqlk",
	  "DUT pixel signal;DUT pixel signal [ToT];linked pixels",
	  16, -0.5, 15.5 );

  TProfile dutlkvst1 =
    TProfile( "dutlkvst1",
	      "track-DUT links vs time;time [s];tracks with DUT links",
	      300, 0, 300, -0.5, 1.5 );
  TProfile dutlkvst2 =
    TProfile( "dutlkvst2",
	      "track-DUT links vs time;time [s];tracks with DUT links",
	      150, 0, 1500, -0.5, 1.5 );
  TProfile dutlkvst3 =
    TProfile( "dutlkvst3",
	      "track-DUT links vs time;time [s];tracks with DUT links",
	      300, 0, 30000, -0.5, 1.5 );
  TProfile dutlkvst5 =
    TProfile( "dutlkvst5",
	      "track-DUT links vs time;time [s];tracks with DUT links",
	      1000, 0, 60000, -0.5, 1.5 );

  TH1I ntrilkHisto = TH1I( "ntrilk", "track - DUT links;track - DUT links;events",
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

  int iev = 0;
  uint64_t evTLU0 = 0;
  const double fTLU = 384E6; // 384 MHz TLU clock
  uint64_t prevTLU = 0;

  do {
    // Get next event:
    DetectorEvent evt = reader->GetDetectorEvent();

    if( evt.IsBORE() )
      eudaq::PluginManager::Initialize(evt);

    bool ldbg = 0;

    if( iev <  0 )
      ldbg = 1;

    if( lev < 100 )
      ldbg = 1;

    uint64_t evTLU = evt.GetTimestamp(); // 384 MHz = 2.6 ns
    if( iev < 2  )
      evTLU0 = evTLU;
    double evsec = (evTLU - evTLU0) / fTLU;
    t1Histo.Fill( evsec );
    t2Histo.Fill( evsec );
    t3Histo.Fill( evsec );
    t4Histo.Fill( evsec );
    t5Histo.Fill( evsec );
    t6Histo.Fill( evsec/3600 );

    double evdt = (evTLU - prevTLU) / fTLU;
    hdtus.Fill( evdt * 1E6 ); // [us]
    hdtms.Fill( evdt * 1E3 ); // [ms]
    prevTLU = evTLU;

    if( iev < 10 || ldbg )
      cout << "scope processing  " << run << "." << iev << "  taken " << evsec << endl;
    else if( iev < 100 && iev%10 == 0 )
      cout << "scope processing  " << run << "." << iev << "  taken " << evsec << endl;
    else if( iev < 1000 && iev%100 == 0 )
      cout << "scope processing  " << run << "." << iev << "  taken " << evsec << endl;
    else if( iev%1000 == 0 )
      cout << "scope processing  " << run << "." << iev << "  taken " << evsec << endl;

    StandardEvent sevt = eudaq::PluginManager::ConvertToStandard(evt);

    if( ldbg ) cout << "planes " << sevt.NumPlanes() << endl;

    vector < cluster > cl[9];

    for( size_t iplane = 0; iplane < sevt.NumPlanes(); ++iplane ) {

      const eudaq::StandardPlane &plane = sevt.GetPlane(iplane);

      std::vector<double> pxl = plane.GetPixels<double>();

      if( ldbg )
	std::cout
	  << "  " << iplane
	  << ": plane " << plane.ID()
	  << " " << plane.Type()
	  << " " << plane.Sensor()
	  << " frames " << plane.NumFrames() // always 32 for BDAQ53
	  << " pivot " << plane.PivotPixel()
	  << " total " << plane.TotalPixels()
	  << " hits " << plane.HitPixels()
	  ;

      int ipl = plane.ID(); // 0 = DUT, 1..6 = Mimosa

      if( ipl < 0 || ipl > 6 ) {
	cout << "event " << iev << " wrong plane number " << ipl << endl;
	continue;
      }

      if( ipl == iDUT )
	dutnframesHisto.Fill( plane.NumFrames() );

      vector <pixel> pb; // for clustering

      for( size_t ipix = 0; ipix < pxl.size(); ++ipix ) {

	if( ldbg ) 
	  std::cout << ": " << plane.GetX(ipix)
		    << "." << plane.GetY(ipix)
		    << "." << plane.GetPixel(ipix) << " ";

	int ix = plane.GetX(ipix); // column
	int iy = plane.GetY(ipix); // row
	int tot = plane.GetPixel(ipix); // ToT 0..15

	// skip hot pixels:

	int ipx = ix*ny[ipl] + iy;
	if( hotset[ipl].count(ipx) ) continue;

	pixel px;
	px.col = ix; // col
	px.row = iy; // row
	px.tot = tot;

	if( ipl == iDUT ) {

	  px.tot += 1; // shift from zero

	  if( fifty ) {
	    px.col = ix; // 50x50
	    px.row = iy;
	  }
	  else {
	    px.col = ix/2; // 100 um
	    if( ix%2 ) 
	      px.row = 2*iy + 1;
	    else
	      px.row = 2*iy + 0;
	  }

	} // DUT

	pb.push_back(px);

	hcol[ipl].Fill( ix );
	hrow[ipl].Fill( iy );
	hmap[ipl]->Fill( ix, iy );

      } // pix

      hnpx[ipl].Fill( pb.size() );

      if( ldbg ) std::cout << std::endl;

      // clustering:

      if( ipl == iDUT )
	cl[ipl] = getClusq( pb );
      else
	cl[ipl] = getClusn( pb );

      if( ldbg ) cout << "    clusters " << cl[ipl].size() << endl;

      hncl[ipl].Fill( cl[ipl].size() );

      for( vector<cluster>::iterator c = cl[ipl].begin(); c != cl[ipl].end(); ++c ) {

	hsiz[ipl].Fill( c->size );
	hncol[ipl].Fill( c->ncol );
	hnrow[ipl].Fill( c->nrow );

      } // cl

    } // planes

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
	  }

	  if( fabs( dxm ) < 0.05 ) {
	    hdridyc.Fill( dym );
	    dridyvsx.Fill( xB, dym );
	    dridyvsty.Fill( slpy, dym );
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

	  vector <double> ux(3);
	  ux[0] = xA;
	  ux[1] = xB;
	  ux[2] = xC;
	  dri.vx = ux;

	  vector <double> uy(3);
	  uy[0] = yA;
	  uy[1] = yB;
	  uy[2] = yC;
	  dri.vy = uy;

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

	  vector <double> ux(3);
	  ux[0] = xA;
	  ux[1] = xB;
	  ux[2] = xC;
	  tri.vx = ux;

	  vector <double> uy(3);
	  uy[0] = yA;
	  uy[1] = yB;
	  uy[2] = yC;
	  tri.vy = uy;

	  triplets.push_back(tri);

	  trixHisto.Fill( avx );
	  triyHisto.Fill( avy );
	  trixyHisto->Fill( avx, avy );
	  tritxHisto.Fill( slpx );
	  trityHisto.Fill( slpy );

	} // cl B

      } // cl C

    } // cl A

    ntriHisto.Fill( triplets.size() );

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // DUT clusters:

    for( vector<cluster>::iterator c = cl[iDUT].begin(); c != cl[iDUT].end(); ++c )

      for( vector<pixel>::iterator px = c->vpix.begin(); px != c->vpix.end(); ++px )
	dutpxqHisto.Fill( px->tot );

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // triplets vs DUT:

    int nm = 0;
    int ntrilk = 0;

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

    for( unsigned int iA = 0; iA < triplets.size(); ++iA ) { // iA = upstream

      double xmA = triplets[iA].xm;
      double ymA = triplets[iA].ym;
      double zmA = triplets[iA].zm;
      double sxA = triplets[iA].sx;
      double syA = triplets[iA].sy;

      if( triplets.size() == 1 ) {
	trix1Histo.Fill( xmA );
	triy1Histo.Fill( ymA );
	trixy1Histo->Fill( xmA, ymA );
      }

      double zA = DUTz - zmA; // z DUT from mid of triplet
      double xA = xmA + sxA * zA; // triplet impact point on DUT
      double yA = ymA + syA * zA;

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // intersect inclined track with tilted DUT plane:

      double zc = (Nz*zA - Ny*ymA - Nx*xmA) / (Nx*sxA + Ny*syA + Nz); // from zmA
      double yc = ymA + syA * zc;
      double xc = xmA + sxA * zc;

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

      double x4 = x3 + DUTalignx; // shift to zero
      double y4 = y3 + DUTaligny; // invert y, shift to zero

      double xmod = fmod( 9.000 + x4, 0.1 ); // [0,0.1] mm
      double ymod = fmod( 9.000 + y4, 0.1 ); // [0,0.1] mm
      if( fifty)
	ymod = fmod( 9.050 + y4, 0.1 );

      bool fiducial = 1;

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

	// driplet at DUT:

	double zB = zc + zmA - zmB; // z from mid of driplet to DUT intersect
	double xB = xmB + sxB * zB; // driplet at DUT
	double yB = ymB + syB * zB;

	// driplet - triplet:

	double dx = xB - xc; // at DUT intersect
	double dy = yB - yc;
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

	  // average driplet and triplet at DUT:

	  double xa = 0.5 * ( xB + xc );
	  double ya = 0.5 * ( yB + yc );

	  sixxyHisto->Fill( xa, ya );

	  // compare slopes:

	  hsixdtx.Fill( dtx );
	  hsixdty.Fill( dty ); // width: 0.3 mrad
	  sixdtvsxav.Fill( xa, dtxy );
	  sixdtvsyav.Fill( ya, dtxy );
	  sixdtvsxyav->Fill( xa, ya, dtxy );
	  sixdtvsx.Fill( x4, dtxy );
	  sixdtvsxy->Fill( x4, y4, dtxy );
	  if( fiducial ) {
	    sixdtvsxm.Fill( xmod, dtxy );
	    sixdtvsym.Fill( ymod, dtxy );
	    sixdtvsxmym->Fill( xmod*1E3, ymod*1E3, dtxy ); // gStyle->SetPalette(55) 105 107
	  }

	  // transform into DUT system: (passive).

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

	  // update:

	  x4 = x8;
	  y4 = y8;
	  xmod = fmod( 9.000 + x8, 0.1 ); // [0,0.1] mm
	  ymod = fmod( 9.000 + y8, 0.1 ); // [0,0.1] mm
	  if( fifty)
	    ymod = fmod( 9.050 + y4, 0.1 );

	} // six match

      } // driplets

      bool fidx = 1;
      if( x4 >  4.7 ) fidx = 0; // rot90
      if( x4 < -4.7 ) fidx = 0;

      bool fidy = 1;
      if( y4 > 3.5 ) fidy = 0;
      if( y4 <-2.2 ) fidy = 0; // 501 with packman cutout

      dutxyHisto->Fill( x4, y4 );
      if( triplets.size() == 1 ) {
	dutxy1Histo->Fill( x4, y4 );
	if( fidy )
	  dutx1Histo.Fill( x4 );
	if( fidx )
	  duty1Histo.Fill( y4 );
      }

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // triplets vs DUT clusters:

      for( vector<cluster>::iterator c = cl[iDUT].begin(); c != cl[iDUT].end(); ++c ) {

	double ccol = c->col;
	double crow = c->row;

	double dutx = ( ccol + 0.5 - nx[iDUT]/2 ) * ptchx[iDUT]; // -3.9..3.9 mm
	double duty = ( crow + 0.5 - ny[iDUT]/2 ) * ptchy[iDUT]; // -4..4 mm

	if( rot90 ) {
	  dutx = ( crow + 0.5 - ny[iDUT]/2 ) * ptchy[iDUT]; // -4..4 mm
	  duty = ( ccol + 0.5 - nx[iDUT]/2 ) * ptchx[iDUT]; // -3.9..3.9 mm
	}

	double Q = c->signal;
	double Q0 = Q*norm;

	int npx = c->size;

	// residuals for pre-alignment:

	dutdxaHisto.Fill( dutx - x3 );
	dutsxaHisto.Fill(-dutx - x3 );
 	dutdyaHisto.Fill( duty - y3 );
 	dutsyaHisto.Fill(-duty - y3 );

	double dutdx = dutx - x4;
	double dutdy = duty - y4;
	if( rot90 ) dutdy = -duty - y4;

	dutdxHisto.Fill( dutdx );
	dutdyHisto.Fill( dutdy );

	if( fabs( dutdy ) < ycutDUT ) {

	  dutxxHisto->Fill( x4, dutx );
	  dutdxcHisto.Fill( dutdx );
	  if(      c->nrow == 1 )
	    dutdxc1Histo.Fill( dutdx );
	  else if( c->nrow == 2 )
	    dutdxc2Histo.Fill( dutdx );

	  dutdxvsx.Fill( x4, dutdx ); // for turn
	  dutdxvsy.Fill( y4, dutdx ); // for rot
	  dutdxvstx.Fill( sxA, dutdx );
	  dutdxvsxm.Fill( xmod, dutdx );
	  dutdxvst5.Fill( evsec, dutdx );

	  dutmadxvsx.Fill( x4, fabs(dutdx) );
	  dutmadxvsxm.Fill( xmod, fabs(dutdx) );
	  dutmadxvstx.Fill( sxA, fabs(dutdx) );
	  dutmadxvsq.Fill( Q0, fabs(dutdx) );
	  if( Q0 > 10 && Q0 < 20 ) // 33483
	    dutdxcqHisto.Fill( dutdx );

	} // cut y

	if( fabs( dutdx ) < xcutDUT ) {

	  dutyyHisto->Fill( y4, duty );
	  dutdy0Histo.Fill( dutdy );
	  dutdycHisto.Fill( dutdy );

	  dutmadyvsq.Fill( Q0, fabs(dutdy) );

	  dutdyvsx.Fill( x4, dutdy ); // for rot
	  dutdyvsy.Fill( y4, dutdy ); // for tilt
	  dutdyvsty.Fill( syA, dutdy );
	  dutdyvsxm.Fill( xmod, dutdy );
	  dutdyvsym.Fill( ymod, dutdy );
	  dutdyvsxmym->Fill( xmod, ymod, dutdy );
	  dutdyvst5.Fill( evsec, dutdy );

	  dutmadyvsx.Fill( x4, fabs(dutdy) );
	  dutmadyvsy.Fill( y4, fabs(dutdy) );
	  dutmadyvstx.Fill( sxA, fabs(dutdy) );
	  dutmadyvsty.Fill( syA, fabs(dutdy) );
	  dutmadyvsxm.Fill( xmod, fabs(dutdy) );
	  dutmadyvsym.Fill( ymod, fabs(dutdy) );
	  dutmadyvst.Fill( evsec, fabs(dutdy) );

	} // cut x

	// cut x and y:

	if( fabs( dutdx ) < xcutDUT &&
	    fabs( dutdy ) < ycutDUT ) {

	  dutqHisto.Fill( Q );

	  dutq0Histo.Fill( Q0 );

	  dutncolHisto.Fill( c->ncol );
	  dutnrowHisto.Fill( c->nrow );
	  dutnpxHisto.Fill( npx );
	  dutncolvsxm.Fill( xmod, c->ncol );
	  dutncolvsym.Fill( ymod, c->ncol );
	  dutnrowvsxm.Fill( xmod, c->nrow );
	  dutnrowvsym.Fill( ymod, c->nrow );
	  dutnpxvsxmym->Fill( xmod, ymod, npx );

	  trixylkHisto->Fill( xA, yA );
	  trixlkHisto.Fill( xA );
	  triylkHisto.Fill( yA );
	  if( triplets.size() == 1 ) {
	    trix1lkHisto.Fill( xA );
	    triy1lkHisto.Fill( yA );
	  }
	  dutlkxmymHisto->Fill( xmod, ymod );
	  dutlkcolHisto.Fill( ccol );
	  dutlkrowHisto.Fill( crow );

	  for( vector<pixel>::iterator px = c->vpix.begin(); px != c->vpix.end(); ++px )
	    dutpxqlkHisto.Fill( px->tot );

	  triplets[iA].lk = 1;
	  nm = 1; // we have a DUT-triplet match in this event
	  ++ntrilk;

	} // DUT link x and y

	if( fabs( dutdx ) < 0.5 && // for eff
	    fabs( dutdy ) < 0.5 ) {

	  dutxylkHisto->Fill( x4, y4 );

	  if( triplets.size() == 1 ) {
	    dutxy1lkHisto->Fill( x4, y4 );
	    if( fidy )
	      dutx1lkHisto.Fill( x4 );
	    if( fidx )
	      duty1lkHisto.Fill( y4 );
	  }

	}

      } // loop DUT clusters

    } // loop triplets

    dutlkvst1.Fill( evsec, nm ); // DUT yield vs time
    dutlkvst2.Fill( evsec, nm );
    dutlkvst3.Fill( evsec, nm );
    dutlkvst5.Fill( evsec, nm );
    ntrilkHisto.Fill( ntrilk );

    ++iev;

  } while( reader->NextEvent() && iev < lev );

  delete reader;

  cout << "done after " << iev << " events" << endl;
  histoFile->Write();
  histoFile->Close();

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // DUT alignment:

  if( dutdxaHisto.GetEntries() > 999 ) {

    double newDUTalignx = DUTalignx;
    double newDUTaligny = DUTaligny;

    if( dutdxaHisto.GetMaximum() > dutsxaHisto.GetMaximum() ) {

      cout << endl << dutdxaHisto.GetTitle()
	   << " bin " << dutdxaHisto.GetBinWidth(1)
	   << endl;
      TF1 * fgp0 = new TF1( "fgp0", "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]", -1, 1 );
      double xpk = dutdxaHisto.GetBinCenter( dutdxaHisto.GetMaximumBin() );
      fgp0->SetParameter( 0, dutdxaHisto.GetMaximum() ); // amplitude
      fgp0->SetParameter( 1, xpk );
      fgp0->SetParameter( 2, dutdxaHisto.GetBinWidth(1) ); // sigma
      fgp0->SetParameter( 3, dutdxaHisto.GetBinContent( dutdxaHisto.FindBin(xpk-1) ) ); // BG
      dutdxaHisto.Fit( "fgp0", "q", "", xpk-1, xpk+1 );
      cout << "Fit Gauss + BG:"
	   << endl << "  A " << fgp0->GetParameter(0)
	   << endl << "mid " << fgp0->GetParameter(1)
	   << endl << "sig " << fgp0->GetParameter(2)
	   << endl << " BG " << fgp0->GetParameter(3)
	   << endl;
      newDUTalignx = fgp0->GetParameter(1);
      //delete fgp0;
    }
    else {
      cout << endl << dutsxaHisto.GetTitle()
	   << " bin " << dutsxaHisto.GetBinWidth(1)
	   << endl;
      TF1 * fgp0 = new TF1( "fgp0", "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]", -1, 1 );
      double xpk = dutsxaHisto.GetBinCenter( dutsxaHisto.GetMaximumBin() );
      fgp0->SetParameter( 0, dutsxaHisto.GetMaximum() ); // amplitude
      fgp0->SetParameter( 1, xpk );
      fgp0->SetParameter( 2, dutsxaHisto.GetBinWidth(1) ); // sigma
      fgp0->SetParameter( 3, dutsxaHisto.GetBinContent( dutsxaHisto.FindBin(xpk-1) ) ); // BG
      dutsxaHisto.Fit( "fgp0", "q", "", xpk-1, xpk+1  );
      cout << "Fit Gauss + BG:"
	   << endl << "  A " << fgp0->GetParameter(0)
	   << endl << "mid " << fgp0->GetParameter(1)
	   << endl << "sig " << fgp0->GetParameter(2)
	   << endl << " BG " << fgp0->GetParameter(3)
	   << endl;
      newDUTalignx = fgp0->GetParameter(1);
      //delete fgp0;
    }

    // y:

    if( dutdyaHisto.GetMaximum() > dutsyaHisto.GetMaximum() ) {
      cout << endl << dutdyaHisto.GetTitle()
	   << " bin " << dutdyaHisto.GetBinWidth(1)
	   << endl;
      TF1 * fgp0 = new TF1( "fgp0", "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]", -1, 1 );
      double xpk = dutdyaHisto.GetBinCenter( dutdyaHisto.GetMaximumBin() );
      fgp0->SetParameter( 0, dutdyaHisto.GetMaximum() ); // amplitude
      fgp0->SetParameter( 1, xpk );
      fgp0->SetParameter( 2, dutdyaHisto.GetBinWidth(1) ); // sigma
      fgp0->SetParameter( 3, dutdyaHisto.GetBinContent( dutdyaHisto.FindBin(xpk-1) ) ); // BG
      dutdyaHisto.Fit( "fgp0", "q", "", xpk-1, xpk+1 );
      cout << "Fit Gauss + BG:"
	   << endl << "  A " << fgp0->GetParameter(0)
	   << endl << "mid " << fgp0->GetParameter(1)
	   << endl << "sig " << fgp0->GetParameter(2)
	   << endl << " BG " << fgp0->GetParameter(3)
	   << endl;
      newDUTaligny = fgp0->GetParameter(1);
      //delete fgp0;
    }
    else {
      cout << endl << dutsyaHisto.GetTitle()
	   << " bin " << dutsyaHisto.GetBinWidth(1)
	   << endl;
      TF1 * fgp0 = new TF1( "fgp0", "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]", -1, 1 );
      double xpk = dutsyaHisto.GetBinCenter( dutsyaHisto.GetMaximumBin() );
      fgp0->SetParameter( 0, dutsyaHisto.GetMaximum() ); // amplitude
      fgp0->SetParameter( 1, xpk );
      fgp0->SetParameter( 2, dutsyaHisto.GetBinWidth(1) ); // sigma
      fgp0->SetParameter( 3, dutsyaHisto.GetBinContent( dutsyaHisto.FindBin(xpk-1) ) ); // BG
      dutsyaHisto.Fit( "fgp0", "q", "", xpk-1, xpk+1 );
      cout << "Fit Gauss + BG:"
	   << endl << "  A " << fgp0->GetParameter(0)
	   << endl << "mid " << fgp0->GetParameter(1)
	   << endl << "sig " << fgp0->GetParameter(2)
	   << endl << " BG " << fgp0->GetParameter(3)
	   << endl;
      newDUTaligny = fgp0->GetParameter(1);
      //delete fgp0;
    }

    // finer alignment:

    cout << endl;

    if( DUTaligniteration > 0 && fabs( newDUTalignx - DUTalignx ) < 0.1 ) {

      if( dutdxHisto.GetEntries() > 999 ) {

	cout << "finer x " << dutdxHisto.GetTitle()
	     << " bin " << dutdxHisto.GetBinWidth(1)
	     << endl;
	TF1 * fgp0 = new TF1( "fgp0", "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]", -1, 1 );
	fgp0->SetParameter( 0, dutdxHisto.GetMaximum() ); // amplitude
	fgp0->SetParameter( 1, dutdxHisto.GetBinCenter( dutdxHisto.GetMaximumBin() ) );
	fgp0->SetParameter( 2, 8*dutdxHisto.GetBinWidth(1) ); // sigma
	fgp0->SetParameter( 3, dutdxHisto.GetBinContent(1) ); // BG
	dutdxHisto.Fit( "fgp0", "q" );
	cout << "Fit Gauss + BG:"
	     << endl << "  A " << fgp0->GetParameter(0)
	     << endl << "mid " << fgp0->GetParameter(1)
	     << endl << "sig " << fgp0->GetParameter(2)
	     << endl << " BG " << fgp0->GetParameter(3)
	     << endl;
	newDUTalignx = DUTalignx + fgp0->GetParameter(1);
	//delete fgp0;

      }

      // dxvsy => rot

      if( ( rot90 || fifty ) && dutdxvsy.GetEntries() > 999 ) {

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
	     << ": extra rot " << fdxvsy->GetParameter(1) << endl;
	DUTrot += fdxvsy->GetParameter(1);

      }

      // dxvsx => turn:

      if( DUTaligniteration > 1 &&
	  dutdxvsx.GetEntries() > 999 &&
	  fabs( DUTturn ) > 0.3
	  ) {

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

      // dxvstx => dz:

      if( ( rot90 || fifty ) && dutdxvstx.GetEntries() > 999 ) {

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

    } // x
    else
      cout << "alignx changed by " << newDUTalignx - DUTalignx << " mm"
	   << ", need more iteration"
	   << endl;

    // y:

    cout << endl;

    if( DUTaligniteration > 0 && fabs( newDUTaligny - DUTaligny ) < 0.1 ) {

      if( dutdyHisto.GetEntries() > 999 ) {

	cout << "finer y " << dutdyHisto.GetTitle()
	     << " bin " << dutdyHisto.GetBinWidth(1)
	     << endl;
	TF1 * fgp0 = new TF1( "fgp0", "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]", -1, 1 );
	fgp0->SetParameter( 0, dutdyHisto.GetMaximum() ); // amplitude
	fgp0->SetParameter( 1, dutdyHisto.GetBinCenter( dutdyHisto.GetMaximumBin() ) );
	fgp0->SetParameter( 2, 5*dutdyHisto.GetBinWidth(1) ); // sigma
	fgp0->SetParameter( 3, dutdyHisto.GetBinContent(1) ); // BG
	dutdyHisto.Fit( "fgp0", "q" );
	cout << "Fit Gauss + BG:"
	     << endl << "  A " << fgp0->GetParameter(0)
	     << endl << "mid " << fgp0->GetParameter(1)
	     << endl << "sig " << fgp0->GetParameter(2)
	     << endl << " BG " << fgp0->GetParameter(3)
	     << endl;
	newDUTaligny = DUTaligny + fgp0->GetParameter(1);
	//delete fgp0;

      }

      // dyvsx => rot

      if( !rot90 && !fifty && dutdyvsx.GetEntries() > 999 ) {

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
	     << ": extra rot " << fdyvsx->GetParameter(1) << endl;
	DUTrot -= fdyvsx->GetParameter(1);

      }

      // dyvsy => tilt:

      if( dutdyvsy.GetEntries() > 999 &&
	  fabs( DUTtilt ) > 3 // else unstable
	  ) {

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

      // dyvsty => dz:

      if( !rot90 && !fifty && dutdyvsty.GetEntries() > 999 ) {

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

    } // y
    else
      cout << "aligny changed by " << newDUTaligny - DUTaligny << " mm"
	   << ", need more iteration"
	   << endl;

    // write new DUT alignment:

    ofstream DUTalignFile( DUTalignFileName.str() );

    DUTalignFile << "# DUT alignment for run " << run << endl;
    ++DUTaligniteration;
    DUTalignFile << "iteration " << DUTaligniteration << endl;
    DUTalignFile << "alignx " << newDUTalignx << endl;
    DUTalignFile << "aligny " << newDUTaligny << endl;
    DUTalignFile << "rot " << DUTrot << endl;
    DUTalignFile << "tilt " << DUTtilt << endl;
    DUTalignFile << "turn " << DUTturn << endl;
    DUTalignFile << "dz " << DUTz - zz[3] << endl;

    DUTalignFile.close();

    cout << endl << "wrote DUT alignment iteration " << DUTaligniteration
	 << " to " << DUTalignFileName.str() << endl
	 << "  alignx " << newDUTalignx << endl
	 << "  aligny " << newDUTaligny << endl
	 << "  rot    " << DUTrot << endl
	 << "  tilt   " << DUTtilt << endl
	 << "  turn   " << DUTturn << endl
	 << "  dz     " << DUTz - zz[3] << endl
      ;

  } // stat
  else
    cout << "not enough for alignment" << endl;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // done

  cout << endl << histoFile->GetName() << endl;

  cout << endl;

  return 0;
}
