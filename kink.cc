
// Daniel Pitzl, DESY, Jun 2019
// telescope analysis with eudaq
// kink-based tracking

// make kink
// kink -g geo_2019_06e.dat -p 0.8 -l 10200 36745
// needs alignment from tele

#include "eudaq/FileReader.hh"
#include "eudaq/PluginManager.hh"
#include "../main/lib/plugins/BDAQ53ConverterPlugin.cc"

#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TProfile2D.h>

#include <sstream> // stringstream
#include <fstream> // filestream
#include <map>
#include <set>
#include <cmath>
#include <time.h> // clock_gettime

using namespace std;
using namespace eudaq;

struct pixel {
  int col;
  int row;
  int nn;
};

struct cluster {
  vector <pixel> vpix;
  int size;
  int ncol, nrow;
  double col, row;
  double mindxy;
  int itr;
};

struct triplet {
  double xm;
  double ym;
  double zm;
  double sx;
  double sy;
  double rxy;
  int ghost;
  vector <double> vx;
  vector <double> vy;
};

//------------------------------------------------------------------------------
vector<cluster> getClus( vector <pixel> pb, int fCluCut = 1 ) // 1 = no gap
{
  // returns clusters with pixel coordinates
  // next-neighbour topological clustering (allows fCluCut-1 empty pixels)

  vector <cluster> vc;
  if( pb.size() == 0 ) return vc;

  int * gone = new int[pb.size()];
  for( unsigned i = 0; i < pb.size(); ++i )
    gone[i] = 0;

  unsigned seed = 0;

  while( seed < pb.size() ) {

    // start a new cluster:

    cluster c;
    c.vpix.push_back( pb[seed] );
    gone[seed] = 1;

    // let it grow as much as possible:

    int growing;
    do{
      growing = 0;
      for( unsigned i = 0; i < pb.size(); ++i ) {
        if( !gone[i] ){ // unused pixel
          for( unsigned int p = 0; p < c.vpix.size(); ++p ) { // vpix in cluster so far
            int dr = c.vpix.at(p).row - pb[i].row;
            int dc = c.vpix.at(p).col - pb[i].col;
	    //if( (   dr>=-fCluCut) && (dr<=fCluCut) 
	    //&& (dc>=-fCluCut) && (dc<=fCluCut) ) { // allow diagonal
            if( ( abs(dr) <= fCluCut && dc == 0 ) ||
		( abs(dc) <= fCluCut && dr == 0 ) ) { // only facing neighbours, same resolution
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
	  ++p->nn;
	  ++q->nn;
	}
    }

    c.col = 0;
    c.row = 0;
    int sumnn = 0;
    int minx = 9999;
    int maxx = 0;
    int miny = 9999;
    int maxy = 0;

    for( vector <pixel>::iterator p = c.vpix.begin(); p != c.vpix.end(); ++p ) {

      int nn = max( 1, p->nn ); // neighbours
      sumnn += nn;
      c.col += p->col * nn;
      c.row += p->row * nn;
      if( p->col > maxx ) maxx = p->col;
      if( p->col < minx ) minx = p->col;
      if( p->row > maxy ) maxy = p->row;
      if( p->row < miny ) miny = p->row;

    }

    //cout << "(cluster with " << c.vpix.size() << " pixels)" << endl;

    //c.col /= c.vpix.size();
    //c.row /= c.vpix.size();
    c.col /= sumnn; // weighted cluster center
    c.row /= sumnn;

    c.size = c.vpix.size();
    c.ncol = maxx-minx+1;
    c.nrow = maxy-miny+1;

    c.mindxy = 999;
    c.itr = 0;

    //c.vpix.clear(); // save space

    vc.push_back(c); // add cluster to vector

    // look for a new seed = unused pixel:

    while( ( ++seed < pb.size() ) && gone[seed] );

  } // while over seeds

  delete [] gone;

  return vc; // vector of clusters

} // getClus

//------------------------------------------------------------------------------
int main( int argc, char* argv[] )
{
  cout << endl << "main " << argv[0] << " called with " << argc << " arguments" << endl;

  if( argc < 4 ) {
    cout << "format: kink -g geo_year_mon.dat run" << endl;
    return 1;
  }

  // run number = last arg

  string runnum( argv[argc-1] );
  int run = atoi( argv[argc-1] );

  cout << endl << "run " << run << endl;

  FileReader * reader;
  if( run < 100 )
    reader = new FileReader( runnum.c_str(), "data/run0000$2R$X");
  else if( run < 1000 )
    reader = new FileReader( runnum.c_str(), "data/run000$3R$X");
  else if( run < 10000 )
    reader = new FileReader( runnum.c_str(), "data/run00$4R$X");
  else if( run < 100000 )
    reader = new FileReader( runnum.c_str(), "data/run0$5R$X");
  else
    reader = new FileReader( runnum.c_str(), "data/run$6R$X");

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // further arguments:

  int lev = 999222111; // last event
  string geoFileName( "geo.dat" );
  double mom = 4.8;

  for( int i = 1; i < argc; ++i ) {

    if( !strcmp( argv[i], "-l" ) )
      lev = atoi( argv[++i] ); // last event

    if( !strcmp( argv[i], "-g" ) )
      geoFileName = argv[++i];

    if( !strcmp( argv[i], "-p" ) )
      mom = atof( argv[++i] ); // momentum

  } // argc

  cout << endl << "p " << mom << " GeV" << endl;

  double ff = 4.8/mom;

  const double ang = 0.0024/sqrt(mom); // [rad] beam divergence

  double X0 = 0.7E-3; // Mimosa
  double scat = 0.0136 * sqrt(X0) / mom * ( 1 + 0.038*log(X0) ); // [rad] Mimosa scattering
  cout << endl << "Mimosa scattering " << scat*1E3 << " mrad" << endl;

  //const double log10 = log(10);

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

  if( geoFile.bad() || ! geoFile.is_open() ) {
    cout << "Error opening " << geoFileName << endl;
    return 1;
  }

  cout << endl << "read geometry from " << geoFileName << endl;

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

      if( ipl < 1 || ipl > 6 ) { // Mimosa 1..6
	//cout << "wrong plane number " << ipl << endl;
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
      midx[ipl] = 0.5 * sizex[ipl]; // mid plane
      midy[ipl] = 0.5 * sizey[ipl]; // mid plane
    }

  } // geo scope

  geoFile.close();

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // Mimosa hot pixels:

  ostringstream hotFileName; // output string stream

  hotFileName << "hot_" << run << ".dat";

  ifstream ihotFile( hotFileName.str() );

  set <int> hotset[9];

  cout << endl;
  if( ihotFile.bad() || ! ihotFile.is_open() ) {
    cout << "no " << hotFileName.str() << ", will be created" << endl;
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

      if( ipl < 1 || ipl > 6 ) {
	//cout << "wrong plane number " << ipl << endl;
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

  for( int ipl = 0; ipl < 6; ++ipl )
    cout << ipl << ": hot " << hotset[ipl].size() << endl;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // telescope alignment:

  int aligniteration = 0;
  double alignx[9];
  double aligny[9];
  double alignz[9];
  double rotx[9];
  double roty[9];

  for( int ipl = 0; ipl < 9; ++ipl ) {

    alignx[ipl] = 0.000; // [mm] same sign as dxAB
    aligny[ipl] = 0.000; // [mm] same sign as dy
    alignz[ipl] = 0.000; // [mm]
    rotx[ipl] = 0.0000; // [rad] rot, same     sign dxvsy
    roty[ipl] = 0.0000; // [rad] rot, opposite sign dyvsx

  }

  ostringstream alignFileName; // output string stream

  alignFileName << "align_" << run << ".dat";

  ifstream ialignFile( alignFileName.str() );

  cout << endl;
  if( ialignFile.bad() || ! ialignFile.is_open() ) {
    cout << "no " << alignFileName.str() << ", will bootstrap" << endl;
    cout << endl;
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

    int ipl = 0;

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

      if( ipl < 1 || ipl > 6 ) {
	//cout << "wrong plane number " << ipl << endl;
	continue;
      }

      double val;
      tokenizer >> val;
      if(      tag == shiftx )
	alignx[ipl] = val;
      else if( tag == shifty )
	aligny[ipl] = val;
      else if( tag == shiftz ) {
	alignz[ipl] = val;
	zz[ipl] += val;
      }
      else if( tag == rotxvsy )
	rotx[ipl] = val;
      else if( tag == rotyvsx )
	roty[ipl] = val;

      // anything else on the line and in the file gets ignored

    } // while getline

  } // alignFile

  ialignFile.close();

  if( aligniteration == 0 ) ff *= 3; // wider binning

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // (re-)create root file:

  ostringstream rootFileName; // output string stream

  rootFileName << "kink" << run << ".root";

  TFile* histoFile = new TFile( rootFileName.str(  ).c_str(  ), "RECREATE" );

  // book histos:

  TH1I hdtus( "dtus", "time between events;time between events [us];events", 100, 0, 1000 );
  TH1I hdtms( "dtms", "time between events;time between events [ms];events", 100, 0, 1000 );

  TH1I t1Histo( "t1", "event time;event time [s];events", 100, 0, 1 );
  TH1I t2Histo( "t2", "event time;event time [s];events", 500, 0, 500 );
  TH1I t3Histo( "t3", "event time;event time [s];events", 150, 0, 1500 );
  TH1I t4Histo( "t4", "event time;event time [s];events", 600, 0, 6000 );
  TH1I t5Histo( "t5", "event time;event time [s];events", 600, 0, 60000 );

  TH1I hnpx0[9];
  TH1I hcol0[9];
  TH1I hrow0[9];
  TH1I hcol[9];
  TH1I hrow[9];

  TH1I hncl[9];
  TH1I hsiz[9];
  TH1I hncol[9];
  TH1I hnrow[9];
  TH1I hmindxy[9];
  TH1I hnpx[9];

  for( int ipl = 1; ipl <= 6; ++ipl ) {

    hnpx0[ipl] = TH1I( Form( "npx0%i", ipl ),
		       Form( "%i pixel per event;all pixels;plane %i events", ipl, ipl ),
		       200, 0, 200 );
    hcol0[ipl] = TH1I( Form( "allcol%i", ipl ),
		       Form( "%i all col;col;%i all pixels", ipl, ipl ), 
		       max( 52, nx[ipl]/4 ), 0, nx[ipl] );
    hrow0[ipl] = TH1I( Form( "allrow%i", ipl ),
		       Form( "%i all row;row;%i all pixels", ipl, ipl ),
		       max( 80, ny[ipl]/2 ), 0, ny[ipl] );

    hcol[ipl] = TH1I( Form( "col%i", ipl ),
		      Form( "%i col;col;%i pixels", ipl, ipl ), 
		      max( 52, nx[ipl]/4 ), 0, nx[ipl] );
    hrow[ipl] = TH1I( Form( "row%i", ipl ),
		      Form( "%i row;row;%i pixels", ipl, ipl ),
		      max( 80, ny[ipl]/2 ), 0, ny[ipl] );

    hncl[ipl] = TH1I( Form( "ncl%i", ipl ),
		      Form( "plane %i cluster per event;cluster;plane %i events", ipl, ipl ),
		      51, -0.5, 50.5 );
    hsiz[ipl] = TH1I( Form( "clsz%i", ipl ),
		      Form( "%i cluster size;pixels/cluster;%i clusters", ipl, ipl ),
		      51, -0.5, 50.5 );
    hncol[ipl] = TH1I( Form( "ncol%i", ipl ), 
		       Form( "%i cluster size x;columns/cluster;%i clusters", ipl, ipl ),
		       21, -0.5, 20.5 );
    hnrow[ipl] = TH1I( Form( "nrow%i", ipl ),
		       Form( "%i cluster size y;rows/cluster;%i clusters", ipl, ipl ),
		       21, -0.5, 20.5 );

    hmindxy[ipl] = TH1I( Form( "mindxy%i", ipl ),
			 Form( "%i cluster isolation;distance to next cluster [mm];%i clusters",
			       ipl, ipl ),
			 200, 0, 20 );
    hnpx[ipl] = TH1I( Form( "npx%i", ipl ),
		      Form( "plane %i cool pixel per event;cool pixels / event;plane %i events", ipl, ipl ),
		      200, 0, 200 );

  } // planes

  double dz = zz[2]-zz[1];
  double rr = 10*ang*dz;

  TH1I hdx21( "dx21", "dx 2-1;dx [mm];2-1 cluster pairs", 100, -rr, rr );
  TH1I hdy21( "dy21", "dy 2-1;dy [mm];2-1 cluster pairs", 100, -rr, rr );
  TH1I hdx21c( "dx21c", "dx 2-1;dx [mm];2-1 cluster pairs", 100, -rr, rr );
  TH1I hdy21c( "dy21c", "dy 2-1;dy [mm];2-1 cluster pairs", 100, -rr, rr );

  TH1I hsx21( "sx21", "slope x 2-1;slope x [mrad];2-1 cluster pairs", 100, -10E3*ang, 10E3*ang );
  TH1I hsy21( "sy21", "slope y 2-1;slope y [mrad];2-1 cluster pairs", 100, -10E3*ang, 10E3*ang );

  TH1I hkx( "kx", "kink x;kink x [mrad];ABC", 300, -30E3*scat, 30E3*scat );
  TH1I hky( "ky", "kink y;kink y [mrad];ABC", 300, -30E3*scat, 30E3*scat );
  TH1I hkxs( "kxs", "kink x;kink x [mrad];straight ABC", 300, -30E3*scat, 30E3*scat );
  TH1I hkys( "kys", "kink y;kink y [mrad];straight ABC", 300, -30E3*scat, 30E3*scat );
  TH1I hkxsc( "kxsc", "kink x;kink x [mrad];straight ABC", 300, -30E3*scat, 30E3*scat );
  TH1I hkysc( "kysc", "kink y;kink y [mrad];straight ABC", 300, -30E3*scat, 30E3*scat );

  TH1I hntri( "ntri", "tracks;tracks;events", 20, -0.5, 19.5 );

  TH1I hdx21i( "dx21i", "dx 2-1;dx [mm];2-1 cluster pairs", 100, -rr, rr );
  TH1I hdy21i( "dy21i", "dy 2-1;dy [mm];2-1 cluster pairs", 100, -rr, rr );
  TH1I hdx21ci( "dx21ci", "dx 2-1;dx [mm];2-1 cluster pairs", 100, -rr, rr );
  TH1I hdy21ci( "dy21ci", "dy 2-1;dy [mm];2-1 cluster pairs", 100, -rr, rr );

  TH1I hsx21i( "sx21i", "slope x 2-1;slope x [mrad];2-1 cluster pairs", 100, -10E3*ang, 10E3*ang );
  TH1I hsy21i( "sy21i", "slope y 2-1;slope y [mrad];2-1 cluster pairs", 100, -10E3*ang, 10E3*ang );

  TH1I hkxi( "kxi", "kink x;kink x [mrad];ABC", 300, -30E3*scat, 30E3*scat );
  TH1I hkyi( "kyi", "kink y;kink y [mrad];ABC", 300, -30E3*scat, 30E3*scat );
  TH1I hkxsi( "kxsi", "kink x;kink x [mrad];straight ABC", 300, -30E3*scat, 30E3*scat );
  TH1I hkysi( "kysi", "kink y;kink y [mrad];straight ABC", 300, -30E3*scat, 30E3*scat );
  TH1I hkxsci( "kxsci", "kink x;kink x [mrad];straight ABC", 300, -30E3*scat, 30E3*scat );
  TH1I hkysci( "kysci", "kink y;kink y [mrad];straight ABC", 300, -30E3*scat, 30E3*scat );

  TH1I hkxs4( "kxs4", "kink x;kink x [mrad];straight ABCD", 300, -30E3*scat, 30E3*scat );
  TH1I hkys4( "kys4", "kink y;kink y [mrad];straight ABCD", 300, -30E3*scat, 30E3*scat );

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // event loop:

  timespec ts;
  clock_gettime( CLOCK_REALTIME, &ts );
  long s0 = ts.tv_sec; // seconds since 1.1.1970
  long f0 = ts.tv_nsec; // nanoseconds
  double zeit1 = 0;
  double zeit2 = 0;

  int iev = 0;
  uint64_t evTLU0 = 0;
  const double fTLU = 384E6; // 384 MHz TLU clock
  uint64_t prevTLU = 0;

  do {

    clock_gettime( CLOCK_REALTIME, &ts );
    long s1 = ts.tv_sec; // seconds since 1.1.1970
    long f1 = ts.tv_nsec; // nanoseconds

    // Get next event:
    DetectorEvent evt = reader->GetDetectorEvent();

    if( evt.IsBORE() ) {
      eudaq::PluginManager::Initialize(evt);
      continue;
    }

    uint64_t evTLU = evt.GetTimestamp(); // 384 MHz = 2.6 ns

    if( iev == 0 )
      evTLU0 = evTLU;

    bool ldbg = 0;

    if( lev < 99 )
      ldbg = 1;

    double evsec = (evTLU - evTLU0) / fTLU;

    t1Histo.Fill( evsec );
    t2Histo.Fill( evsec );
    t3Histo.Fill( evsec );
    t4Histo.Fill( evsec );
    t5Histo.Fill( evsec );

    double evdt = (evTLU - prevTLU) / fTLU;
    hdtus.Fill( evdt * 1E6 ); // [us]
    hdtms.Fill( evdt * 1E3 ); // [ms]
    prevTLU = evTLU;

    if( iev < 10 )
      cout << "kink processing  " << run << "." << iev << "  taken " << evsec << endl;
    else if( iev < 100 && iev%10 == 0 )
      cout << "kink processing  " << run << "." << iev << "  taken " << evsec << endl;
    else if( iev < 1000 && iev%100 == 0 )
      cout << "kink processing  " << run << "." << iev << "  taken " << evsec << endl;
    else if( iev%1000 == 0 )
      cout << "kink processing  " << run << "." << iev << "  taken " << evsec << endl;

    StandardEvent sevt = eudaq::PluginManager::ConvertToStandard(evt);

    string MIM("MIMOSA26");
    int mpl = 1; // Mimosa planes

    vector <cluster> cl[9]; // Mimosa planes

    for( size_t ipl = 0; ipl < sevt.NumPlanes(); ++ipl ) {

      const eudaq::StandardPlane &plane = sevt.GetPlane(ipl);

      if( ldbg )
	std::cout
	  << "plane " << plane.ID()
	  << " " << plane.Type()
	  << " " << plane.Sensor()
	  << " frames " << plane.NumFrames()
	  << " pivot " << plane.PivotPixel()
	  << " total " << plane.TotalPixels()
	  << " hits " << plane.HitPixels()
	  ;

      if( plane.Sensor() != MIM ) continue;

      std::vector<double> pxl = plane.GetPixels<double>();

      hnpx0[mpl].Fill( pxl.size() );

      vector <pixel> pb; // for clustering

      for( size_t ipix = 0; ipix < pxl.size(); ++ipix ) {

	if( ldbg ) 
	  std::cout << ": " << plane.GetX(ipix)
		    << " " << plane.GetY(ipix)
		    << " " << plane.GetPixel(ipix);

	int ix = plane.GetX(ipix); // col pixel index
	int iy = plane.GetY(ipix); // row pixel index
	
	hcol0[mpl].Fill( ix );
	hrow0[mpl].Fill( iy );

	int ipx = ix*ny[mpl] + iy;

	// fill pixel block for clustering:

	if( hotset[mpl].count(ipx) )
	  continue; // skip hot

	hcol[mpl].Fill( ix );
	hrow[mpl].Fill( iy );

	pixel px;
	px.col = ix;
	px.row = iy;
	px.nn = 1; // init neighbours for best resolution
	pb.push_back(px);

	if( pb.size() == 999 ) {
	  cout << "pixel buffer overflow in plane " << mpl
	       << ", event " << iev
	       << endl;
	  break;
	}

      } // pix

      hnpx[mpl].Fill( pb.size() );

      if( ldbg ) std::cout << std::endl;

      // clustering:

      cl[mpl] = getClus( pb );

      if( ldbg ) cout << "A clusters " << cl[mpl].size() << endl;

      hncl[mpl].Fill( cl[mpl].size() );

      for( vector<cluster>::iterator cA = cl[mpl].begin(); cA != cl[mpl].end(); ++cA ) {

	hsiz[mpl].Fill( cA->size );
	hncol[mpl].Fill( cA->ncol );
	hnrow[mpl].Fill( cA->nrow );

	// cluster isolation:

	double mindxy = 99.9;

	for( vector<cluster>::iterator cD = cl[mpl].begin(); cD != cl[mpl].end(); ++cD ) {
	  if( cD == cA ) continue; // self
	  double dx = (cD->col - cA->col)*ptchx[mpl];
	  double dy = (cD->row - cA->row)*ptchy[mpl];
	  double dxy = sqrt( dx*dx + dy*dy );
	  if( dxy < mindxy ) mindxy = dxy;
	}
	cA->mindxy = mindxy;
	hmindxy[mpl].Fill( mindxy ); // mean 4.5 mm

      } // cl A

      ++mpl;
      if( mpl > 6 ) break;

    } // planes

    clock_gettime( CLOCK_REALTIME, &ts );
    long s2 = ts.tv_sec; // seconds since 1.1.1970
    long f2 = ts.tv_nsec; // nanoseconds
    zeit1 += s2 - s1 + ( f2 - f1 ) * 1e-9; // read and cluster

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // make tracks, mark used clusters:

    double zA = zz[1] + alignz[1];
    double zB = zz[2] + alignz[2];
    double zC = zz[3] + alignz[3];

    int itr = 0;

    for( vector<cluster>::iterator cA = cl[1].begin(); cA != cl[1].end(); ++cA ) {

      double xA = cA->col*ptchx[1] - alignx[1];
      double yA = cA->row*ptchy[1] - aligny[1];
      double xmid = xA - midx[1];
      double ymid = yA - midy[1];
      xA = xmid - ymid*rotx[1];
      yA = ymid + xmid*roty[1];

      for( vector<cluster>::iterator cB = cl[2].begin(); cB != cl[2].end(); ++cB ) {

	double xB = cB->col*ptchx[2] - alignx[2];
	double yB = cB->row*ptchy[2] - aligny[2];
	double xmid = xB - midx[2];
	double ymid = yB - midy[2];
	xB = xmid - ymid*rotx[2];
	yB = ymid + xmid*roty[2];

	double dx = xB - xA;
	double dy = yB - yA;
	double dz = zB - zA;

	hdx21.Fill( dx );
	hdy21.Fill( dy );

	if( fabs(dy) < 3 * ang * dz ) // beam divergence
	  hdx21c.Fill( dx );

	if( fabs(dx) < 3 * ang * dz ) // beam divergence
	  hdy21c.Fill( dy );

	double sxAB = dx/dz;
	double syAB = dy/dz;

	hsx21.Fill( sxAB*1E3 );
	hsy21.Fill( syAB*1E3 );

	bool straight = 1;
	if( fabs(sxAB) > 2*ang ) straight = 0;
	if( fabs(syAB) > 2*ang ) straight = 0;

	for( vector<cluster>::iterator cC = cl[3].begin(); cC != cl[3].end(); ++cC ) {

	  double xC = cC->col*ptchx[3] - alignx[3];
	  double yC = cC->row*ptchy[3] - aligny[3];
	  double xmid = xC - midx[3];
	  double ymid = yC - midy[3];
	  xC = xmid - ymid*rotx[3];
	  yC = ymid + xmid*roty[3];

	  double dx = xC - xB;
	  double dy = yC - yB;
	  double dz = zC - zB;

	  double sxBC = dx/dz;
	  double syBC = dy/dz;

	  double kx = sxBC - sxAB;
	  double ky = syBC - syAB;

	  hkx.Fill( kx*1E3 );
	  hky.Fill( ky*1E3 );

	  if( straight ) {

	    hkxs.Fill( kx*1E3 );
	    hkys.Fill( ky*1E3 );
	    if( fabs(ky) < 3*scat )
	      hkxsc.Fill( kx*1E3 );
	    if( fabs(kx) < 3*scat )
	      hkysc.Fill( ky*1E3 );

	    if( fabs(ky) < 3*scat && fabs(kx) < 3*scat ) {
	      ++itr;
	      cA->itr = itr;
	      cB->itr = itr;
	      cC->itr = itr;
	    }

	  } // straight AB

	} // cC

      } // cB

    } // cA

    hntri.Fill( itr );

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // kink tracking:

    double zD = zz[4] + alignz[4];

    for( vector<cluster>::iterator cA = cl[1].begin(); cA != cl[1].end(); ++cA ) {

      double xA = cA->col*ptchx[1] - alignx[1];
      double yA = cA->row*ptchy[1] - aligny[1];
      double xmid = xA - midx[1];
      double ymid = yA - midy[1];
      xA = xmid - ymid*rotx[1];
      yA = ymid + xmid*roty[1];

      for( vector<cluster>::iterator cB = cl[2].begin(); cB != cl[2].end(); ++cB ) {

	if( cB->itr != 0 && cB->itr != cA->itr ) continue;

	double xB = cB->col*ptchx[2] - alignx[2];
	double yB = cB->row*ptchy[2] - aligny[2];
	double xmid = xB - midx[2];
	double ymid = yB - midy[2];
	xB = xmid - ymid*rotx[2];
	yB = ymid + xmid*roty[2];

	double dx = xB - xA;
	double dy = yB - yA;
	double dz = zB - zA;

	hdx21i.Fill( dx );
	hdy21i.Fill( dy );

	if( fabs(dy) < 3 * ang * dz ) // beam divergence
	  hdx21ci.Fill( dx );

	if( fabs(dx) < 3 * ang * dz ) // beam divergence
	  hdy21ci.Fill( dy );

	double sxAB = dx/dz;
	double syAB = dy/dz;

	hsx21i.Fill( sxAB*1E3 );
	hsy21i.Fill( syAB*1E3 );

	bool straight = 1;
	if( fabs(sxAB) > 2*ang ) straight = 0;
	if( fabs(syAB) > 2*ang ) straight = 0;

	for( vector<cluster>::iterator cC = cl[3].begin(); cC != cl[3].end(); ++cC ) {

	  double xC = cC->col*ptchx[3] - alignx[3];
	  double yC = cC->row*ptchy[3] - aligny[3];
	  double xmid = xC - midx[3];
	  double ymid = yC - midy[3];
	  xC = xmid - ymid*rotx[3];
	  yC = ymid + xmid*roty[3];

	  double dx = xC - xB;
	  double dy = yC - yB;
	  double dz = zC - zB;

	  double sxBC = dx/dz;
	  double syBC = dy/dz;

	  double kxB = sxBC - sxAB;
	  double kyB = syBC - syAB;

	  hkxi.Fill( kxB*1E3 );
	  hkyi.Fill( kyB*1E3 );

	  if( straight ) {

	    hkxsi.Fill( kxB*1E3 );
	    hkysi.Fill( kyB*1E3 );
	    if( fabs(kyB) < 3*scat )
	      hkxsci.Fill( kxB*1E3 );
	    if( fabs(kxB) < 3*scat )
	      hkysci.Fill( kyB*1E3 );

	    for( vector<cluster>::iterator cD = cl[4].begin(); cD != cl[4].end(); ++cD ) {

	      double xD = cD->col*ptchx[4] - alignx[4];
	      double yD = cD->row*ptchy[4] - aligny[4];
	      double xmid = xD - midx[4];
	      double ymid = yD - midy[4];
	      xD = xmid - ymid*rotx[4];
	      yD = ymid + xmid*roty[4];

	      double dx = xD - xC;
	      double dy = yD - yC;
	      double dz = zD - zC;

	      double sxCD = dx/dz;
	      double syCD = dy/dz;

	      double kxC = sxCD - sxBC;
	      double kyC = syCD - syBC;

	      if( fabs(kxC) < 3*scat &&
		  fabs(kyC) < 3*scat ) {

		hkxs4.Fill( kxB*1E3 );
		hkys4.Fill( kyB*1E3 );

	      }

	    } // cD

	  } // straight AB

	} // cC

      } // cB

    } // cA

    clock_gettime( CLOCK_REALTIME, &ts );
    long s3 = ts.tv_sec; // seconds since 1.1.1970
    long f3 = ts.tv_nsec; // nanoseconds
    zeit2 += s3 - s2 + ( f3 - f2 ) * 1e-9; // read and cluster

    ++iev;

  } while( reader->NextEvent() && iev < lev );

  delete reader;

  clock_gettime( CLOCK_REALTIME, &ts );
  long s9 = ts.tv_sec; // seconds since 1.1.1970
  long f9 = ts.tv_nsec; // nanoseconds

  cout << "done after " << iev << " events"
       << " in " << s9 - s0 + ( f9 - f0 ) * 1e-9 << " s"
       << " (read and cluster " << zeit1 << " s, tracking " << zeit2 << " s)"
       << endl;

  histoFile->Write();
  histoFile->Close();

  cout << endl << histoFile->GetName() << endl;

  cout << endl;

  return 0;
}
