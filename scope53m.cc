
// Daniel Pitzl, DESY, Jun 2018
// telescope analysis with RD53A
// module in front: only triplets

// make scope53
// scope53 33095
// needs runs.dat
// needs geo_201x_yy.dat
// needs align_33095.dat from tele
// uses hot_33485.dat from tele
// uses alignDUT_33485.dat
// uses alignMOD_33485.dat

#include "eudaq/FileReader.hh"
#include "eudaq/PluginManager.hh"
#include "../main/lib/plugins/BDAQ53ConverterPlugin.cc"

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
  int frm;
  bool pivot;
};

struct cluster {
  vector <pixel> vpix; // Armin Burgmeier: list
  int size;
  int ncol, nrow, nfrm;
  double col, row;
  int signal;
  double mindxy;
};

struct triplet {
  double xm;
  double ym;
  double zm;
  double sx;
  double sy;
  bool lk;
  double ttdmin;
  unsigned iA;
  unsigned iB;
  unsigned iC;
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
    int minf = 999;
    int maxf = 0;

    for( vector<pixel>::iterator p = c.vpix.begin(); p != c.vpix.end(); ++p ) {

      int nn = max( 1, p->tot ); // neighbours
      sumnn += nn;
      c.col += p->col * nn;
      c.row += p->row * nn;

      if( p->col > maxx ) maxx = p->col;
      if( p->col < minx ) minx = p->col;
      if( p->row > maxy ) maxy = p->row;
      if( p->row < miny ) miny = p->row;
      if( p->frm > maxf ) maxf = p->frm;
      if( p->frm < minf ) minf = p->frm;

    }

    //cout << "(cluster with " << c.vpix.size() << " pixels)" << endl;

    c.col /= sumnn;
    c.row /= sumnn;
    c.signal = sumnn;
    c.ncol = maxx-minx+1;
    c.nrow = maxy-miny+1;
    c.nfrm = maxf-minf+1;
    c.mindxy = 999;

    vc.push_back(c); // add cluster to vector

    // look for a new seed = used pixel:

    while( ( ++seed < pb.size() ) && gone[seed] );

  } // while over seeds

  delete[] gone;

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
    int minf = 999;
    int maxf = 0;

    for( vector<pixel>::iterator p = c.vpix.begin(); p != c.vpix.end(); ++p ) {

      double Qpix = p->tot;

      sumQ += Qpix;

      c.col += p->col*Qpix;
      c.row += p->row*Qpix;

      if( p->col > maxx ) maxx = p->col;
      if( p->col < minx ) minx = p->col;
      if( p->row > maxy ) maxy = p->row;
      if( p->row < miny ) miny = p->row;
      if( p->frm > maxf ) maxf = p->frm;
      if( p->frm < minf ) minf = p->frm;

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
    c.nfrm = maxf-minf+1;
    c.mindxy = 999;

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

  int fev = 0; // 1st event
  int lev = 999222111; // last event

  for( int i = 1; i < argc; ++i ) {

    if( !strcmp( argv[i], "-f" ) )
      fev = atoi( argv[++i] ); // 1st event

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
    string TURN( "turn" );
    string FIFTY( "fifty" );
    string ROT90( "rot90" );
    bool found = 0;

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
	<< "  beam " << pbeam << " GeV" << endl
	<< "  geo file " << geoFileName << endl
	<< "  nominal DUT turn " << DUTturn0 << " deg" << endl
	<< "  DUT chip " << chip0 << endl
	<< "  fifty " << fifty << endl
	<< "  rot90 " << rot90 << endl
	<< "  modrun " << modrun << endl
	;
    else {
      cout << "run " << run << " not found in runs.dat" << endl;
      return 1;
    }

  } // runsFile

  runsFile.close();

  const double fTLU = 384E6; // 384 MHz TLU clock
  /*
    effvsxy->Draw("colz")
    TArc c( 3.0, -0.5, 3.0 ); // circle
    c.SetFillStyle(0);
    c.SetLineWidth(3);
    c.Draw("same");
    //c.SetY1(-0.5);
    c.SetR1(3.5); // ax
    c.SetR2(2.5); // ay
    TArc m( 3.0, -0.5, 0.1 ); // mid
    m.SetFillStyle(1000);
    m.SetFillColor(1);
    m.Draw("same");
  */

  double xbeam =  3.0; // 509 Lin
  double ybeam = -0.5;
  double rbeam = 3.0;

  if( run >= 35150 ) {
    xbeam = -3.0; // 509 Syn
    ybeam = -0.5;
    rbeam = 4.0;
  }
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

  TFile* histoFile = new TFile( rootFileName.str(  ).c_str(  ), "RECREATE" );

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

  if( DUTaligniteration <= 1 ) {
    DUTtilt = DUTtilt0;
    DUTturn = DUTturn0;
  }

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

  double qL = 10; // 33483
  double qR = 20;

  if( run >= 33314 && run <= 33333 ) { // 501  KRUM_CURR 10
    qL = 15;
    qR = 29;
  }

  if( run >= 33475 && run <= 33483 ) { // 504  KRUM_CURR 30
    qL = 10;
    qR = 20;
  }

  if( run >= 34001 ) { // 524  KRUM_CURR 32
    qL =  7;
    qR = 20;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // DUT dead pixels:

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
	    if( col%2 )
	      iy = 2*row + 0; // sensor 25
	    else
	      iy = 2*row + 1;
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

  if( iDUThotFile.bad() || ! iDUThotFile.is_open() ) {
    cout << "no " << DUThotFileName.str() << endl;
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
	  if( ix%2 )
	    row = 2*iy + 0; // sensor 25
	  else
	    row = 2*iy + 1;
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
  cout << "try to open  " << modfileA;
  ifstream Astream( modfileA.c_str() );
  if( !Astream ) {
    cout << " : failed " << endl;
    modrun = 0; // flag
  }
  cout << " : succeed " << endl;

  string modfileB = "mod/run" + to_string(modrun) + "B.out";
  cout << "try to open  " << modfileB;
  ifstream Bstream( modfileB.c_str() );
  if( !Bstream ) {
    cout << " : failed " << endl;
    modrun = 0;
  }
  cout << " : succeed " << endl;

  string sl;
  int mev = 0;
  if( fev ) cout << "MOD skip " << fev << endl;
  while( mev < fev ) {
    getline( Astream, sl ); // fast forward
    getline( Bstream, sl );
    ++mev;
  }

  if( modrun == 100 ) { // Jun 2018
    getline( Astream, sl ); // read one line for sync
    getline( Bstream, sl );
  }

  if( modrun == 106 ) { // Jun 2018
    getline( Astream, sl ); // read one line for sync
    getline( Bstream, sl );
  }

  if( modrun == 116 ) { // Jun 2018
    getline( Astream, sl ); // read one line for sync
    getline( Bstream, sl );
  }

  if( modrun == 163 ) { // Dec 2018
    getline( Astream, sl ); // read one line for sync
    getline( Bstream, sl );
    getline( Astream, sl ); // read one line for sync
    getline( Bstream, sl );
  }

  if( modrun >= 244 && modrun <= 249 ) { // Feb 2019
    getline( Astream, sl ); // read one line for sync
    getline( Bstream, sl );
    getline( Astream, sl ); // read one line for sync
    getline( Bstream, sl );
  }

  if( modrun >= 252 && modrun <= 256 ) { // Feb 2019
    getline( Astream, sl ); // read one line for sync
    getline( Bstream, sl );
  }

  int iMOD = 7;

  int MODaligniteration = 0;
  double MODalignx = -2.5;
  double MODaligny = -0.4;
  double MODrot = 0.007;
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
  // book histos:

  double f = 4.8/pbeam;

  TH1I t1Histo( "t1", "event time;event time [s];events / 10 ms", 100, 0, 1 );
  TH1I t2Histo( "t2", "event time;event time [s];events / s", 300, 0, 300 );
  TH1I t3Histo( "t3", "event time;event time [s];events / 10 s", 150, 0, 1500 );
  TH1I t4Histo( "t4", "event time;event time [s];events /10 s", 600, 0, 6000 );
  TH1I t5Histo( "t5", "event time;event time [s];events / 100 s", 600, 0, 60000 );
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
		      Form( "%i pixel per event;pixels;plane %i events", ipl, ipl ),
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

  hmap[8] = new TH2I( Form( "map%i", 8 ),
		      Form( "%i map;col;row;plane %i pixels", 8, 8 ),
		      400, 0, 400, 192, 0, 192 );

  TH1I dutpxcol0Histo( "dutpxcol0",
		      "DUT pixel column before threshold;DUT pixel column;pixels",
		      400, 0, 400 );

  TH1I synpxqHisto( "synpxq",
		    "Sync pixel signal;Sync pixel signal [ToT];Sync pixels",
		    16, -0.5, 15.5 );
  TH1I linpxqHisto( "linpxq",
		    "Lin pixel signal;Lin pixel signal [ToT];Lin pixels",
		    16, -0.5, 15.5 );
  TH1I difpxqHisto( "difpxq",
		    "Diff pixel signal;Diff pixel signal [ToT];Diff pixels",
		    16, -0.5, 15.5 );

  TProfile dutpxqvsx( "dutpxqvsx",
		      "DUT pixel signal vs x;column;<pixel signal> [ToT]",
		      400, 0, 400, -1, 16 );
  TProfile2D * dutpxqvsxy = new
    TProfile2D( "dutpxqvsxy",
		"DUT pixel signal map;column;row;<pixel signal> [ToT]",
		400, 0, 400, 192, 0, 192, -1, 16 );

  TH1I dutpxbcHisto( "dutpxbc",
		    "DUT pixel BC;DUT pixel BC;DUT pixels without masking",
		    32, -0.5, 31.5 );
  TH1I synpxbcHisto( "synpxbc",
		    "Sync pixel BC;Sync pixel BC;Sync pixels",
		    32, -0.5, 31.5 );
  TH1I linpxbcHisto( "linpxbc",
		    "Lin pixel BC;Lin pixel BC;Lin pixels",
		    32, -0.5, 31.5 );
  TH1I difpxbcHisto( "difpxbc",
		    "Diff pixel BC;Diff pixel BC;Diff pixels",
		    32, -0.5, 31.5 );
  TProfile dutpxbcvsx( "dutpxbcvsx",
		      "DUT pixel BC vs x;column;<pixel BC>",
		      400, 0, 400, -1, 32 );
  TProfile2D * dutpxbcvsxy = new
    TProfile2D( "dutpxbcvsxy",
		"DUT pixel BC map;column;row;<pixel BC>",
		400, 0, 400, 192, 0, 192, -1, 32 );
  TProfile dutpxbcvst5( "dutpxbcvst5",
			"DUT pixel BC vs time;time [s];<pixel BC>",
			1100, 0, 66000, -0.5, 31.5 );

  TH1I dutpxcol9Histo( "dutpxcol9",
		      "DUT pixel column after threshold;DUT pixel column;pixels",
		      400, 0, 400 );

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
  TH1I tritxHisto( "tritx", "triplet slope x;slope x [rad];triplets",
		   100, -0.005*f, 0.005*f );
  TH1I trityHisto( "trity", "triplet slope y;slope y [rad];triplets",
		   100, -0.005*f, 0.005*f );

  TH1I ntriHisto( "ntri", "triplets;triplets;events", 51, -0.5, 50.5 );

  // MOD vs triplets:

  TH1I ttdminmod1Histo( "ttdminmod1",
		     "triplet isolation at MOD;triplet at MOD min #Delta_{xy} [mm];triplet pairs",
		     100, 0, 1 );
  TH1I ttdminmod2Histo( "ttdminmod2",
		     "triplet isolation at MOD;triplet at MOD min #Delta_{xy} [mm];triplet pairs",
		     150, 0, 15 );

  TH1I modxHisto( "modx",
		  "MOD x;MOD cluster x [mm];MOD clusters",
		  700, -35, 35 );
  TH1I modyHisto( "mody",
		  "MOD y;MOD cluster y [mm];MOD clusters",
		  200, -10, 10 );

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
  TH1I moddxcqHisto( "moddxcq",
		     "MOD - triplet x Landau peak;MOD cluster - triplet #Deltax [mm];Landau peak MOD clusters",
		     500, -0.5, 0.5 );
  TProfile moddxvsx( "moddxvsx",
		     "MOD #Deltax vs x;x track [mm];<cluster - triplet #Deltax> [mm]",
		     216, -32.4, 32.4, -2.5, 2.5 );
  TProfile moddxvsy( "moddxvsy",
		     "MOD #Deltax vs y;y track [mm];<cluster - triplet #Deltax> [mm]",
		     160, -8, 8, -2.5, 2.5 );
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

  TH1I modlkxBHisto( "modlkxb",
		     "linked triplet at MOD x;triplet x at MOD [mm];linked triplets",
		     216, -32.4, 32.4 );
  TH1I modlkyBHisto( "modlkyb",
		     "linked triplet at MOD y;triplet y at MOD [mm];linked triplets",
		     160, -8, 8 );
  TH1I modlkxHisto( "modlkx",
		    "linked triplet at MOD x;triplet x at MOD [mm];linked triplets",
		    216, -32.4, 32.4 );
  TH1I modlkyHisto( "modlky",
		    "linked triplet at MOD y;triplet y at MOD [mm];linked triplets",
		    160, -8, 8 );

  TH1I modlkcolHisto( "modlkcol",
		      "MOD linked col;MOD linked col;linked MOD cluster",
		      216, 0, 432 );
  TH1I modlkrowHisto( "modlkrow",
		      "MOD linked row;MOD linked row;linked MOD cluster",
		      182, 0, 182 );

  TH1I trinfrmlkHisto( "trinfrmlk",
		     "triplet linked frames;frames;linked triplet cluster",
		     2, 0.5, 2.5 );
  TH1I tripxfrmlkHisto( "tripxfrmlk",
			"triplet linked pixel frame;frame;linked triplet pixels",
			2, -0.5, 1.5 );
  TH1I tripxpivlkHisto( "tripxpivlk",
			"triplet linked pixel pivot;pivot;linked triplet pixels",
			2, -0.5, 1.5 );

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
		      26000, 0, 26E6, -0.5, 1.5 );
  TProfile modlkvsev1( "modlkvsev1",
		      "triplet-MOD links vs events;events;triplets with MOD links / 100",
		      100, 755e3, 765E3, -0.5, 1.5 );
  TProfile modlkvsev2( "modlkvsev2",
		      "triplet-MOD links vs events;events;triplets with MOD links / 100",
		      100, 570e3, 580E3, -0.5, 1.5 );

  TH1I ntrimodHisto( "ntrimod", "triplet - MOD links;triplet - MOD links;events",
		    11, -0.5, 10.5 );

  // DUT clusters:

  TH2I * dutmodxxHisto = new
    TH2I( "dutmodxx", "Mod vs DUT x-x;x_{DUT} [mm];x_{MOD} [mm];cluster pairs",
	  220, -11, 11, 660, -33, 33 );
  TH2I * dutmodyyHisto = new
    TH2I( "dutmodyy", "Mod vs DUT y-y;y_{DUT} [mm];y_{MOD} [mm];cluster pairs",
	  120, -6, 6, 160, -8, 8 );

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
		   500, -0.5, 0.5 );
  TH1I dutdyHisto( "dutdy",
		   "DUT - track dy;DUT cluster - track #Deltay [mm];DUT clusters",
		   500, -0.5, 0.5 );

  TH2I * dutxxHisto = new
    TH2I( "dutxx", "tracks vs DUT in x;track x [mm];DUT cluster x [mm];track-cluster pairs",
	  500, -5, 5, 500, -5, 5 );

  TH1I dutdxcHisto( "dutdxc",
		    "DUT - track dx;DUT cluster - track #Deltax [mm];DUT clusters",
		    500, -0.5, 0.5 );
  TH1I lindxcHisto( "lindxc",
		    "Lin - track dx;Lin cluster - track #Deltax [mm];Lin clusters",
		    200, -0.2, 0.2 );
  TH1I difdxcHisto( "difdxc",
		    "Dif - track dx;Dif cluster - track #Deltax [mm];Dif clusters",
		    200, -0.2, 0.2 );
  TH1I dutdxc1Histo( "dutdxc1",
		     "DUT - track dx;DUT cluster - track #Deltax [mm];DUT 1-px clusters",
		     200, -0.2, 0.2 );
  TH1I dutdxc2Histo( "dutdxc2",
		     "DUT - track dx;DUT cluster - track #Deltax [mm];DUT 2-px clusters",
		     200, -0.2, 0.2 );

  double limx = 0.1;
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
  TProfile dutdxvst2( "dutdxvst2",
		      "DUT #Deltax vs time;time [s];<DUT #Deltax> [mm/5s]",
		      200, 0, 1000, -limx, limx );
  TProfile dutdxvst5( "dutdxvst5",
		      "DUT #Deltax vs time;time [s];<DUT #Deltax> [mm]",
		      1000, 0, 60000, -limx, limx );

  TProfile dutmadxvsx( "dutmadxvsx",
		       "DUT MAD(#Deltax) vs x;x track [mm];MAD(#Deltax) [mm]",
		       200, -10, 10, 0, limx );
  TProfile dutmadxvsxm( "dutmadxvsxm",
			"DUT MAD(#Deltax) vs xmod;x track mod 100 [#mum];MAD(#Deltax) [mm]",
			50, 0, 100, 0, limx );
  TProfile dutmadxvstx( "dutmadxvstx",
			"DUT MAD(#Deltax) vs #theta_{x};x track slope [rad];MAD(#Deltax) [mm]",
			80, -0.002, 0.002, 0, limx );
  TProfile dutmadxvsq( "dutmadxvsq",
		       "DUT MAD(#Deltax) vs Q;cluster signal [ToT];MAD(#Deltax) [mm]",
		       80, 0, 80, 0, limx );
  TH1I dutdxcqHisto( "dutdxcq",
		     "DUT - track dx Landau peak;DUT cluster - track #Deltax [mm];Landau peak DUT clusters",
		     500, -0.2, 0.2 );

  TH2I * dutyyHisto = new
    TH2I( "dutyy", "tracks vs DUT in y;track y [mm];DUT cluster y [mm];track-cluster pairs",
	  500, -5, 5, 500, -5, 5 );

  TH1I dutdy0Histo( "dutdy0",
		    "DUT - track dy;DUT cluster - track #Deltay [mm];DUT clusters",
		    500, -5, 5 ); // binning like dutyy
  TH1I dutdycHisto( "dutdyc",
		    "DUT - track dy;DUT cluster - track #Deltay [mm];DUT clusters",
		    500, -0.5, 0.5 );
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
		      1000, 0, 60000, -0.1, 0.1 );

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

  TH1I synqHisto( "synq",
		  "SYN linked clusters;SYN cluster signal [ToT];linked SYN clusters",
		  80, 0, 80 );
  TH1I synq0Histo( "synq0",
		   "SYN linked clusters;SYN normal cluster signal [ToT];linked SYN clusters",
		   80, 0, 80 );
  TProfile synqvsx( "synqvsx",
		    "SYN cluster signal vs x;x track [mm];SYN <cluster signal> [ToT]",
		    200, -10, 10 );
  TProfile synqvsy( "synqvsy",
		    "SYN cluster signal vs y;y track [mm];SYN <cluster signal> [ToT]",
		    100, -5, 5 );
  TProfile2D * synqvsxy = new
    TProfile2D( "synqvsxy",
		"SYN cluster signal vs xy;x track [mm];y track [mm];SYN <cluster signal> [ToT]",
		200, -10, 10, 100, -5, 5 );
  TProfile synqvsr( "synqvsr",
		    "SYN cluster signal vs Rbeam;R track from beam [mm];SYN <cluster signal> [ToT]",
		    100, 0, 10 );
  TH1I synqbeamHisto( "synqbeam",
		  "SYN linked clusters beam;SYN cluster signal [ToT];linked SYN beam clusters",
		  80, 0, 80 );

  TH1I synnpxHisto( "synnpx",
		    "SYN linked cluster size;SYN cluster size [pixels];linked SYN clusters",
		    20, 0.5, 20.5 );
  TH1I synncolHisto( "synncol",
		     "SYN linked cluster cols;SYN cluster size [columns];linked SYN clusters",
		     20, 0.5, 20.5 );
  TH1I synnrowHisto( "synnrow",
		     "SYN linked cluster rows;SYN cluster size [rows];linked SYN clusters",
		     20, 0.5, 20.5 );

  TH1I synnrow1Histo( "synnrow1",
		      "SYN linked 1-col cluster rows;SYN cluster size [rows];linked SYN 1-col clusters",
		      20, 0.5, 20.5 );
  TH1I synnrow1eveHisto( "synnrow1eve",
			 "SYN linked even seed row 1-col cluster rows;SYN cluster size [rows];linked SYN even seed row 1-col clusters",
			 20, 0.5, 20.5 );
  TH1I synnrow1oddHisto( "synnrow1odd",
			 "SYN linked odd seed row 1-col cluster rows;SYN cluster size [rows];linked SYN odd seed row 1-col clusters",
			 20, 0.5, 20.5 );

  TH1I synnpxfHisto( "synnpxf",
		    "SYN linked 1-BC cluster size;SYN cluster size [pixels];linked 1-BC SYN clusters",
		     20, 0.5, 20.5 );
  TH1I synncolfHisto( "synncolf",
		      "SYN linked 1-BC cluster cols;SYN cluster size [columns];linked 1-BC SYN clusters",
		      20, 0.5, 20.5 );
  TH1I synnrowfHisto( "synnrowf",
		      "SYN linked 1-BC cluster rows;SYN cluster size [rows];linked 2-BC SYN clusters",
		      20, 0.5, 20.5 );

  TH1I synnpxffHisto( "synnpxff",
		      "SYN linked 2-BC cluster size;SYN cluster size [pixels];linked 2-BC SYN clusters",
		      20, 0.5, 20.5 );
  TH1I synncolffHisto( "synncolff",
		       "SYN linked 2-BC cluster cols;SYN cluster size [columns];linked 2-BC SYN clusters",
		       20, 0.5, 20.5 );
  TH1I synnrowffHisto( "synnrowff",
		       "SYN linked 2-BC cluster rows;SYN cluster size [rows];linked 2-BC SYN clusters",
		       20, 0.5, 20.5 );

  TProfile2D * synnpxvsxmym = new
    TProfile2D( "synnpxvsxmym",
		"SYN cluster size vs xmod ymod;x track mod 100 [#mum];y track mod 100 [#mum];SYN <cluster size> [pixels]",
		50, 0, 100, 50, 0, 100, 0, 20 );
  TProfile synncolvsxm( "synncolvsxm",
			"SYN cluster size vs xmod;x track mod 100 [#mum];SYN <cluster size> [columns]",
			50, 0, 100, 0, 20 );
  TProfile synncolvsym( "synncolvsym",
			"SYN cluster size vs ymod;y track mod 100 [#mum];SYN <cluster size> [columns]",
			50, 0, 100, 0, 20 );

  TProfile2D * synnrowvsxmym = new
    TProfile2D( "synnrowvsxmym",
		"SYN cluster size vs xmod ymod;x track mod 100 [#mum];y track mod 100 [#mum];SYN <cluster size> [rows]",
		50, 0, 100, 50, 0, 100, 0, 20 );
  TProfile synnrowvsxm( "synnrowvsxm",
			"SYN cluster size vs xmod;x track mod 100 [#mum];SYN <cluster size> [rows]",
			50, 0, 100, 0, 20 );
  TProfile synnrowvsxm5( "synnrowvsxm5",
			 "SYN cluster size vs xmod;x track mod 100 [#mum];SYN <cluster size> [rows]",
			 50, 0, 50, 0, 20 );
  TProfile synnrowvsym( "synnrowvsym",
			"Syn cluster size vs ymod;y track mod 100 [#mum];Syn <cluster size> [rows]",
			50, 0, 100, 0, 20 );

  TProfile2D * synqvsxmym = new
    TProfile2D( "synqvsxmym",
		"SYN cluster signal vs xmod ymod;x track mod 100 [#mum];y track mod 100 [#mum];SYN <cluster signal> [ToT]",
		50, 0, 100, 50, 0, 100, 0, 200 );
  TProfile synqvsxm( "synqvsxm",
		     "SYN cluster signal vs xmod;x track mod 100 [#mum];SYN <cluster signal> [ToT]",
		     50, 0, 100, 0, 100 );
  TProfile synqvsxm5( "synqvsxm5",
		      "SYN cluster signal vs xmod;x track mod 50 [#mum];SYN <cluster signal> [ToT]",
		      50, 0, 50, 0, 100 );
  TProfile synqvsym( "synqvsym",
		     "SYN cluster signal vs ymod;y track mod 100 [#mum];SYN <cluster signal> [ToT]",
		     50, 0, 100, 0, 100 );

  TH1I linqHisto( "linq",
		  "LIN linked clusters;LIN cluster signal [ToT];linked LIN clusters",
		  80, 0, 80 );
  TH1I linq0Histo( "linq0",
		   "LIN linked clusters;LIN normal cluster signal [ToT];linked LIN clusters",
		   80, 0, 80 );
  TProfile linqvsx( "linqvsx",
		    "LIN cluster signal vs x;x track [mm];LIN <cluster signal> [ToT]",
		    200, -10, 10 );
  TProfile linqvsy( "linqvsy",
		    "LIN cluster signal vs y;y track [mm];LIN <cluster signal> [ToT]",
		    100, -5, 5 );
  TProfile2D * linqvsxy = new
    TProfile2D( "linqvsxy",
		"LIN cluster signal vs xy;x track [mm];y track [mm];LIN <cluster signal> [ToT]",
		200, -10, 10, 100, -5, 5 );
  TProfile linqvsr( "linqvsr",
		    "LIN cluster signal vs Rbeam;R track from beam [mm];LIN <cluster signal> [ToT]",
		    100, 0, 10 );

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
  TProfile linncolvsxm( "linncolvsxm",
			"LIN cluster size vs xmod;x track mod 100 [#mum];LIN <cluster size> [columns]",
			50, 0, 100, 0, 20 );
  TProfile linncolvsym( "linncolvsym",
			"LIN cluster size vs ymod;y track mod 100 [#mum];LIN <cluster size> [columns]",
			50, 0, 100, 0, 20 );

  TProfile2D * linnrowvsxmym = new
    TProfile2D( "linnrowvsxmym",
		"LIN cluster size vs xmod ymod;x track mod 100 [#mum];y track mod 100 [#mum];LIN <cluster size> [rows]",
		50, 0, 100, 50, 0, 100, 0, 20 );
  TProfile linnrowvsxm( "linnrowvsxm",
			"LIN cluster size vs xmod;x track mod 100 [#mum];LIN <cluster size> [rows]",
			50, 0, 100, 0, 20 );
  TProfile linnrowvsxm5( "linnrowvsxm5",
			 "LIN cluster size vs xmod;x track mod 100 [#mum];LIN <cluster size> [rows]",
			 50, 0, 50, 0, 20 );
  TProfile linnrowvsym( "linnrowvsym",
			"Lin cluster size vs ymod;y track mod 100 [#mum];Lin <cluster size> [rows]",
			50, 0, 100, 0, 20 );

  TProfile2D * linqvsxmym = new
    TProfile2D( "linqvsxmym",
		"LIN cluster signal vs xmod ymod;x track mod 100 [#mum];y track mod 100 [#mum];LIN <cluster signal> [ToT]",
		50, 0, 100, 50, 0, 100, 0, 200 );
  TProfile linqvsxm( "linqvsxm",
		     "LIN cluster signal vs xmod;x track mod 100 [#mum];LIN <cluster signal> [ToT]",
		     50, 0, 100, 0, 100 );
  TProfile linqvsxm5( "linqvsxm5",
		      "LIN cluster signal vs xmod;x track mod 50 [#mum];LIN <cluster signal> [ToT]",
		      50, 0, 50, 0, 100 );
  TProfile linqvsym( "linqvsym",
		     "LIN cluster signal vs ymod;y track mod 100 [#mum];LIN <cluster signal> [ToT]",
		     50, 0, 100, 0, 100 );

  TH1I difqHisto( "difq",
		  "DIF linked clusters;DIF cluster signal [ToT];linked DIF clusters",
		  80, 0, 80 );
  TH1I difq0Histo( "difq0",
		   "DIF linked clusters;DIF normal cluster signal [ToT];linked DIF clusters",
		   80, 0, 80 );

  TH1I difnpxHisto( "difnpx",
		    "DIF linked cluster size;DIF cluster size [pixels];linked DIF clusters",
		    20, 0.5, 20.5 );
  TH1I difncolHisto( "difncol",
		     "DIF linked cluster cols;DIF cluster size [columns];linked DIF clusters",
		     20, 0.5, 20.5 );
  TH1I difnrowHisto( "difnrow",
		     "DIF linked cluster rows;DIF cluster size [rows];linked DIF clusters",
		     20, 0.5, 20.5 );

  TH1I difnrow1Histo( "difnrow1",
		      "DIF linked 1-col cluster rows;DIF cluster size [rows];linked DIF 1-col clusters",
		      20, 0.5, 20.5 );
  TH1I difnrow1eveHisto( "difnrow1eve",
			 "DIF linked even seed row 1-col cluster rows;DIF cluster size [rows];linked DIF even seed row 1-col clusters",
			 20, 0.5, 20.5 );
  TH1I difnrow1oddHisto( "difnrow1odd",
			 "DIF linked odd seed row 1-col cluster rows;DIF cluster size [rows];linked DIF odd seed row 1-col clusters",
			 20, 0.5, 20.5 );

  TH1I difnpxfHisto( "difnpxf",
		    "DIF linked 1-BC cluster size;DIF cluster size [pixels];linked 1-BC DIF clusters",
		     20, 0.5, 20.5 );
  TH1I difncolfHisto( "difncolf",
		      "DIF linked 1-BC cluster cols;DIF cluster size [columns];linked 1-BC DIF clusters",
		      20, 0.5, 20.5 );
  TH1I difnrowfHisto( "difnrowf",
		      "DIF linked 1-BC cluster rows;DIF cluster size [rows];linked 2-BC DIF clusters",
		      20, 0.5, 20.5 );

  TH1I difnpxffHisto( "difnpxff",
		      "DIF linked 2-BC cluster size;DIF cluster size [pixels];linked 2-BC DIF clusters",
		      20, 0.5, 20.5 );
  TH1I difncolffHisto( "difncolff",
		       "DIF linked 2-BC cluster cols;DIF cluster size [columns];linked 2-BC DIF clusters",
		       20, 0.5, 20.5 );
  TH1I difnrowffHisto( "difnrowff",
		       "DIF linked 2-BC cluster rows;DIF cluster size [rows];linked 2-BC DIF clusters",
		       20, 0.5, 20.5 );

  TProfile2D * difnpxvsxmym = new
    TProfile2D( "difnpxvsxmym",
		"DIF cluster size vs xmod ymod;x track mod 100 [#mum];y track mod 100 [#mum];DIF <cluster size> [pixels]",
		50, 0, 100, 50, 0, 100, 0, 20 );
  TProfile difncolvsxm( "difncolvsxm",
			"DIF cluster size vs xmod;x track mod 100 [#mum];DIF <cluster size> [columns]",
			50, 0, 100, 0, 20 );
  TProfile difncolvsym( "difncolvsym",
			"DIF cluster size vs ymod;y track mod 100 [#mum];DIF <cluster size> [columns]",
			50, 0, 100, 0, 20 );

  TProfile2D * difnrowvsxmym = new
    TProfile2D( "difnrowvsxmym",
		"DIF cluster size vs xmod ymod;x track mod 100 [#mum];y track mod 100 [#mum];DIF <cluster size> [rows]",
		50, 0, 100, 50, 0, 100, 0, 20 );
  TProfile difnrowvsxm( "difnrowvsxm",
			"DIF cluster size vs xmod;x track mod 100 [#mum];DIF <cluster size> [rows]",
			50, 0, 100, 0, 20 );
  TProfile difnrowvsym( "difnrowvsym",
			"Dif cluster size vs ymod;y track mod 100 [#mum];Dif <cluster size> [rows]",
			50, 0, 100, 0, 20 );

  TProfile2D * difqvsxmym = new
    TProfile2D( "difqvsxmym",
		"DIF cluster signal vs xmod ymod;x track mod 100 [#mum];y track mod 100 [#mum];DIF <cluster signal> [ToT]",
		50, 0, 100, 50, 0, 100, 0, 200 );
  TProfile difqvsxm( "difqvsxm",
		     "DIF cluster signal vs xmod;x track mod 100 [#mum];DIF <cluster signal> [ToT]",
		     50, 0, 100, 0, 100 );
  TProfile difqvsxm5( "difqvsxm5",
		      "DIF cluster signal vs xmod;x track mod 50 [#mum];DIF <cluster signal> [ToT]",
		      50, 0, 50, 0, 100 );
  TProfile difqvsym( "difqvsym",
		     "DIF cluster signal vs ymod;y track mod 100 [#mum];DIF <cluster signal> [ToT]",
		     50, 0, 100, 0, 100 );

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
		     16, -0.5, 15.5 );
  TH1I dutpxqbeamHisto( "dutpxqbeam",
		     "DUT beam pixel signal;DUT pixel signal [ToT];linked beam pixels",
		     16, -0.5, 15.5 );
  TH1I dutpxq1Histo( "dutpxq1",
		     "DUT 1-col pixel signal;DUT pixel signal [ToT];linked 1-col pixels",
		     16, -0.5, 15.5 );
  TH1I dutpxq2Histo( "dutpxq2",
		     "DUT 2-col pixel signal;DUT pixel signal [ToT];linked 2-col pixels",
		     16, -0.5, 15.5 );

  TH1I dutpxqeveHisto( "dutpxqeve",
		       "DUT even pixel signal;DUT pixel signal [ToT];linked even pixels",
		       16, -0.5, 15.5 );
  TH1I dutpxqoddHisto( "dutpxqodd",
		       "DUT odd pixel signal;DUT pixel signal [ToT];linked odd pixels",
		       16, -0.5, 15.5 );

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
		     16, -0.5, 15.5, 16, -0.5, 15.5 );
  TH2I dutpxqq0Histo( "dutpxqq0",
		      "DUT even-odd row correlation;DUT pixel signal [ToT];DUT pixel signal [ToT];linked even-odd row pixels",
		      16, -0.5, 15.5, 16, -0.5, 15.5 );
  TH2I dutpxqq1Histo( "dutpxqq1",
		      "DUT odd-even row correlation;DUT pixel signal [ToT];DUT pixel signal [ToT];linked odd-even row pixels",
		      16, -0.5, 15.5, 16, -0.5, 15.5 );
  TH1I dutpxq1eveHisto( "dutpxq1eve",
			"DUT 1-col even row 1st pixel signal;DUT pixel signal [ToT];linked 1-col even row 1st pixels",
			16, -0.5, 15.5 );
  TH1I dutpxq1oddHisto( "dutpxq1odd",
			"DUT 1-col odd row 1st pixel signal;DUT pixel signal [ToT];linked 1-col odd row 1st pixels",
			16, -0.5, 15.5 );
  TH1I dutpxq2eveHisto( "dutpxq2eve",
			"DUT 1-col even row 2nd pixel signal;DUT pixel signal [ToT];linked 1-col even row 2nd pixels",
			16, -0.5, 15.5 );
  TH1I dutpxq2oddHisto( "dutpxq2odd",
			"DUT 1-col odd row 2nd pixel signal;DUT pixel signal [ToT];linked 1-col odd row 2nd pixels",
			16, -0.5, 15.5 );

  TH1I dutcolszHisto( "dutcolsz", "DUT column size;DUT column size [rows];columns inside linked clusters",
		      20, 0.5, 20.5 );
  TProfile dutcolszvsc( "dutcolszvsc", "DUT column size vs column;cluster column;DUT <column size> [rows]",
			40, 0.5, 40.5 );
  TH1I dutcolqHisto( "dutcolq", "DUT column signal;DUT column signal [ToT];columns inside linked clusters",
		     80, 0.5, 80.5 );
  TProfile dutcolqvsc( "dutcolqvsc", "DUT column signal vs column;cluster column;DUT <column signal> [ToT]",
		       40, 0.5, 40.5 );

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

  TH1I synqseedHisto( "synqseed",
		      "SYN seed row signal;SYN seed row signal [ToT];synked SYN clusters",
		      16, -0.5, 15.5 );
  TH1I synqpairHisto( "synqpair",
		      "SYN pair row signal;SYN pair row signal [ToT];synked SYN clusters",
		      16, -0.5, 15.5 );

  TProfile synqseedvsym( "synqseedvsym",
		      "SYN seed row signal vs ymod;y track mod 50 [#mum];1-col SYN <seed row signal> [ToT]",
		      50, 0, 50, -1, 99 );
  TProfile synqpairvsym( "synqpairvsym",
		      "SYN pair row signal vs ymod;y track mod 50 [#mum];1-col SYN <pair row signal> [ToT]",
		      50, 0, 50, -1, 99 );
  TH2I synqymHisto( "synqym",
		    "SYN cluster signal vs ymod;y track mod 50 [#mum];SYN <row signal> [ToT]",
		    25, 0, 50, 31, -15.5, 15.5 );

  TProfile synpqvsym( "synpqvsym",
		      "SYN signal covariance;y track mod 50 [#mum];SYN < ToT_{i-1} #upoint ToT_{i} >",
		      50, 0, 50, -1, 999 );
  TProfile synppvsym( "synppvsym",
		      "SYN signal variance;y track mod 50 [#mum];SYN < ToT_{i-1} #upoint ToT_{i-i} >",
		      50, 0, 50, -1, 999 );
  TProfile synqqvsym( "synqqvsym",
		      "SYN signal variance;y track mod 50 [#mum];SYN < ToT_{i} #upoint ToT_{i} >",
		      50, 0, 50, -1, 999 );

  TH1I synrowminHisto( "synrowmin",
			"SYN 1st pixel row;SYN 1st pixel row;synked 1-col clusters",
			nbr, 0, nbr );
  TH1I synrowmaxHisto( "synrowmax",
			"SYN lst pixel row;SYN lst pixel row;synked 1-col clusters",
			nbr, 0, nbr );
  TH1I synrowmin2Histo( "synrowmin2",
			"SYN 1st pixel row;SYN 1st pixel row;synked 1-col 2-row clusters",
			nbr, 0, nbr );
  TH1I synrowmax2Histo( "synrowmax2",
			"SYN lst pixel row;SYN lst pixel row;synked 1-col 2-row clusters",
			nbr, 0, nbr );

  TH1I linqseedHisto( "linqseed",
		      "LIN seed row signal;LIN seed row signal [ToT];linked LIN clusters",
		      16, -0.5, 15.5 );
  TH1I linqpairHisto( "linqpair",
		      "LIN pair row signal;LIN pair row signal [ToT];linked LIN clusters",
		      16, -0.5, 15.5 );

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
  TH1I linrowmin2Histo( "linrowmin2",
			"LIN 1st pixel row;LIN 1st pixel row;linked 1-col 2-row clusters",
			nbr, 0, nbr );
  TH1I linrowmax2Histo( "linrowmax2",
			"LIN lst pixel row;LIN lst pixel row;linked 1-col 2-row clusters",
			nbr, 0, nbr );

  TH1I difqseedHisto( "difqseed",
		      "DIF seed row signal;DIF seed row signal [ToT];linked DIF clusters",
		      16, -0.5, 15.5 );
  TH1I difqpairHisto( "difqpair",
		      "DIF pair row signal;DIF pair row signal [ToT];linked DIF clusters",
		      16, -0.5, 15.5 );

  TProfile difqseedvsym( "difqseedvsym",
			 "DIF seed row signal vs ymod;y track mod 50 [#mum];1-col DIF <seed row signal> [ToT]",
			 50, 0, 50, -1, 99 );
  TProfile difqpairvsym( "difqpairvsym",
			 "DIF pair row signal vs ymod;y track mod 50 [#mum];1-col DIF <pair row signal> [ToT]",
			 50, 0, 50, -1, 99 );
  TH2I difqymHisto( "difqym",
		    "DIF cluster signal vs ymod;y track mod 50 [#mum];DIF <row signal> [ToT]",
		    25, 0, 50, 31, -15.5, 15.5 );

  TProfile difpqvsym( "difpqvsym",
		      "DIF signal covariance;y track mod 50 [#mum];DIF < ToT_{i-1} #upoint ToT_{i} >",
		      50, 0, 50, -1, 999 );
  TProfile difppvsym( "difppvsym",
		      "DIF signal variance;y track mod 50 [#mum];DIF < ToT_{i-1} #upoint ToT_{i-i} >",
		      50, 0, 50, -1, 999 );
  TProfile difqqvsym( "difqqvsym",
		      "DIF signal variance;y track mod 50 [#mum];DIF < ToT_{i} #upoint ToT_{i} >",
		      50, 0, 50, -1, 999 );

  TH1I difrowminHisto( "difrowmin",
			"DIF seed pixel row;DIF seed pixel row;linked 1-col clusters",
			nbr, 0, nbr );
  TH1I difrowmaxHisto( "difrowmax",
			"DIF lst pixel row;DIF lst pixel row;linked 1-col clusters",
			nbr, 0, nbr );
  TH1I difrowmin2Histo( "difrowmin2",
			"DIF seed pixel row;DIF seed pixel row;linked 1-col 2-row clusters",
			nbr, 0, nbr );
  TH1I difrowmax2Histo( "difrowmax2",
			"DIF lst pixel row;DIF lst pixel row;linked 1-col 2-row clusters",
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

  TH1I dutpxqbcHisto[32];
  for( int frm = 0; frm < 32; ++frm )
    dutpxqbcHisto[frm] =
      TH1I( Form( "dutpxqbc%i", frm ),
	    Form( "DUT pixel signal in BC %i;DUT pixel signal in BC %i [ToT];linked pixels", frm, frm ),
	    16, -0.5, 15.5 );

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
		   2*nbx, -11, 11, -1, 2 ); // bin = col
  TProfile effvsx4( "effvsx4",
		   "DUT efficiency vs x;x track at DUT [mm];efficiency",
		   4*nbx, -11, 11, -1, 2 ); // bin = col
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
		     "DUT efficiency vs events;events;efficiency",
		     300, 0, 30*1000, -1, 2 );
  TProfile effvsev2( "effvsev2",
		     "DUT efficiency vs events;events;efficiency",
		     160, 0, 160*1000, -1, 2 );

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
      cout << "scope53 processing  " << run << "." << iev << "  taken " << evsec << endl;
    else if( iev < 100 && iev%10 == 0 )
      cout << "scope53 processing  " << run << "." << iev << "  taken " << evsec << endl;
    else if( iev < 1000 && iev%100 == 0 )
      cout << "scope53 processing  " << run << "." << iev << "  taken " << evsec << endl;
    else if( iev%1000 == 0 )
      cout << "scope53 processing  " << run << "." << iev << "  taken " << evsec << endl;

    StandardEvent sevt = eudaq::PluginManager::ConvertToStandard(evt);

    if( ldbg ) cout << "planes " << sevt.NumPlanes() << endl;

    vector < cluster > cl[9];

    for( size_t iplane = 0; iplane < sevt.NumPlanes(); ++iplane ) {

      const eudaq::StandardPlane &plane = sevt.GetPlane(iplane);

      if( ldbg )
	std::cout
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

      if( ipl < 0 || ipl > 6 ) {
	cout << "event " << iev << " wrong plane number " << ipl << endl;
	continue;
      }

      hpivot[ipl].Fill( plane.PivotPixel() );
      hnpx[ipl].Fill( plane.HitPixels() );
      hnframes[ipl].Fill( plane.NumFrames() ); // 32

      vector <pixel> pb; // for clustering

      // loop over frames, then pixels per frame

      for( unsigned frm = 0; frm < plane.NumFrames(); ++frm ) 

	for( size_t ipix = 0; ipix < plane.HitPixels( frm ); ++ipix ) {

	  if( ldbg ) 
	    std::cout << ": " << plane.GetX(ipix,frm)
		      << "." << plane.GetY(ipix,frm)
		      << "." << plane.GetPixel(ipix,frm) << " ";

	  int ix = plane.GetX(ipix,frm); // column
	  int iy = plane.GetY(ipix,frm); // row
	  int tot = plane.GetPixel(ipix,frm); // ToT 0..15

	  if( ipl == iDUT ) {
	    dutpxbcHisto.Fill( frm ); // before hot pixel masking
	    hmap[8]->Fill( ix, iy ); // before masking
	  }

	  // skip hot pixels:

	  int ipx = ix*ny[ipl] + iy;
	  if( hotset[ipl].count(ipx) ) continue;

	  pixel px;
	  px.col = ix; // ROC col
	  px.row = iy; // row
	  px.tot = tot;
	  px.frm = frm;
	  px.pivot = plane.GetPivot(ipix,frm);

	  if( ipl == iDUT ) {

	    dutpxcol0Histo.Fill( ix+0.5 );

	    px.tot += 1; // shift from zero

	    if( ix < 128 )
	      synpxqHisto.Fill( px.tot );
	    else if( ix < 264 )
	      linpxqHisto.Fill( px.tot );
	    else
	      difpxqHisto.Fill( px.tot );

	    dutpxqvsx.Fill( ix, px.tot );
	    dutpxqvsxy->Fill( ix, iy, px.tot );

	    if( ix < 128 )
	      synpxbcHisto.Fill( frm );
	    else if( ix < 264 )
	      linpxbcHisto.Fill( frm );
	    else
	      difpxbcHisto.Fill( frm );
	    dutpxbcvsx.Fill( ix, frm );
	    dutpxbcvsxy->Fill( ix, iy, frm );
	    dutpxbcvst5.Fill( evsec, frm ); // long run 34135: not stable

	    int thr = 0;
	    if( chip0 == -550 ) {
	      if( ix < 264 ) // lin
		//thr = 2; // hard enough: linnrowvsym
		thr = 3; // uniform 1st cluster row
	      else
		//thr = 2; // not enough
		//thr = 5; // better
		thr = 7; // uniform 1st cluster row
	    }
	    if( run == 35088 )
	      //thr = 1; // 35088 dxc 26
	      thr = 2; // 35088 dxc 27
	    //thr = 3; // 35088 dxc 30

	    if( px.tot < thr ) continue; // offline threshold

	    dutpxcol9Histo.Fill( ix+0.5 );

	    if( !fifty ) { // 100x25 from ROC to sensor:
	      px.col = ix/2; // 100 um
	      if( ix%2 ) 
		px.row = 2*iy + 0; // different from R4S
	      else
		px.row = 2*iy + 1; // see ed53 for shallow angle
	    }

	  } // DUT

	  pb.push_back(px);

	  hcol[ipl].Fill( ix ); // ROC
	  hrow[ipl].Fill( iy );
	  hmap[ipl]->Fill( ix, iy ); // after thr

	} // pix

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

    } // planes

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
	  tritxHisto.Fill( slpx );
	  trityHisto.Fill( slpy );

	} // cl B

      } // cl C

    } // cl A

    ntriHisto.Fill( triplets.size() );

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // triplets vs MOD and DUT:

    int nmdm = 0;
    int ntrimod = 0;

    double xcutMOD = 0.15;
    double ycutMOD = 0.15;

    int nmtd = 0;
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

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // triplets vs MOD clusters:

      bool ltrimod = 0;

      for( vector<cluster>::iterator c = cl[iMOD].begin(); c != cl[iMOD].end(); ++c ) {

	double ccol = c->col;
	double crow = c->row;
	double modx = ( ccol + 0.5 - nx[iMOD]/2 ) * ptchx[iMOD]; // -33..33 mm
	double mody = ( crow + 0.5 - ny[iMOD]/2 ) * ptchy[iMOD]; // -8..8 mm
	double q = c->signal;
	double q0 = q*normm;

	modxHisto.Fill( modx );
	modyHisto.Fill( mody );

	// residuals for pre-alignment:

	modsxaHisto.Fill( modx + x3m ); // peak
	moddxaHisto.Fill( modx - x3m ); // 

	modsyaHisto.Fill( mody + y3m ); // 
	moddyaHisto.Fill( mody - y3m ); // peak

	double moddx = modx - x4m;
	double moddy = mody - y4m;

	moddxHisto.Fill( moddx );
	moddyHisto.Fill( moddy );

	if( fabs( moddx ) < xcutMOD ) {

	  moddycHisto.Fill( moddy );
	  moddyvsx.Fill( -x3m, moddy ); // for rot
	  moddyvsy.Fill( y2m, moddy ); // for tilt
	  moddyvsty.Fill( syA, moddy );
	  moddyvst5.Fill( evsec, moddy );

	}

	if( fabs( moddy ) < ycutMOD ) {

	  moddxcHisto.Fill( moddx );
	  moddxvsx.Fill( -x1m, moddx ); // for turn
	  moddxvsy.Fill( y3m, moddx ); // for rot
	  moddxvstx.Fill( sxA, moddx );
	  moddxvst5.Fill( evsec, moddx );

	}

	if( fabs( moddx ) < xcutMOD &&
	    fabs( moddy ) < ycutMOD ) {

	  modnpxHisto.Fill( c->size );
	  modqHisto.Fill( q );
	  modq0Histo.Fill( q0 );
	  modnpxvsxmym->Fill( xmodm*1E3, ymodm*1E3, c->size );

	}

	if( fabs( moddx ) < xcutMOD &&
	    fabs( moddy ) < ycutMOD ) {

	  modlkxBHisto.Fill( xB );
	  modlkyBHisto.Fill( yB );
	  modlkxHisto.Fill( x4m );
	  modlkyHisto.Fill( y4m );
	  modlkcolHisto.Fill( ccol );
	  modlkrowHisto.Fill( crow );

	  //if( crow > 80 ) cout << "B link " << crow << " ev " << iev << endl;
	  //if( crow < 80 ) cout << "A link " << crow << " ev " << iev << endl;
	  //if( iev > 5100200 && crow < 80 ) cout << "link " << crow << " ev " << iev << endl; // A

	  triplets[iA].lk = 1;
	  nmdm = 1; // we have a MOD-triplet match in this event
	  ++ntrimod;
	  ltrimod = 1;

	  for( vector<pixel>::iterator px = c->vpix.begin(); px != c->vpix.end(); ++px ) {
	    tripxfrmlkHisto.Fill( px->frm ); // 0 or 1
	    tripxpivlkHisto.Fill( px->pivot ); // 0 or 1, cannot reduce pile-up
	  }

	} // MOD link x and y

      } // MOD

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

      double x4 = x3 + DUTalignx; // shift to zero
      double y4 = y3 + DUTaligny; // invert y, shift to zero

      double xmod = fmod( 9.000 + x4, 0.100 ); // [0,0.100] mm
      double ymod = fmod( 9.000 + y4, 0.100 ); // [0,0.100] mm
      double xmod2 = fmod( 9.000 + x4, 0.200 ); // [0,0.200] mm
      double ymod2 = fmod( 9.000 + y4, 0.200 ); // [0,0.200] mm
      double xmod5 = fmod( 9.000 + x4, 0.050 ); // [0,0.050] mm
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

      dutxyHisto->Fill( x4, y4 );

      int sect; // Sync, Lin, Diff

      if( x4 < -3.6 ) // Sync 0..127
	sect = 0;

      else if( x4 < 3.2 ) // Lin 128..263
	sect = 1;

      else // Diff 264..399
	sect = 2;

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

	//if( ccol > 200 && ccol < 201 ) cout << "col " << ccol << "  " << dutx << endl;
	//if( crow > 95.5 && crow < 96.5 ) cout << "row " << crow << "  " << duty << endl;

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

	if( fabs(dutdx) < fabs(dxmin) ) dxmin = dutdx;
	if( fabs(dutdy) < fabs(dymin) ) dymin = dutdy;

	if( fabs( dutdy ) < ycutDUT ) {

	  dutxxHisto->Fill( x4, dutx );

	  dutdxcHisto.Fill( dutdx );
	  if( sect == 1 )
	    lindxcHisto.Fill( dutdx );
	  else if( sect == 2 )
	    difdxcHisto.Fill( dutdx );

	  if(      c->nrow == 1 )
	    dutdxc1Histo.Fill( dutdx );
	  else if( c->nrow == 2 )
	    dutdxc2Histo.Fill( dutdx );

	  dutdxvsx.Fill( x4, dutdx ); // for turn
	  dutdxvsy.Fill( y4, dutdx ); // for rot
	  dutdxvstx.Fill( sxA, dutdx );
	  if( x4 < 0 )
	    dutdxvstx0.Fill( sxA, dutdx );
	  else
	    dutdxvstx1.Fill( sxA, dutdx );
	  dutdxvsxm.Fill( xmod*1E3, dutdx );
	  dutdxvst2.Fill( evsec, dutdx );
	  dutdxvst5.Fill( evsec, dutdx );

	  dutmadxvsx.Fill( x4, fabs(dutdx) );
	  dutmadxvsxm.Fill( xmod*1E3, fabs(dutdx) );
	  dutmadxvstx.Fill( sxA, fabs(dutdx) );
	  dutmadxvsq.Fill( Q0, fabs(dutdx) );
	  if( Q0 > qL && Q0 < qR )
	    dutdxcqHisto.Fill( dutdx );

	} // cut y

	if( fabs( dutdx ) < xcutDUT ) {

	  dutyyHisto->Fill( y4, duty );
	  dutdy0Histo.Fill( dutdy );
	  dutdycHisto.Fill( dutdy );

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

	// cut x and y:

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

	  if( sect == 0 ) {

	    synqHisto.Fill( Q );
	    synq0Histo.Fill( Q0 );
	    synqvsx.Fill( x4, Q0 );
	    synqvsy.Fill( y4, Q0 );
	    synqvsxy->Fill( x4, y4, Q0 );
	    synqvsr.Fill( drbeam, Q0 );
	    if( drbeam < rbeam )
	      synqbeamHisto.Fill( Q );

	    synnpxHisto.Fill( npx );
	    synncolHisto.Fill( c->ncol );
	    synnrowHisto.Fill( c->nrow );

	    if( c->ncol == 1 ) {
	      synnrow1Histo.Fill( c->nrow );
	      vector<pixel>::iterator px = c->vpix.begin(); // 1st pixel
	      if( px->row%2 )
		synnrow1oddHisto.Fill( c->nrow );
	      else
		synnrow1eveHisto.Fill( c->nrow );
	    }

	    if( c->nfrm == 1 ) {
	      synnpxfHisto.Fill( npx );
	      synncolfHisto.Fill( c->ncol );
	      synnrowfHisto.Fill( c->nrow );
	    }
	    else {
	      synnpxffHisto.Fill( npx );
	      synncolffHisto.Fill( c->ncol );
	      synnrowffHisto.Fill( c->nrow );
	    }

	    synnpxvsxmym->Fill( xmod*1E3, ymod*1E3, npx );
	    synncolvsxm.Fill( xmod*1E3, c->ncol );
	    synncolvsym.Fill( ymod*1E3, c->ncol );

	    synnrowvsxmym->Fill( xmod*1E3, ymod*1E3, c->nrow );
	    synnrowvsxm.Fill( xmod*1E3, c->nrow );
	    synnrowvsxm5.Fill( xmod5*1E3, c->nrow );
	    if( !fifty && !rot90 && xmod > 0.010 && xmod < 0.090 )
	      synnrowvsym.Fill( ymod*1E3, c->nrow );

	    synqvsxmym->Fill( xmod*1E3, ymod*1E3, Q );
	    synqvsxm.Fill( xmod*1E3, Q );
	    synqvsxm5.Fill( xmod5*1E3, Q );
	    synqvsym.Fill( ymod*1E3, Q );

	  } // Syn sect

	  if( sect == 1 ) {

	    linqHisto.Fill( Q );
	    linq0Histo.Fill( Q0 );
	    linqvsx.Fill( x4, Q0 );
	    linqvsy.Fill( y4, Q0 );
	    linqvsxy->Fill( x4, y4, Q0 );
	    linqvsr.Fill( drbeam, Q0 );

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

	    if( c->nfrm == 1 ) {
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
	    linncolvsxm.Fill( xmod*1E3, c->ncol );
	    linncolvsym.Fill( ymod*1E3, c->ncol );

	    linnrowvsxmym->Fill( xmod*1E3, ymod*1E3, c->nrow );
	    linnrowvsxm.Fill( xmod*1E3, c->nrow );
	    linnrowvsxm5.Fill( xmod5*1E3, c->nrow );
	    if( !fifty && !rot90 && xmod > 0.010 && xmod < 0.090 )
	      linnrowvsym.Fill( ymod*1E3, c->nrow );

	    linqvsxmym->Fill( xmod*1E3, ymod*1E3, Q );
	    linqvsxm.Fill( xmod*1E3, Q );
	    linqvsxm5.Fill( xmod5*1E3, Q );
	    linqvsym.Fill( ymod*1E3, Q );

	  } // Lin sect

	  else if( sect == 2 ) {

	    difqHisto.Fill( Q );
	    difq0Histo.Fill( Q0 );

	    difnpxHisto.Fill( npx );
	    difncolHisto.Fill( c->ncol );
	    difnrowHisto.Fill( c->nrow );

	    if( c->ncol == 1 ) {
	      difnrow1Histo.Fill( c->nrow );
	      vector<pixel>::iterator px = c->vpix.begin(); // 1st pixel
	      if( px->row%2 )
		difnrow1oddHisto.Fill( c->nrow );
	      else
		difnrow1eveHisto.Fill( c->nrow ); // twice more, similar shape
	    }

	    if( c->nfrm == 1 ) {
	      difnpxfHisto.Fill( npx );
	      difncolfHisto.Fill( c->ncol );
	      difnrowfHisto.Fill( c->nrow );
	    }
	    else {
	      difnpxffHisto.Fill( npx );
	      difncolffHisto.Fill( c->ncol );
	      difnrowffHisto.Fill( c->nrow );
	    }

	    difnpxvsxmym->Fill( xmod*1E3, ymod*1E3, npx );
	    difncolvsxm.Fill( xmod*1E3, c->ncol );
	    difncolvsym.Fill( ymod*1E3, c->ncol );

	    difnrowvsxmym->Fill( xmod*1E3, ymod*1E3, c->nrow );
	    difnrowvsxm.Fill( xmod*1E3, c->nrow );
	    if( !fifty && !rot90 && xmod > 0.010 && xmod < 0.090 )
	      difnrowvsym.Fill( ymod*1E3, c->nrow );

	    difqvsxmym->Fill( xmod*1E3, ymod*1E3, Q );
	    difqvsxm.Fill( xmod*1E3, Q );
	    difqvsxm5.Fill( xmod5*1E3, Q );
	    difqvsym.Fill( ymod*1E3, Q );

	  } // Diff-sect

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

	  for( vector<pixel>::iterator px = c->vpix.begin(); px != c->vpix.end(); ++px ) {

	    int icol = px->col; // sensor
	    int irow = px->row;
	    int itot = px->tot;
	    int ifrm = px->frm;

	    dutpxcolHisto.Fill( icol+0.5 ); // Dif: 8-fold pattern
	    dutpxrowHisto.Fill( irow+0.5 );

	    if( icol < colmin ) colmin = icol;
	    if( icol > colmax ) colmax = icol;
	    if( irow < rowmin ) rowmin = irow;
	    if( irow > rowmax ) rowmax = irow;

	    if( px == c->vpix.begin() ) { // 1st px in readout order

	      dutpxrow1Histo.Fill( irow+0.5 ); // even-odd pattern in 25x100
	      dutpxcol1Histo.Fill( icol+0.5 ); // 8-pattern in 50x50

	      if( c->ncol == 1 ) {

		dutpxrow11Histo.Fill( irow+0.5 ); // 70% even

		if( c->nfrm == 1 ) { // 1-frame clusters: pixels read out sequentially

		  dutpxrow111Histo.Fill( irow+0.5 ); // 75% even

		  if( c->nrow == 2 )
		    dutpxrow1112Histo.Fill( irow+0.5 ); // 100% even as it should be for 25x100

		  if( irow%2 )
		    dutpxcol111oddHisto.Fill( icol+0.5 );
		  else
		    dutpxcol111eveHisto.Fill( icol+0.5 ); // same shape

		} // frm

	      } // ncol

	      if( c->nrow == 1 )
		dutpxcol11Histo.Fill( icol+0.5 );

	    } // 1st px

	    ++colsz[icol];
	    colq[icol] += itot;
	    rowq[irow] += itot;

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
	    dutpxqbcHisto[ifrm].Fill( itot ); // charge per BC
	    if( c->size == 1 )
	      dutpxbclk1pHisto.Fill( ifrm );
	    if( c->ncol == 1 )
	      dutpxbclk1cHisto.Fill( ifrm );
	    if( c->nrow == 1 )
	      dutpxbclk1rHisto.Fill( ifrm );

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
	    dutrowpqmap->Fill( kcol+0.5, krow+0.5, rowq.at(krow-1) * rowq.at(krow) );
	    dutrowppmap->Fill( kcol+0.5, krow+0.5, rowq.at(krow-1) * rowq.at(krow-1) );
	    dutrowqqmap->Fill( kcol+0.5, krow+0.5, rowq.at(krow) * rowq.at(krow) );
	  }
	  if( kcol > 0 ) { // p = k-1, q = k
	    dutcolpqmap->Fill( kcol+0.5, krow+0.5, colq.at(kcol-1) * colq.at(kcol) );
	    dutcolppmap->Fill( kcol+0.5, krow+0.5, colq.at(kcol-1) * colq.at(kcol-1) );
	    dutcolqqmap->Fill( kcol+0.5, krow+0.5, colq.at(kcol) * colq.at(kcol) );
	  }

	  if( c->ncol == 1 ) { // for 100x25

	    if( sect == 0 ) { // Syn

	      // seed row and pair:

	      synqseedHisto.Fill( rowq.at(krow) );
	      synqpairHisto.Fill( rowq.at(lrow) );

	      synqseedvsym.Fill( ymod5*1E3, rowq.at(krow) );
	      synqpairvsym.Fill( ymod5*1E3, rowq.at(lrow) );
	      synqymHisto.Fill( ymod5*1E3, -rowq.at(lrow) );
	      synqymHisto.Fill( ymod5*1E3,  rowq.at(krow) );

	      // correlation profile: p = l, q = k

	      synpqvsym.Fill( ymod5*1E3, rowq.at(lrow) * rowq.at(krow) );
	      synppvsym.Fill( ymod5*1E3, rowq.at(lrow) * rowq.at(lrow) );
	      synqqvsym.Fill( ymod5*1E3, rowq.at(krow) * rowq.at(krow) );

	      synrowminHisto.Fill( rowmin+0.5 ); // 35150: strong even-odd
	      synrowmaxHisto.Fill( rowmax+0.5 );

	      if( c->nrow == 2 ) {
		synrowmin2Histo.Fill( rowmin+0.5 ); // mostly even
		synrowmax2Histo.Fill( rowmax+0.5 );
	      }

	    } // Syn

	    else if( sect == 1 ) { // Lin

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

	      linrowminHisto.Fill( rowmin+0.5 );
	      linrowmaxHisto.Fill( rowmax+0.5 );

	      if( c->nrow == 2 ) {
		linrowmin2Histo.Fill( rowmin+0.5 ); // mostly even
		linrowmax2Histo.Fill( rowmax+0.5 );
	      }

	    } // Lin

	    else if( sect == 2 ) { // Diff

	      // seed row:

	      difqseedHisto.Fill( rowq.at(krow) );
	      difqpairHisto.Fill( rowq.at(lrow) );

	      difqseedvsym.Fill( ymod5*1E3, rowq.at(krow) );
	      difqpairvsym.Fill( ymod5*1E3, rowq.at(lrow) );
	      difqymHisto.Fill( ymod5*1E3, -rowq.at(lrow) );
	      difqymHisto.Fill( ymod5*1E3,  rowq.at(krow) );

	      // correlation profile: p = k-1, q = k

	      difpqvsym.Fill( ymod5*1E3, rowq.at(lrow) * rowq.at(krow) );
	      difppvsym.Fill( ymod5*1E3, rowq.at(lrow) * rowq.at(lrow) );
	      difqqvsym.Fill( ymod5*1E3, rowq.at(krow) * rowq.at(krow) );

	      difrowminHisto.Fill( rowmin+0.5 );
	      difrowmaxHisto.Fill( rowmax+0.5 );

	      if( c->nrow == 2 ) {
		difrowmin2Histo.Fill( rowmin+0.5 );
		difrowmax2Histo.Fill( rowmax+0.5 );
	      }

	    } // Diff

	  } // 1-col

	  triplets[iA].lk = 1;
	  nmtd = 1; // we have a DUT-triplet match in this event
	  ++ntrilk;
	  ++nlk;

	} // DUT link x and y

	// for eff: nearest pixel

	for( unsigned ipx = 0; ipx < c->vpix.size(); ++ipx ) {

	  double px = ( c->vpix[ipx].col + 0.5 - nx[iDUT]/2 ) * ptchx[iDUT]; // mm
	  double py = ( c->vpix[ipx].row + 0.5 - ny[iDUT]/2 ) * ptchy[iDUT]; // mm
	  if( rot90 ) {
	    px = ( c->vpix[ipx].row + 0.5 - ny[iDUT]/2 ) * ptchy[iDUT]; // mm
	    py =-( c->vpix[ipx].col + 0.5 - nx[iDUT]/2 ) * ptchx[iDUT]; // mm
	  }
	  dutpxxyHisto->Fill( px, py );
	  double pdx = px - x4; // triplet extrapol
	  double pdy = py - y4;
	  double pdxy = sqrt( pdx*pdx + pdy*pdy );
	  if( pdxy < pdmin ) pdmin = pdxy;

	} // pix

      } // loop DUT clusters

      int nm[99] = {0};
      for( int iw = 1; iw < 99; ++iw )
	if( pdmin < iw*0.010 ) // 10 um bins
	  nm[iw] = 1; // eff

      if( nm[49] )
	dutxylkHisto->Fill( x4, y4 ); // tracks with cluster link

      if( nm[20]
	  && y4 < 4.7 && y4 > -4.7
	  && x4 <  8  && x4 > -1.5 // fiducial
	  ) {
	dutnlkHisto.Fill( nlk );
	dutdxminHisto.Fill( dxmin );
	dutdyminHisto.Fill( dymin );
      }
      if( shrt )
	dutnlkshrHisto.Fill( nlk ); // nlk = 2 => split clusters

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // DUT efficiency vs isolated MOD-linked fiducial tracks:

      if( ltrimod ) {

	bool fidx = 1;
	bool fidy = 1;
	if( chip0 == 501 ) { // rot90 Lin
	  if( x4 >  4.7 ) fidx = 0;
	  if( x4 < -4.7 ) fidx = 0;
	  if( y4 >  3.5 ) fidy = 0;
	  if( y4 < -2.1 ) fidy = 0; // packman cutout
	}
	if( chip0 == 504 ) { // straight Lin
	  if( x4 >  3.2 ) fidx = 0;
	  if( x4 < -3.5 ) fidx = 0;
	  if( y4 >  4.7 ) fidy = 0;
	  if( y4 < -4.7 ) fidy = 0;
	}
	if( chip0 == 509 ) { // straight Lin
	  if( x4 >  3.2 ) fidx = 0;
	  if( x4 < -3.5 ) fidx = 0;
	  if( y4 >  4.7 ) fidy = 0;
	  if( y4 < -4.7 ) fidy = 0;
	}
	if( chip0 == 509 && run >= 35150 && run <= 99999 ) { // straight Syn
	  fidx = 1;
	  fidy = 1;
	  if( x4 < -9.9 ) fidx = 0;
	  if( x4 > -3.7 ) fidx = 0;
	  if( y4 >  4.7 ) fidy = 0;
	  if( y4 < -4.7 ) fidy = 0;
	}
	if( chip0 == 524 && run < 34090 ) { // straight Lin
	  if( x4 >  3.2 ) fidx = 0;
	  if( x4 < -3.5 ) fidx = 0;
	  if( y4 >  4.7 ) fidy = 0;
	  if( y4 < -4.7 ) fidy = 0;
	}
	if( chip0 == 524 ) { // straight Lin+Diff
	  if( x4 >  9.9 ) fidx = 0;
	  if( x4 < -3.5 ) fidx = 0;
	  if( y4 >  4.7 ) fidy = 0;
	  if( y4 < -4.7 ) fidy = 0;
	}
	if( chip0 == 531 ) { // straight Lin+Diff
	  if( x4 >  9.9 ) fidx = 0;
	  if( x4 < -3.5 ) fidx = 0;
	  if( y4 >  4.7 ) fidy = 0;
	  if( y4 < -4.7 ) fidy = 0;
	}
	if( chip0 == 520 ) { // straight Lin+Diff
	  if( x4 >  9.9 ) fidx = 0;
	  if( x4 < -3.5 ) fidx = 0;
	  if( y4 >  4.7 ) fidy = 0;
	  if( y4 < -4.7 ) fidy = 0;
	}
	if( chip0 == 550 ) { // straight Lin+Diff
	  if( x4 >  9.9 ) fidx = 0;
	  if( x4 < -3.5 ) fidx = 0;
	  if( y4 >  4.7 ) fidy = 0;
	  if( y4 < -4.7 ) fidy = 0;
	}
	if( chip0 == 1470 ) { // straight Lin
	  if( x4 >  3.2 ) fidx = 0;
	  if( x4 < -3.5 ) fidx = 0;
	  if( y4 >  4.7 ) fidy = 0;
	  if( y4 < -4.7 ) fidy = 0;
	}

	// from track x, y (at DUT) to sensor col, row:

	int kcol = ( x4 / ptchx[iDUT] + 0.5*nx[iDUT] ); // straight 50x50
	int krow = ( y4 / ptchy[iDUT] + 0.5*ny[iDUT] );
	if( rot90 ) { // 25x100
	  kcol = (-y4 / ptchx[iDUT] + 0.5*nx[iDUT] ); // 0..398 even
	  krow = ( x4 / ptchy[iDUT] + 0.5*ny[iDUT] ); // 0..191 full
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

	    // efficiency around dead pixels:

	    for( int idead : deadset ) { // C++11 range-based for loop over all dead px
	      int dcol = idead/384;
	      int drow = idead%384; // sensor-based
	      double px = ( dcol + 0.5 - nx[iDUT]/2 ) * ptchx[iDUT]; // mm
	      double py = ( drow + 0.5 - ny[iDUT]/2 ) * ptchy[iDUT]; // mm
	      double pdx = px - x4; // triplet extrapol
	      double pdy = py - y4;
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

      } // sixlk

    } // loop triplets

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
  histoFile->Write();
  histoFile->Close();

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // MOD alignment:

  if( moddxaHisto.GetEntries() > 9999 ) {

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
	     << ": extra rot " << fdyvsx->GetParameter(1) << endl;
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
	 << "  alignx " << newMODalignx << endl
	 << "  aligny " << newMODaligny << endl
	 << "  rot    " << MODrot << endl
	 << "  tilt   " << MODtilt << endl
	 << "  turn   " << MODturn << endl
	 << "  dz     " << MODz - zz[1] << endl
      ;

  } // MOD

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // DUT alignment:

  cout << endl
       << "DUT efficiency " << 100*effvst5.GetMean(2) << "%"
       << " from " << effvst5.GetEntries() << " events"
       << endl;

  if( dutdxaHisto.GetEntries() > 999 ) {

    double newDUTalignx = DUTalignx0;
    double newDUTaligny = DUTaligny0;

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

    if( DUTaligniteration > 0 && fabs( newDUTalignx - DUTalignx0 ) < 0.1 ) {

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
	newDUTalignx = DUTalignx0 + fgp0->GetParameter(1);
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
	  fabs( DUTturn ) > 0.3 // [deg]
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
      cout << "alignx changed by " << newDUTalignx - DUTalignx0 << " mm"
	   << ", need more iteration"
	   << endl;

    // y:

    cout << endl;

    if( DUTaligniteration > 0 && fabs( newDUTaligny - DUTaligny0 ) < 0.1 ) {

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
	newDUTaligny = DUTaligny0 + fgp0->GetParameter(1);
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
	  fabs( DUTtilt ) > 0.3 // [deg]
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
      cout << "aligny changed by " << newDUTaligny - DUTaligny0 << " mm"
	   << ", need more iteration"
	   << endl;

    // write new DUT alignment:

    cout << endl
	 << "DUT alignment iteration " << DUTaligniteration << endl
	 << "  alignx " << newDUTalignx << endl
	 << "  aligny " << newDUTaligny << endl
	 << "  rot    " << DUTrot << endl
	 << "  tilt   " << DUTtilt << endl
	 << "  turn   " << DUTturn << endl
	 << "  dz     " << DUTz - zz[3] << endl
      ;

    cout << "update DUT alignment file? (y/n)" << endl;
    string ans{"n"};
    string YES{"y"};

    if( fabs(DUTturn) < 44 )
      cin >> ans;

    if( ans == YES ) {

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

      cout << " to " << DUTalignFileName.str() << endl;

    }
    else
      cout << "no" << endl;

  } // stat
  else
    cout << "not enough for alignment" << endl;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // done

  cout << endl << histoFile->GetName() << endl;

  cout << endl;

  return 0;
}
