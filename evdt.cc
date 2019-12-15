
// Daniel Pitzl, DESY, Apr 2016, Jun 2019
// telescope 2-D event display using ROOT

// needs runs.dat

// make evdt
// evdt 36663

#include "eudaq/FileReader.hh"
#include "eudaq/PluginManager.hh"

#include <TApplication.h>
#include <TSystem.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TGraph.h>
#include <TH2I.h>
#include <TF1.h>
#include <TLine.h>

#include <sstream> // stringstream
#include <fstream> // filestream
#include <set>
#include <cmath> // fabs
#include <unistd.h> // usleep
#include <sys/ioctl.h> // kbhit

using namespace std;
using namespace eudaq;

struct pixel {
  int col;
  int row;
  int adc;
  double q;
  int ord;
  bool big;
};

struct cluster {
  vector <pixel> vpix; // Armin Burgmeier: list
  int size;
  int ncol, nrow;
  double col, row;
  double charge;
  bool big;
};

struct triplet {
  double xm;
  double ym;
  double zm;
  double sx;
  double sy;
  bool lk;
  double ttdmin;
};

//------------------------------------------------------------------------------
bool kbhit()
{
  usleep(1000); // [us]
  int byteswaiting;
  ioctl( 0, FIONREAD, &byteswaiting );
  return byteswaiting > 0;
}

//------------------------------------------------------------------------------
vector <cluster> getClus( vector <pixel> pb, int fCluCut = 1 ) // 1 = no gap
{
  // returns clusters with local coordinates
  // decodePixels should have been called before to fill pixel buffer pb 
  // simple clusterization
  // cluster search radius fCluCut ( allows fCluCut-1 empty pixels)

  vector<cluster> vc;
  if( pb.size() == 0 ) return vc;

  vector <bool> gone( pb.size() ); // initialized to zero

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
        if( !gone[i] ) { // unused pixel
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
    c.big = 0;
    int minx = 999;
    int maxx = 0;
    int miny = 999;
    int maxy = 0;

    for( vector<pixel>::iterator p = c.vpix.begin();  p != c.vpix.end();  ++p ) {
      double Qpix = p->q; // calibrated [Vcal]

      c.charge += Qpix;
      sumQ += Qpix;
      c.col += (*p).col*Qpix;
      c.row += (*p).row*Qpix;
      if( p->big ) c.big = 1;
      if( p->col > maxx ) maxx = p->col;
      if( p->col < minx ) minx = p->col;
      if( p->row > maxy ) maxy = p->row;
      if( p->row < miny ) miny = p->row;
    }

    //cout << "(cluster with " << c.vpix.size() << " pixels)" << endl;

    if( ! ( c.charge == 0 ) ) {
      c.col /= sumQ;
      c.row /= sumQ;
    }
    else {
      c.col = (*c.vpix.begin()).col;
      c.row = (*c.vpix.begin()).row;
      cout << "GetClus: cluster with zero charge" << endl;
    }

    c.ncol = maxx-minx+1;
    c.nrow = maxy-miny+1;

    vc.push_back(c); // add cluster to vector

    // look for a new seed = used pixel:

    while( ( ++seed < pb.size() ) && gone[seed] );

  } // while over seeds

  // nothing left, return clusters

  return vc;

} // getClus

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

  int lev = 9; // last event displayed

  for( int i = 1; i < argc; ++i ) {

    if( !strcmp( argv[i], "-l" ) )
      lev = atoi( argv[++i] ); // last event

  } // argc

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // runs.dat:

  cout << endl;

  string geoFileName( "geo.dat" );
  double pbeam = 5.6;

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

      if( tag == GEO ) {
	tokenizer >> geoFileName;
	continue;
      }

      if( tag == GeV ) {
	tokenizer >> pbeam;
	continue;
      }

      // anything else on the line and in the file gets ignored

    } // while getline

    if( found )
      cout 
	<< "settings for run " << run << ":" << endl
	<< "  beam " << pbeam << " GeV" << endl
	<< "  geo file " << geoFileName << endl
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

      if( ipl < 0 || ipl >= 9 ) {
	cout << "wrong plane number " << ipl << endl;
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
  // alignment:

  int aligniteration = 0;
  double alignx[9];
  double aligny[9];
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
  // can there be instructions between if and else ?
  else {

    cout << "read alignment from " << alignFileName.str() << endl;

    string hash( "#" );
    string iteration( "iteration" );
    string plane( "plane" );
    string shiftx( "shiftx" );
    string shifty( "shifty" );
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

      if( ipl < 0 || ipl >= 9 ) {
	cout << "wrong plane number " << ipl << endl;
	continue;
      }

      double val;
      tokenizer >> val;
      if(      tag == shiftx )
	alignx[ipl] = val;
      else if( tag == shifty )
	aligny[ipl] = val;
      else if( tag == rotxvsy )
	rotx[ipl] = val;
      else if( tag == rotyvsx )
	roty[ipl] = val;

      // anything else on the line and in the file gets ignored

    } // while getline

  } // alignFile

  ialignFile.close();

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // hot pixels:

  ostringstream hotFileName; // output string stream

  hotFileName << "hot_" << run << ".dat";

  ifstream ihotFile( hotFileName.str() );

  set <int> hotset[9];

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
      cout << line << endl;

      if( line.empty() ) continue;

      stringstream tokenizer( line );
      string tag;
      tokenizer >> tag; // leading white space is suppressed
      if( tag.substr(0,1) == hash ) // comments start with #
	continue;

      if( tag == plane )
	tokenizer >> ipl;

      if( ipl < 0 || ipl >= 6 ) {
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
  // ROOT:

  cout << "ROOT application..." << endl;

  TApplication theApp( "comet", &argc, argv );

  gStyle->SetTextFont( 62 ); // 62 = Helvetica bold
  gStyle->SetTextAlign( 11 );

  gStyle->SetTickLength( -0.02, "x" ); // tick marks outside
  gStyle->SetTickLength( -0.01, "y" );
  gStyle->SetTickLength( -0.01, "z" );

  gStyle->SetLabelOffset( 0.022, "x" );
  gStyle->SetLabelOffset( 0.013, "y" );
  gStyle->SetLabelOffset( 0.022, "z" );

  gStyle->SetTitleOffset( 1.6, "x" );
  gStyle->SetTitleOffset( 1.6, "y" );
  gStyle->SetTitleOffset( 1.7, "z" );

  gStyle->SetLabelFont( 62, "X" );
  gStyle->SetLabelFont( 62, "Y" );
  gStyle->SetLabelFont( 62, "z" );

  gStyle->SetTitleFont( 62, "X" );
  gStyle->SetTitleFont( 62, "Y" );
  gStyle->SetTitleFont( 62, "z" );

  gStyle->SetTitleBorderSize( 0 ); // no frame around global title
  gStyle->SetTitleAlign( 13 ); // 13 = left top align
  gStyle->SetTitleX( 0.12 ); // global title
  gStyle->SetTitleY( 0.98 ); // global title

  gStyle->SetLineWidth( 1 ); // frames
  gStyle->SetHistLineColor( 4 ); // 4=blau
  gStyle->SetHistLineWidth( 3 );
  gStyle->SetHistFillColor( 5 ); // 5=gelb
  //  gStyle->SetHistFillStyle(4050); // 4050 = half transparent
  gStyle->SetHistFillStyle( 1001 ); // 1001 = solid

  gStyle->SetFrameLineWidth( 2 );

  // statistics box:

  gStyle->SetOptStat( 0 );
  gStyle->SetStatFormat( "8.6g" ); // more digits, default is 6.4g
  gStyle->SetStatFont( 42 ); // 42 = Helvetica normal
  //  gStyle->SetStatFont(62); // 62 = Helvetica bold
  gStyle->SetStatBorderSize( 1 ); // no 'shadow'

  gStyle->SetStatX( 0.80 );
  gStyle->SetStatY( 0.95 );

  gStyle->SetPalette( 1 ); // rainbow colors

  gStyle->SetHistMinimumZero(  ); // no zero suppression

  gStyle->SetOptDate( 0 );

  TApplication app( "app", 0, 0 );

  TCanvas * c1 = new TCanvas( "c1", "RD53A event display", 900, 500 );
  TCanvas * c2 = new TCanvas( "c2", "RD53A event display", 900, 500 );

  c1->SetBottomMargin( 0.08 );
  c1->SetLeftMargin( 0.03 );
  c1->SetRightMargin( 0.10 );

  c2->SetBottomMargin( 0.08 );
  c2->SetLeftMargin( 0.03 );
  c2->SetRightMargin( 0.10 );

  gPad->Update(  ); // required

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
  int nevd = 0;
  uint64_t evTLU0 = 0;
  const double fTLU = 384E6; // 384 MHz TLU clock

  bool more = 1; // event displays
  string Q{"q"};

  do {
    // Get next event:
    DetectorEvent evt = reader->GetDetectorEvent();

    if( evt.IsBORE() ) {
      eudaq::PluginManager::Initialize(evt);
      continue;
    }

    ++iev;

    bool ldbg = 0;

    if( iev == 1 )
      ldbg = 1;

    if( lev < 100 )
      ldbg = 1;

    uint64_t evTLU = evt.GetTimestamp(); // 384 MHz = 2.6 ns
    if( iev < 2  )
      evTLU0 = evTLU;
    double evsec = (evTLU - evTLU0) / fTLU;

    if( iev < 10 )
      cout << "scope processing  " << iev << "  taken " << evsec << endl;
    else if( iev < 100 && iev%10 == 0 )
      cout << "scope processing  " << iev << "  taken " << evsec << endl;
    else if( iev < 1000 && iev%100 == 0 )
      cout << "scope processing  " << iev << "  taken " << evsec << endl;
    else if( iev%1000 == 0 )
      cout << "scope processing  " << iev << "  taken " << evsec << endl;

    StandardEvent sevt = eudaq::PluginManager::ConvertToStandard(evt);

    vector <cluster> cl[9];
    int ncl = 0;

    for( size_t iplane = 0; iplane < sevt.NumPlanes(); ++iplane ) {

      const eudaq::StandardPlane &plane = sevt.GetPlane(iplane);

      std::vector<double> pxl = plane.GetPixels<double>();

      if( ldbg ) cout << "plane " << plane.ID() << ": ";

      // /home/pitzl/eudaq/main/include/eudaq/CMSPixelHelper.hh

      int ipl = plane.ID();

      vector <pixel> pb; // for clustering

      for( size_t ipix = 0; ipix < pxl.size(); ++ipix ) {

	if( ldbg ) 
	  cout << plane.GetX(ipix)
	       << " " << plane.GetY(ipix)
	       << " " << plane.GetPixel(ipix) << " ";

	int ix = plane.GetX(ipix); // global column 0..415
	int iy = plane.GetY(ipix); // global row 0..159
	int adc = plane.GetPixel(ipix); // ADC 0..255

	// skip hot pixels:

	int ipx = ix*ny[ipl] + iy;
	if( hotset[ipl].count(ipx) ) continue;

	double q = adc;

	// fill pixel block for clustering:

	pixel px;
	px.col = ix; // col
	px.row = iy; // row
	px.adc = adc;
	px.q = q;
	px.ord = pb.size(); // readout order
	px.big = 0;
	pb.push_back(px);

      } // pix

      if( ldbg ) cout << endl;

      // clustering:

      cl[ipl] = getClus( pb );

      if( ldbg ) cout << "  z " << int(zz[ipl]+0.5)
		      << ", clusters " << cl[ipl].size()
		      << endl;

      ncl += cl[ipl].size();

    } // planes

    // event selection:

    if( ncl < 6 ) continue;

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // make driplets 3+5-4:

    vector <triplet> driplets;

    double driCut = 0.05; // [mm]

    double dz46 = zz[6] - zz[4];

    for( vector<cluster>::iterator cA = cl[4].begin(); cA != cl[4].end(); ++cA ) {

      double xA = cA->col*ptchx[4] - alignx[4];
      double yA = cA->row*ptchy[4] - aligny[4];
      double xmid = xA - midx[4];
      double ymid = yA - midy[4];
      xA = xmid - ymid*rotx[4];
      yA = ymid + xmid*roty[4];

      for( vector<cluster>::iterator cC = cl[6].begin(); cC != cl[6].end(); ++cC ) {

	double xC = cC->col*ptchx[6] - alignx[6];
	double yC = cC->row*ptchy[6] - aligny[6];
	double xmid = xC - midx[6];
	double ymid = yC - midy[6];
	xC = xmid - ymid*rotx[6];
	yC = ymid + xmid*roty[6];

	double avx = 0.5 * ( xA + xC ); // mid
	double avy = 0.5 * ( yA + yC );
	double avz = 0.5 * ( zz[4] + zz[6] ); // mid z
 
	double slpx = ( xC - xA ) / dz46; // slope x
	double slpy = ( yC - yA ) / dz46; // slope y

	// middle plane B = 5:

	for( vector<cluster>::iterator cB = cl[5].begin(); cB != cl[5].end(); ++cB ) {

	  double xB = cB->col*ptchx[5] - alignx[5];
	  double yB = cB->row*ptchy[5] - aligny[5];
	  double xmid = xB - midx[5];
	  double ymid = yB - midy[5];
	  xB = xmid - ymid*rotx[5];
	  yB = ymid + xmid*roty[5];

	  // interpolate track to B:

	  double dz = zz[5] - avz;
	  double xk = avx + slpx * dz; // driplet at k
	  double yk = avy + slpy * dz;

	  double dx3 = xB - xk;
	  double dy3 = yB - yk;

	  // telescope driplet cuts:

	  if( fabs(dx3) > driCut ) continue;
	  if( fabs(dy3) > driCut ) continue;

	  triplet dri;
	  dri.xm = avx;
	  dri.ym = avy;
	  dri.zm = avz;
	  dri.sx = slpx;
	  dri.sy = slpy;
	  dri.lk = 0;
	  dri.ttdmin = 99.9; // isolation [mm]
	  driplets.push_back(dri);

	} // cl B

      } // cl C

    } // cl A

    // event selection:
    //if( driplets.size() > 4 ) continue;
    //if( driplets.size() < 6 ) continue;

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // make triplets 2+0-1:

    vector <triplet> triplets;

    double triCut = 0.05; // [mm]

    double dz13 = zz[3] - zz[1]; // from 1 to 3 in z

    for( vector<cluster>::iterator cA = cl[1].begin(); cA != cl[1].end(); ++cA ) {

      double xA = cA->col*ptchx[1] - alignx[1];
      double yA = cA->row*ptchy[1] - aligny[1];
      double xmid = xA - midx[1];
      double ymid = yA - midy[1];
      xA = xmid - ymid*rotx[1];
      yA = ymid + xmid*roty[1];

      for( vector<cluster>::iterator cC = cl[3].begin(); cC != cl[3].end(); ++cC ) {

	double xC = cC->col*ptchx[3] - alignx[3];
	double yC = cC->row*ptchy[3] - aligny[3];
	double xmid = xC - midx[3];
	double ymid = yC - midy[3];
	xC = xmid - ymid*rotx[3];
	yC = ymid + xmid*roty[3];

	double avx = 0.5 * ( xA + xC ); // mid
	double avy = 0.5 * ( yA + yC );
	double avz = 0.5 * ( zz[1] + zz[3] ); // mid z
 
	double slpx = ( xC - xA ) / dz13; // slope x
	double slpy = ( yC - yA ) / dz13; // slope y

	// middle plane B = 2:

	for( vector<cluster>::iterator cB = cl[2].begin(); cB != cl[2].end(); ++cB ) {

	  double xB = cB->col*ptchx[2] - alignx[2];
	  double yB = cB->row*ptchy[2] - aligny[2];
	  double xmid = xB - midx[2];
	  double ymid = yB - midy[2];
	  xB = xmid - ymid*rotx[2];
	  yB = ymid + xmid*roty[2];

	  // interpolate track to B:

	  double dz = zz[2] - avz;
	  double xk = avx + slpx * dz; // triplet at k
	  double yk = avy + slpy * dz;

	  double dx3 = xB - xk;
	  double dy3 = yB - yk;

	  // telescope triplet cuts:

	  if( fabs(dx3) > triCut ) continue;
	  if( fabs(dy3) > triCut ) continue;

	  triplet tri;
	  tri.xm = avx;
	  tri.ym = avy;
	  tri.zm = avz;
	  tri.sx = slpx;
	  tri.sy = slpy;
	  tri.lk = 0;
	  tri.ttdmin = 99.9; // isolation [mm]
	  triplets.push_back(tri);

	} // cl B

      } // cl C

    } // cl A

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // book event histos:

    ++nevd;

    cout << "display " << nevd
	 << " with " << ncl << " clusters, "
	 << triplets.size() << " triplets, "
	 << driplets.size() << " driplets"
	 << endl;

    double dz = zz[6] - zz[1];

    TH2D xzview( Form( "ev%ixz", nevd ),
		 Form( "x-z display %i trigger %i;z [mm];x [mm];hit", nevd, (int) iev ),
		 120, -0.1*dz, 1.1*dz, 220, -11, 11 );

    c1->cd();
    xzview.Draw("");

    TH2D yzview( Form( "ev%iyz", nevd ),
		 Form( "y-z display %i trigger %i;z [mm];y [mm];hit", nevd, (int) iev ),
		 120, -0.1*dz, 1.1*dz, 220, -11, 11 );

    c2->cd();
    yzview.Draw("");

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // triplets:

    vector <TLine> xlines; // in main scope
    vector <TLine> ylines; // in main scope

    double dz12 = zz[2] - zz[1];
    double z1 = zz[1] - 0.5 * dz12;
    double z3 = zz[3] + 0.5 * dz12;

    for( unsigned int iA = 0; iA < triplets.size(); ++iA ) { // iA = upstream

      double avx = triplets[iA].xm;
      double avy = triplets[iA].ym;
      double avz = triplets[iA].zm;
      double slx = triplets[iA].sx;
      double sly = triplets[iA].sy;

      double d1 = z1 - avz;
      double x1 = avx + slx * d1;
      double y1 = avy + sly * d1;

      double d3 = z3 - avz;
      double x3 = avx + slx * d3;
      double y3 = avy + sly * d3;

      TLine lx( z1, x1, z3, x3 );
      lx.SetLineColor(1);
      xlines.push_back(lx);

      TLine ly( z1, y1, z3, y3 );
      ly.SetLineColor(1);
      ylines.push_back(ly);

    } // loop triplets iA

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // driplets:

    double dz56 = zz[6] - zz[5];
    double z4 = zz[4] - 0.5 * dz56;
    double z6 = zz[6] + 0.5 * dz56;

    for( unsigned int jB = 0; jB < driplets.size(); ++jB ) { // jB = downstream

      double avx = driplets[jB].xm;
      double avy = driplets[jB].ym;
      double avz = driplets[jB].zm;
      double slx = driplets[jB].sx;
      double sly = driplets[jB].sy;

      double d4 = z4 - avz;
      double x4 = avx + slx * d4;
      double y4 = avy + sly * d4;

      double d6 = z6 - avz;
      double x6 = avx + slx * d6;
      double y6 = avy + sly * d6;

      TLine lx( z4, x4, z6, x6 );
      lx.SetLineColor(8);
      xlines.push_back(lx);

      TLine ly( z4, y4, z6, y6 );
      ly.SetLineColor(8);
      ylines.push_back(ly);

    } // loop driplets jB

    c1->cd();
    for( size_t ii = 0; ii < xlines.size(); ++ii )
      xlines[ii].Draw("same");

    c2->cd();
    for( size_t ii = 0; ii < ylines.size(); ++ii )
      ylines[ii].Draw("same");

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // Mimosa hits:

    vector <double> vx(ncl);
    vector <double> vy(ncl);
    vector <double> vz(ncl);

    int icl = 0;

    for( int ipl = 1; ipl <= 6; ++ipl ) {

      for( vector<cluster>::iterator cA = cl[ipl].begin(); cA != cl[ipl].end(); ++cA ) {

	double xA = cA->col*ptchx[ipl] - alignx[ipl];
	double yA = cA->row*ptchy[ipl] - aligny[ipl];
	double xmid = xA - midx[ipl];
	double ymid = yA - midy[ipl];
	xA = xmid - ymid*rotx[ipl];
	yA = ymid + xmid*roty[ipl];
	double zA = zz[ipl];
	vx.at(icl) = xA;
	vy.at(icl) = yA;
	vz.at(icl) = zA;
	++icl;

      } // cl

    } // planes

    c1->cd();

    TGraph gxz( icl, &vz[0], &vx[0] );
    gxz.SetMarkerColor(4);
    gxz.SetMarkerStyle(20);
    gxz.SetMarkerSize(0.8);
    gxz.Draw("P"); // without axis option: overlay

    c1->Update();

    c2->cd();

    TGraph gyz( icl, &vz[0], &vy[0] );
    gyz.SetMarkerColor(4);
    gyz.SetMarkerStyle(20);
    gyz.SetMarkerSize(0.8);
    gyz.Draw("P"); // without axis option: overlay

    c2->Update();

    cout << "enter any key, q to stop" << endl;

    while( !kbhit() )
      gSystem->ProcessEvents(); // ROOT

    string any;
    cin >> any;
    if( any == Q )
      more = 0;

  } while( more && reader->NextEvent() );

  delete reader;

  cout << "done after " << iev << " events" << endl;

  cout << endl;

  return 0;
}
